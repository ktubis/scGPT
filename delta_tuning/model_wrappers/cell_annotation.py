import os
import sys
import numpy as np
from torch.utils.data import Dataset, DataLoader
from torch import nn
from scipy.sparse import issparse
import torch
from typing import Dict, Tuple
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
import time
import warnings
from sklearn.model_selection import train_test_split
import wandb

sys.path.insert(0, "../")
from scgpt.tokenizer import tokenize_and_pad_batch, random_mask_value
from scgpt.model import TransformerModel
import scgpt


INPUT_LAYER = "X_binned"
MASK_RATIO = 0.0
CLS = True
MASK_VALUE = -1
LR = 1e-4
SCHEDULE_INTERVAL = 100
EPS = 1e-8
SCHEDULE_RATIO = 0.9
BATCH_SIZE = 32
EVAL_BATCH_SIZE = 64
LOG_INTERVAL = 100
RETRAINED_MODELS_DIR = "retrained_models/"



# TODO: move to outer file
class SeqDataset(Dataset):
    def __init__(self, data: Dict[str, torch.Tensor]):
        self.data = data

    def __len__(self):
        return self.data["gene_ids"].shape[0]

    def __getitem__(self, idx):
        return {k: v[idx] for k, v in self.data.items()}
    

def prepare_dataloader(
    data_pt: Dict[str, torch.Tensor],
    batch_size: int,
    shuffle: bool = False,
    drop_last: bool = False,
    num_workers: int = 0,
) -> DataLoader:
    if num_workers == 0:
        num_workers = min(len(os.sched_getaffinity(0)), batch_size // 2)

    dataset = SeqDataset(data_pt)

    data_loader = DataLoader(
        dataset=dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        drop_last=drop_last,
        num_workers=num_workers,
        pin_memory=True,
    )
    return data_loader
    

def init_wandb(lr, model_name, epochs, batch_size, schedule_ratio, schedule_interval, seed):
    config = {
        "learning_rate": lr,
        "model_name": model_name,
        "epochs": epochs,
        "batch_size": batch_size,
        "schedule_ratio": schedule_ratio,
        "schedule_interval": schedule_interval,
        "seed": seed,
    }
    wandb.init(
        config=config,
        project="cell_annotation",
        reinit=True,
        settings=wandb.Settings(start_method="fork"))


class CellAnnotationModelWrapper():

    def __init__(self, model_path, pad_value, vocab, config_dict, num_batches, num_celltypes, max_seq_len, lr=LR, log_dir="cell_annotation_logs/",
                 mask_value=MASK_VALUE, mask_ratio=MASK_RATIO, model_name="awesome_model", wandb=False):
        
        self.model = TransformerModel(ntoken=len(vocab), 
                            num_batch_labels=num_batches,
                            n_cls=num_celltypes, 
                            vocab=vocab, 
                            **config_dict)
        
        if model_path is not None:
            self.load_model(model_path)

        self.vocab = vocab
        self.mask_ratio = mask_ratio
        self.mask_value = mask_value
        self.pad_value = pad_value
        self.pad_token = config_dict["pad_token"]
        self.criterion = nn.CrossEntropyLoss()
        self.lr = lr
        self.schedule_interval = SCHEDULE_INTERVAL
        self.eps = EPS
        self.schedule_ratio = SCHEDULE_RATIO
        self.batch_size = BATCH_SIZE
        self.eval_batch_size = EVAL_BATCH_SIZE
        self.log_interval = LOG_INTERVAL
        self.logger = scgpt.logger
        scgpt.utils.add_file_handler(self.logger, log_dir + "{model_name}.log")
        self.max_seq_len = max_seq_len
        self.model_name = model_name
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.model.to(self.device)
        self.model_name = model_name
        self.wandb = wandb


    def load_model(self, model_path):
        try:
            if not torch.cuda.is_available():
                self.model.load_state_dict(torch.load(model_path, map_location=torch.device('cpu')), strict=False)
            else:
                self.model.load_state_dict(torch.load(model_path))
        except:
            print("Error loading model, trying to load only matching parameters")
            # only load params that are in the model and match the size
            model_dict = self.model.state_dict()
            pretrained_dict = torch.load(model_path)
            pretrained_dict = {
                k: v
                for k, v in pretrained_dict.items()
                if k in model_dict and v.shape == model_dict[k].shape
            }
            model_dict.update(pretrained_dict)
            self.model.load_state_dict(model_dict)


    def _prepare_data_for_train(self, train_data, valid_data,
                     train_batch_labels, valid_batch_labels,
                     train_celltype_labels, valid_celltype_labels,
                     gene_ids) -> Tuple[Dict[str, torch.Tensor]]:
        
        print(train_data.shape)
        print(train_data)

        tokenized_train = tokenize_and_pad_batch(
            train_data,
            gene_ids,
            max_len=self.max_seq_len,
            vocab=self.vocab,
            pad_token=self.pad_token,
            pad_value=self.pad_value,
            append_cls=True,  # append <cls> token at the beginning
        )

        tokenized_valid = tokenize_and_pad_batch(
            valid_data,
            gene_ids,
            max_len=self.max_seq_len,
            vocab=self.vocab,
            pad_token=self.pad_token,
            pad_value=self.pad_value,
            append_cls=True,
        )

        input_gene_ids_train, input_gene_ids_valid = (
            tokenized_train["genes"],
            tokenized_valid["genes"],
        )

        masked_values_train = random_mask_value(
            tokenized_train["values"],
            mask_ratio=self.mask_ratio,
            mask_value=self.mask_value,
            pad_value=self.pad_value,
        )

        masked_values_valid = random_mask_value(
            tokenized_valid["values"],
            mask_ratio=self.mask_ratio,
            mask_value=self.mask_value,
            pad_value=self.pad_value,
        )

        input_values_train, input_values_valid = masked_values_train, masked_values_valid
        target_values_train, target_values_valid = (
            tokenized_train["values"],
            tokenized_valid["values"],
        )

        tensor_batch_labels_train = torch.from_numpy(train_batch_labels).long()
        tensor_batch_labels_valid = torch.from_numpy(valid_batch_labels).long()

        tensor_celltype_labels_train = torch.from_numpy(train_celltype_labels).long()
        tensor_celltype_labels_valid = torch.from_numpy(valid_celltype_labels).long()

        train_data_pt = {
            "gene_ids": input_gene_ids_train,
            "values": input_values_train,
            "target_values": target_values_train,
            "batch_labels": tensor_batch_labels_train,
            "celltype_labels": tensor_celltype_labels_train,
        }
        valid_data_pt = {
            "gene_ids": input_gene_ids_valid,
            "values": input_values_valid,
            "target_values": target_values_valid,
            "batch_labels": tensor_batch_labels_valid,
            "celltype_labels": tensor_celltype_labels_valid,
        }

        return train_data_pt, valid_data_pt

    
    def _evaluate(self, loader: DataLoader, epoch=0, return_raw: bool = False) -> float:
        """
        Evaluate the model on the evaluation data.
        """
        self.model.eval()
        total_loss = 0.0
        total_error = 0.0
        total_num = 0
        predictions = []
        with torch.no_grad():
            for batch_data in loader:
                input_gene_ids = batch_data["gene_ids"].to(self.device)
                input_values = batch_data["values"].to(self.device)
                celltype_labels = batch_data["celltype_labels"].to(self.device)

                src_key_padding_mask = input_gene_ids.eq(self.vocab[self.pad_token])
                with torch.cuda.amp.autocast(enabled=True):
                    output_dict = self.model(
                        input_gene_ids,
                        input_values,
                        src_key_padding_mask=src_key_padding_mask,
                        CLS=True,
                    )
                    output_values = output_dict["cls_output"]
                    loss = self.criterion(output_values, celltype_labels)

                total_loss += loss.item() * len(input_gene_ids)
                accuracy = (output_values.argmax(1) == celltype_labels).sum().item()
                total_error += (1 - accuracy / len(input_gene_ids)) * len(input_gene_ids)
                total_num += len(input_gene_ids)
                preds = output_values.argmax(1).cpu().numpy()
                predictions.append(preds)

        if self.wandb:
            wandb.log(
                {
                    "valid/mse": total_loss / total_num,
                    "valid/err": total_error / total_num,
                    "epoch": epoch,
                }
            )

        if return_raw:
            return np.concatenate(predictions, axis=0)

        return total_loss / total_num, total_error / total_num



    def test(self, adata: DataLoader, eval_batch_size=EVAL_BATCH_SIZE) -> float:
        all_counts = (
            adata.layers[INPUT_LAYER].A
            if issparse(adata.layers[INPUT_LAYER])
            else adata.layers[INPUT_LAYER]
        )

        celltypes_labels = adata.obs["celltype_id"].tolist()  # make sure count from 0
        celltypes_labels = np.array(celltypes_labels)

        batch_ids = adata.obs["batch_id"].tolist()
        batch_ids = np.array([int(batch_id) for batch_id in batch_ids])

        genes = adata.var.index.tolist()
        gene_ids = np.array(self.vocab(genes), dtype=int)

        tokenized_test = tokenize_and_pad_batch(
            all_counts,
            gene_ids,
            max_len=self.max_seq_len,
            vocab=self.vocab,
            pad_token=self.pad_token,
            pad_value=self.pad_value,
            append_cls=True,  # append <cls> token at the beginning
        )

        input_values_test = random_mask_value(
            tokenized_test["values"],
            mask_ratio=self.mask_ratio,
            mask_value=self.mask_value,
            pad_value=self.pad_value,
        )

        test_data_pt = {
            "gene_ids": tokenized_test["genes"],
            "values": input_values_test,
            "target_values": tokenized_test["values"],
            "batch_labels": torch.from_numpy(batch_ids).long(),
            "celltype_labels": torch.from_numpy(celltypes_labels).long(),
        }

        test_loader = DataLoader(
            dataset=SeqDataset(test_data_pt),
            batch_size=eval_batch_size,
            shuffle=False,
            drop_last=False,
            num_workers=min(len(os.sched_getaffinity(0)), eval_batch_size // 2),
            pin_memory=True,
        )

        predictions = self._evaluate(
            loader=test_loader,
            return_raw=True,
        )

        # compute accuracy, precision, recall, f1
        accuracy = accuracy_score(celltypes_labels, predictions)
        precision = precision_score(celltypes_labels, predictions, average="macro")
        recall = recall_score(celltypes_labels, predictions, average="macro")
        macro_f1 = f1_score(celltypes_labels, predictions, average="macro")

        results = {
            "test/accuracy": accuracy,
            "test/precision": precision,
            "test/recall": recall,
            "test/macro_f1": macro_f1,
        }

        return predictions, celltypes_labels, results

    def _train_step(self, data_loader: DataLoader, optimizer, scheduler, scaler, epoch):

        self.model.train()
        total_loss = 0.0
        total_cls = 0.0
        total_error = 0.0
        start_time = time.time()

        num_batches = len(data_loader)
        for batch, batch_data in enumerate(data_loader):
            input_gene_ids = batch_data["gene_ids"].to(self.device)
            input_values = batch_data["values"].to(self.device)
            celltype_labels = batch_data["celltype_labels"].to(self.device)

            src_key_padding_mask = input_gene_ids.eq(self.vocab[self.pad_token])
            with torch.cuda.amp.autocast(enabled=True):
                output_dict = self.model(
                    input_gene_ids,
                    input_values,
                    src_key_padding_mask=src_key_padding_mask,
                    batch_labels=None,
                    CLS=True,
                    do_sample=False,
                )

                loss = 0.0
                metrics_to_log = {}
                loss_cls = self.criterion(output_dict["cls_output"], celltype_labels)
                loss = loss + loss_cls
                metrics_to_log.update({"train/cls": loss_cls.item()})

                error_rate = 1 - (
                    (output_dict["cls_output"].argmax(1) == celltype_labels)
                    .sum()
                    .item()
                ) / celltype_labels.size(0)

            self.model.zero_grad()
            if torch.cuda.is_available():
                loss = loss.cuda()
            scaler.scale(loss).backward()
            scaler.unscale_(optimizer)
            with warnings.catch_warnings(record=True) as w:
                warnings.filterwarnings("always")
                torch.nn.utils.clip_grad_norm_(
                    self.model.parameters(),
                    1.0,
                    error_if_nonfinite=False if scaler.is_enabled() else True,
                )
                if len(w) > 0:
                    self.logger.warning(
                        f"Found infinite gradient. This may be caused by the gradient "
                        f"scaler. The current scale is {scaler.get_scale()}. This warning "
                        "can be ignored if no longer occurs after autoscaling of the scaler."
                    )
            scaler.step(optimizer)
            scaler.update()

            if self.wandb:
                wandb.log({"train/loss": cur_loss, "train/cls": cur_cls, "train/err": cur_error})

            total_loss += loss.item()
            total_cls += loss_cls.item() if CLS else 0.0
            total_error += error_rate
            if batch % self.log_interval == 0 and batch > 0:
                lr = scheduler.get_last_lr()[0]
                ms_per_batch = (time.time() - start_time) * 1000 / self.log_interval
                cur_loss = total_loss / self.log_interval
                cur_cls = total_cls / self.log_interval
                cur_error = total_error / self.log_interval
                self.logger.info(
                    f"| epoch {epoch:3d} | {batch:3d}/{num_batches:3d} batches | "
                    f"lr {lr:05.4f} | ms/batch {ms_per_batch:5.2f} | "
                    f"loss {cur_loss:5.2f} | "
                    + (f"cls {cur_cls:5.2f} | " if CLS else "")
                    + (f"err {cur_error:5.2f} | " if CLS else "")
                )
                total_loss = 0
                total_cls = 0
                total_error = 0
                start_time = time.time()


    def train(self, num_epochs, adata, seed, adata_test1, adata_test2):

        if self.wandb:
            init_wandb(self.lr, self.model_name, num_epochs, self.batch_size, self.schedule_ratio, self.schedule_interval, seed)
            wandb.watch(self.model)

        all_counts = (
            adata.layers[INPUT_LAYER].A
            if issparse(adata.layers[INPUT_LAYER])
            else adata.layers[INPUT_LAYER]
        )

        celltypes_labels = adata.obs["celltype_id"].tolist()  # make sure count from 0
        celltypes_labels = np.array(celltypes_labels)

        batch_ids = adata.obs["batch_id"].tolist()
        batch_ids = np.array(batch_ids)

        genes = adata.var.index.tolist()
        gene_ids = np.array(self.vocab(genes), dtype=int)

        (
            train_data,
            valid_data,
            train_celltype_labels,
            valid_celltype_labels,
            train_batch_labels,
            valid_batch_labels,
        ) = train_test_split(
            all_counts, celltypes_labels, batch_ids, test_size=0.1, shuffle=True
        )

        best_val_loss = float("inf")
        #define_wandb_metrcis()

        optimizer = torch.optim.Adam(
            self.model.parameters(), lr=self.lr, eps=self.eps
        )
        scheduler = torch.optim.lr_scheduler.StepLR(
            optimizer, self.schedule_interval, gamma=self.schedule_ratio
        )

        scaler = torch.cuda.amp.GradScaler(enabled=True)

        for group in optimizer.param_groups:
            for p in group['params']:
                if not p.is_cuda:
                    print(f"Found optimizer parameter on CPU: {p.shape}")
                    #p.data = p.data.cuda()
                    if p.grad is not None:
                        p.grad.data = p.grad.data.cuda()

        for epoch in range(1, num_epochs + 1):
            epoch_start_time = time.time()
            train_data_pt, valid_data_pt = self._prepare_data_for_train(
                train_data, valid_data,
                train_batch_labels, valid_batch_labels,
                train_celltype_labels, valid_celltype_labels,
                gene_ids
            )

            train_loader = prepare_dataloader(
                train_data_pt,
                batch_size=self.batch_size,
                shuffle=False,
                drop_last=False,
            )
            valid_loader = prepare_dataloader(
                valid_data_pt,
                batch_size=self.eval_batch_size,
                shuffle=False,
                drop_last=False,
            )

            self._train_step(train_loader, optimizer, scheduler, scaler, epoch)
            val_loss, val_err = self._evaluate(loader=valid_loader, epoch=epoch)
            _, _, test_results1 = self.test(adata_test1, eval_batch_size=self.eval_batch_size)
            _, _, test_results2 = self.test(adata_test2, eval_batch_size=self.eval_batch_size)
            self.logger.info(
                f"Test results on Muraro: {test_results1}"
            )
            self.logger.info(
                f"Test results on Xin: {test_results2}"
            )

            elapsed = time.time() - epoch_start_time
            self.logger.info("-" * 89)
            self.logger.info(
                f"| end of epoch {epoch:3d} | time: {elapsed:5.2f}s | "
                f"valid loss/mse {val_loss:5.4f} | err {val_err:5.4f}"
            )
            self.logger.info("-" * 89)

            if val_loss < best_val_loss:
                best_val_loss = val_loss
                self.logger.info(f"Best model with score {best_val_loss:5.4f}")
                torch.save(self.model.state_dict(), RETRAINED_MODELS_DIR + self.model_name + '.pth')

            scheduler.step()

        #TODO: save the model. Also add intermediate saving steps.

    def finetune_cls_decoder(self):
        for param in self.model.parameters():
            param.requires_grad = False
        
        for param in self.model.cls_decoder.parameters():
            param.requires_grad = True
