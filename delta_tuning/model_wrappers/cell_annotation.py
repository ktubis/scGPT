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
from collections import namedtuple
from transformers import get_linear_schedule_with_warmup

from arcitecture.transformer_wrapper import copy_original_model

sys.path.insert(0, "../")
from scgpt.tokenizer import tokenize_and_pad_batch, random_mask_value
from scgpt.model import TransformerModel
import scgpt
import matplotlib.pyplot as plt


INPUT_LAYER = "X_binned"
MASK_RATIO = 0.0
CLS = True
MASK_VALUE = -1
LR = 1e-4
SCHEDULE_INTERVAL = 20
EPS = 1e-8
SCHEDULE_RATIO = 0.99
BATCH_SIZE = 32
EVAL_BATCH_SIZE = 64
LOG_INTERVAL = 100
RETRAINED_MODELS_DIR = "retrained_models/"
EARLY_STOPPING_EPOCHS_AVG = 10
EARLY_STOPPING_PATIENCE = 3
LR_FINDER_LOG_DIR = "cell_annotation_logs/lr_finder/"
FIND_LR_PERIOD = 600
FIND_LR_GAMMA = 3

TrainTestSplitResults = namedtuple(
    "TrainTestSplitResults",
    [
        "train_data",
        "valid_data",
        "train_celltype_labels",
        "valid_celltype_labels",
        "train_batch_labels",
        "valid_batch_labels",
    ],
)


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

def move_optimizer_params_to_cuda(optimizer):
        for group in optimizer.param_groups:
            for p in group['params']:
                if not p.is_cuda:
                    print(f"Found optimizer parameter on CPU: {p.shape}")
                    p.data = p.data.cuda()
                    if p.grad is not None:
                        p.grad.data = p.grad.data.cuda()
                        

class CellAnnotationModelWrapper():

    def __init__(self, model_path, pad_value, vocab, config_dict, num_batches, num_celltypes, max_seq_len, lr=LR, batch_size=BATCH_SIZE, eval_batch_size=BATCH_SIZE, log_dir="cell_annotation_logs/",
                 mask_value=MASK_VALUE, mask_ratio=MASK_RATIO, model_name="awesome_model", log_wandb=False, schedule_interval=SCHEDULE_INTERVAL, schedule_ratio=SCHEDULE_RATIO):
        
        self.model = TransformerModel(ntoken=len(vocab), 
                            num_batch_labels=num_batches,
                            n_cls=num_celltypes, 
                            vocab=vocab,
                            seq_len=max_seq_len,
                            train_batch_size=batch_size,
                            **config_dict)
                
        if model_path is not None:
            self.load_model(model_path)

        self.vocab = vocab
        self.mask_ratio = mask_ratio
        self.pad_token = config_dict["pad_token"]
        self.criterion = nn.CrossEntropyLoss()
        self.lr = lr
        self.schedule_interval = schedule_interval
        self.eps = EPS
        self.schedule_ratio = schedule_ratio
        self.batch_size = batch_size
        self.eval_batch_size = eval_batch_size
        self.log_interval = LOG_INTERVAL
        self.logger = scgpt.logger
        scgpt.utils.add_file_handler(self.logger, log_dir + f"{model_name}.log")
        self.max_seq_len = max_seq_len
        self.model_name = model_name
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.model.to(self.device)
        self.log_wandb = log_wandb

        if config_dict["input_emb_style"] == "category":
            self.mask_value = config_dict["n_bins"] + 1
            self.pad_value = config_dict["n_bins"]  # for padding gene expr values
        else:
            self.mask_value = mask_value
            self.pad_value = pad_value            

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
            copy_original_model(pretrained_dict, model_dict)
            pretrained_dict = {
                k: v
                for k, v in pretrained_dict.items()
                if k in model_dict and v.shape == model_dict[k].shape
            }
            model_dict.update(pretrained_dict)
            self.model.load_state_dict(model_dict)


    def _prepare_data_for_train(self, split_data: TrainTestSplitResults, gene_ids) -> Tuple[Dict[str, torch.Tensor]]:

        tokenized_train = tokenize_and_pad_batch(
            split_data.train_data,
            gene_ids,
            max_len=self.max_seq_len,
            vocab=self.vocab,
            pad_token=self.pad_token,
            pad_value=self.pad_value,
            append_cls=True,  # append <cls> token at the beginning
        )

        tokenized_valid = tokenize_and_pad_batch(
            split_data.valid_data,
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

        tensor_batch_labels_train = torch.from_numpy(split_data.train_batch_labels).long()
        tensor_batch_labels_valid = torch.from_numpy(split_data.valid_batch_labels).long()

        tensor_celltype_labels_train = torch.from_numpy(split_data.train_celltype_labels).long()
        tensor_celltype_labels_valid = torch.from_numpy(split_data.valid_celltype_labels).long()

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
    
    # Should be called every epoch, because it randomly selects the genes for each cell
    def _get_train_valid_data_per_epoch(self, split_data: TrainTestSplitResults, gene_ids):
        train_data_pt, valid_data_pt = self._prepare_data_for_train(split_data, gene_ids)

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
        return train_loader, valid_loader

    
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
                        input_ids=input_gene_ids,
                        inputs_embeds=input_values,
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
        total_loss_in_epoch = 0.0

        num_batches = len(data_loader)
        epoch_loss = 0
        for batch, batch_data in enumerate(data_loader):
            input_gene_ids = batch_data["gene_ids"].to(self.device)
            input_values = batch_data["values"].to(self.device)
            celltype_labels = batch_data["celltype_labels"].to(self.device)

            src_key_padding_mask = input_gene_ids.eq(self.vocab[self.pad_token])
            with torch.cuda.amp.autocast(enabled=True):
                output_dict = self.model(
                    input_ids=input_gene_ids,
                    inputs_embeds=input_values,
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

            if self.log_wandb:
                wandb.log({"train/loss": loss.item(), "train/err": error_rate})

            total_loss += loss.item()
            total_loss_in_epoch += loss.item()
            epoch_loss += loss.item()
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
                    f"lr {lr:05.8f} | ms/batch {ms_per_batch:5.2f} | "
                    f"loss {cur_loss:5.2f} | "
                    + (f"cls {cur_cls:5.2f} | " if CLS else "")
                    + (f"err {cur_error:5.2f} | " if CLS else "")
                )     
                total_loss = 0
                total_cls = 0
                total_error = 0
                start_time = time.time()
            scheduler.step()
        return total_loss_in_epoch / num_batches


    # Need to be called before training, once
    def _get_split_data(self, adata):
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

        split_data = TrainTestSplitResults(*train_test_split(
            all_counts, celltypes_labels, batch_ids, test_size=0.1, shuffle=True
        ))
        return split_data, gene_ids


    def train(self, num_epochs, adata, seed, adata_test, find_lr=False, warm_up_epochs=0, early_stop=False):

        split_data, gene_ids = self._get_split_data(adata)
        best_val_loss = float("inf")

        optimizer = torch.optim.Adam(
            self.model.parameters(), lr=self.lr, eps=self.eps
        )
        if find_lr:
            assert warm_up_epochs == 0, "Warm up epochs are not supported in lr finding mode."
            # If in the lr finding mode.
            scheduler = torch.optim.lr_scheduler.StepLR(
                optimizer, FIND_LR_PERIOD, gamma=FIND_LR_GAMMA
            )
        else:
            scheduler = get_linear_schedule_with_warmup(
                optimizer,
                num_warmup_steps=warm_up_epochs * np.ceil(len(split_data.train_data) / self.batch_size),
                num_training_steps=num_epochs * np.ceil(len(split_data.train_data) / self.batch_size),
            )
            #scheduler = torch.optim.lr_scheduler.StepLR(
            #    optimizer, self.schedule_interval, gamma=self.schedule_ratio
            #)
        move_optimizer_params_to_cuda(optimizer)
        scaler = torch.cuda.amp.GradScaler(enabled=True)

        # for early stopping
        last_test_losses = np.zeros(EARLY_STOPPING_EPOCHS_AVG)
        early_stopping_counter = 0
        prev_epoch_loss = np.inf
        curr_test_loss_avg = np.inf
        for epoch in range(1, num_epochs + 1):
            epoch_start_time = time.time()
            train_loader, valid_loader = self._get_train_valid_data_per_epoch(split_data, gene_ids)
            epoch_loss = self._train_step(train_loader, optimizer, scheduler, scaler, epoch)
            if find_lr:
                with open(LR_FINDER_LOG_DIR + self.model_name, 'a') as f:
                    f.write(f"{scheduler.get_last_lr()[0]} {epoch_loss}\n")
            val_loss, val_err = self._evaluate(loader=valid_loader, epoch=epoch)
            _, _, test_results = self.test(adata_test, eval_batch_size=self.eval_batch_size)

            self.logger.info(
                f"Test results: {test_results}"
            )

            elapsed = time.time() - epoch_start_time
            self.logger.info("-" * 89)
            self.logger.info(
                f"| end of epoch {epoch:3d} | time: {elapsed:5.2f}s | "
                f"valid loss/mse {val_loss:5.4f} | err {val_err:5.4f}"
            )
            self.logger.info("-" * 89)

            # Early stopping mechanism
            prev_test_loss_avg = curr_test_loss_avg
            last_test_losses[(epoch - 1) % EARLY_STOPPING_EPOCHS_AVG] = test_results["test/macro_f1"]
            curr_test_loss_avg = np.sum(last_test_losses) / np.minimum(EARLY_STOPPING_EPOCHS_AVG, epoch)
            if curr_test_loss_avg >= prev_test_loss_avg:
                early_stopping_counter += 1
            else:
                early_stopping_counter = 0
            if early_stop and early_stopping_counter >= EARLY_STOPPING_PATIENCE:
                # return the value for optuna hyperparameter search
                self.logger.info(
                    f"Early stopping at epoch {epoch} with score {curr_test_loss_avg:5.4f}"
                )
                return curr_test_loss_avg

            if find_lr and (epoch_loss > prev_epoch_loss):
                with open(LR_FINDER_LOG_DIR + self.model_name, 'r') as f:
                    lines = f.readlines()
                lrs = []
                losses = []
                for line in lines:
                    lr, loss = line.split()
                    lrs.append(lr)
                    losses.append(loss)
                plt.plot(lrs, losses)
                plt.savefig(LR_FINDER_LOG_DIR + self.model_name)
                return
            prev_epoch_loss = epoch_loss

            if val_loss < best_val_loss:
                best_val_loss = val_loss
                self.logger.info(f"Best model with score {best_val_loss:5.4f}")
                torch.save(self.model.state_dict(), RETRAINED_MODELS_DIR + self.model_name + '.pth')

            if self.log_wandb:
                wandb.log(
                    {
                        "train_loss": epoch_loss,
                        "best_val_loss": best_val_loss,
                        "val_loss": val_loss,
                        "test_loss": test_results["test/macro_f1"],
                        "epoch": epoch,
                    }
                )

            #scheduler.step()


    def finetune_cls_decoder(self):
        for param in self.model.parameters():
            param.requires_grad = False
        
        for param in self.model.cls_decoder.parameters():
            param.requires_grad = True

        trainable_params = sum(p.numel() for p in self.model.parameters() if p.requires_grad) / sum(p.numel() for p in self.model.parameters())
        print(f"Trainable parameters: {trainable_params}")

    def train_transformer_modules(self, train_modules: list):
        """
        Freeze all of the modules except the specified transformer modules.
        :param train_modules: list of module names to freeze
        """

        self.finetune_cls_decoder()  # make sure cls decoder is trainable

        modules_to_train = ["transformer_encoder.layers." + str(i) for i in train_modules]

        for name, module in self.model.named_modules():
            if name.startswith(tuple(modules_to_train)):
                for param in module.parameters():
                    param.requires_grad = True
                self.logger.info(f"Un-freezing module {name}")

        trainable_params = sum(p.numel() for p in self.model.parameters() if p.requires_grad) / sum(p.numel() for p in self.model.parameters())
        print(f"Trainable parameters after unfreezing transformers: {trainable_params}")
