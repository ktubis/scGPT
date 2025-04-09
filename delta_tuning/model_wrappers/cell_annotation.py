import os
import sys
import numpy as np
from torch.utils.data import Dataset, DataLoader
from torch import nn
from scipy.sparse import issparse
import torch
from typing import Dict, Tuple
import json
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
import time

sys.path.insert(0, "../")
from scgpt.tokenizer import tokenize_and_pad_batch, random_mask_value
from scgpt.model import TransformerModel
from sklearn.model_selection import train_test_split


INPUT_LAYER = "X_binned"
MASK_RATIO = 0.0
CLS = True
MASK_VALUE = -1



# TODO: move to outer file
class SeqDataset(Dataset):
    def __init__(self, data: Dict[str, torch.Tensor]):
        self.data = data

    def __len__(self):
        return self.data["gene_ids"].shape[0]

    def __getitem__(self, idx):
        return {k: v[idx] for k, v in self.data.items()}
    

class CellAnnotationModelWrapper():

    def __init__(self, model_path, pad_value, vocab, config_dict, num_batches, num_celltypes,
                 mask_value=MASK_VALUE, mask_ratio=MASK_RATIO):
        
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
                     train_celltype_labels, valid_celltype_labels) -> Tuple[Dict[str, torch.Tensor]]:

        genes = train_data.var.index.tolist()
        gene_ids = np.array(self.vocab(genes), dtype=int)

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
        input_values_train, input_values_valid = tokenized_train, tokenized_valid
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

    
    def _evaluate(self, loader: DataLoader, return_raw: bool = False) -> float:
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
                input_gene_ids = batch_data["gene_ids"]
                input_values = batch_data["values"]
                celltype_labels = batch_data["celltype_labels"]

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

        if return_raw:
            return np.concatenate(predictions, axis=0)

        return total_loss / total_num, total_error / total_num



    def test(self, adata: DataLoader, max_seq_len, eval_batch_size) -> float:
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
            max_len=max_seq_len,
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

    def train(self, adata: DataLoader, max_seq_len, train_batch_size):
        # train test split

        all_counts = (
            adata.layers[INPUT_LAYER].A
            if issparse(adata.layers[INPUT_LAYER])
            else adata.layers[INPUT_LAYER]
        )

        celltypes_labels = adata.obs["celltype_id"].tolist()  # make sure count from 0
        celltypes_labels = np.array(celltypes_labels)

        batch_ids = adata.obs["batch_id"].tolist()
        num_batch_types = len(set(batch_ids))
        batch_ids = np.array(batch_ids)

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

        train_data_pt, valid_data_pt = self._prepare_data_for_train(
            train_data, valid_data,
            train_batch_labels, valid_batch_labels,
            train_celltype_labels, valid_celltype_labels)

        """
        Train the model for one epoch.
        """
        self.model.train()
        (
            total_loss,
            total_mse,
            total_cls,
            total_zero_log_prob,
        ) = (0.0, 0.0, 0.0, 0.0)
        total_error = 0.0
        start_time = time.time()

        num_batches = len(loader)
        for batch, batch_data in enumerate(loader):
            input_gene_ids = batch_data["gene_ids"]
            input_values = batch_data["values"]
            target_values = batch_data["target_values"]
            batch_labels = batch_data["batch_labels"]
            celltype_labels = batch_data["celltype_labels"]

            src_key_padding_mask = input_gene_ids.eq(vocab[pad_token])
            with torch.cuda.amp.autocast(enabled=config.amp):
                output_dict = model(
                    input_gene_ids,
                    input_values,
                    src_key_padding_mask=src_key_padding_mask,
                    batch_labels=batch_labels if INPUT_BATCH_LABELS or config.DSBN else None,
                    CLS=True,
                    do_sample=do_sample_in_train,
                    #generative_training=False
                )

                loss = 0.0
                metrics_to_log = {}
                if CLS:
                    loss_cls = criterion_cls(output_dict["cls_output"], celltype_labels)
                    loss = loss + loss_cls
                    metrics_to_log.update({"train/cls": loss_cls.item()})

                    error_rate = 1 - (
                        (output_dict["cls_output"].argmax(1) == celltype_labels)
                        .sum()
                        .item()
                    ) / celltype_labels.size(0)

            self.model.zero_grad()
            scaler.scale(loss).backward()
            scaler.unscale_(optimizer)
            with warnings.catch_warnings(record=True) as w:
                warnings.filterwarnings("always")
                torch.nn.utils.clip_grad_norm_(
                    model.parameters(),
                    1.0,
                    error_if_nonfinite=False if scaler.is_enabled() else True,
                )
                if len(w) > 0:
                    logger.warning(
                        f"Found infinite gradient. This may be caused by the gradient "
                        f"scaler. The current scale is {scaler.get_scale()}. This warning "
                        "can be ignored if no longer occurs after autoscaling of the scaler."
                    )
            scaler.step(optimizer)
            scaler.update()

            #wandb.log(metrics_to_log)

            total_loss += loss.item()
            total_mse += loss_mse.item() if MLM else 0.0
            total_cls += loss_cls.item() if CLS else 0.0
            total_cce += loss_cce.item() if CCE else 0.0
            total_mvc += loss_mvc.item() if MVC else 0.0
            total_ecs += loss_ecs.item() if ECS else 0.0
            total_dab += loss_dab.item() if DAB else 0.0
            total_adv_E += loss_adv_E.item() if ADV else 0.0
            total_adv_D += loss_adv_D.item() if ADV else 0.0
            total_zero_log_prob += loss_zero_log_prob.item() if explicit_zero_prob else 0.0
            total_mvc_zero_log_prob += (
                loss_mvc_zero_log_prob.item() if MVC and explicit_zero_prob else 0.0
            )
            total_error += error_rate
            if batch % log_interval == 0 and batch > 0:
                lr = scheduler.get_last_lr()[0]
                ms_per_batch = (time.time() - start_time) * 1000 / log_interval
                cur_loss = total_loss / log_interval
                cur_mse = total_mse / log_interval
                cur_cls = total_cls / log_interval if CLS else 0.0
                cur_cce = total_cce / log_interval if CCE else 0.0
                cur_mvc = total_mvc / log_interval if MVC else 0.0
                cur_ecs = total_ecs / log_interval if ECS else 0.0
                cur_dab = total_dab / log_interval if DAB else 0.0
                cur_adv_E = total_adv_E / log_interval if ADV else 0.0
                cur_adv_D = total_adv_D / log_interval if ADV else 0.0
                cur_zero_log_prob = (
                    total_zero_log_prob / log_interval if explicit_zero_prob else 0.0
                )
                cur_mvc_zero_log_prob = (
                    total_mvc_zero_log_prob / log_interval
                    if MVC and explicit_zero_prob
                    else 0.0
                )
                cur_error = total_error / log_interval
                # ppl = math.exp(cur_loss)
                logger.info(
                    f"| epoch {epoch:3d} | {batch:3d}/{num_batches:3d} batches | "
                    f"lr {lr:05.4f} | ms/batch {ms_per_batch:5.2f} | "
                    f"loss {cur_loss:5.2f} | "
                    + (f"mse {cur_mse:5.2f} | mre {cur_error:5.2f} |" if MLM else "")
                    + (f"cls {cur_cls:5.2f} | " if CLS else "")
                    + (f"err {cur_error:5.2f} | " if CLS else "")
                    + (f"cce {cur_cce:5.2f} |" if CCE else "")
                    + (f"mvc {cur_mvc:5.2f} |" if MVC else "")
                    + (f"ecs {cur_ecs:5.2f} |" if ECS else "")
                    + (f"dab {cur_dab:5.2f} |" if DAB else "")
                    + (f"adv_E {cur_adv_E:5.2f} |" if ADV else "")
                    + (f"adv_D {cur_adv_D:5.2f} |" if ADV else "")
                    + (f"nzlp {cur_zero_log_prob:5.2f} |" if explicit_zero_prob else "")
                    + (
                        f"mvc_nzlp {cur_mvc_zero_log_prob:5.2f} |"
                        if MVC and explicit_zero_prob
                        else ""
                    )
                )
                total_loss = 0
                total_mse = 0
                total_cls = 0
                total_cce = 0
                total_mvc = 0
                total_ecs = 0
                total_dab = 0
                total_adv_E = 0
                total_adv_D = 0
                total_zero_log_prob = 0
                total_mvc_zero_log_prob = 0
                total_error = 0
                start_time = time.time()

    def run_training(self, num_epochs):
        best_val_loss = float("inf")
        best_avg_bio = 0.0
        best_model = None
        #define_wandb_metrcis()

        for epoch in range(1, num_epochs + 1):
            epoch_start_time = time.time()
            train_data_pt, valid_data_pt = prepare_data(sort_seq_batch=per_seq_batch_sample)
            train_loader = prepare_dataloader(
                train_data_pt,
                batch_size=batch_size,
                shuffle=False,
                intra_domain_shuffle=True,
                drop_last=False,
            )
            valid_loader = prepare_dataloader(
                valid_data_pt,
                batch_size=eval_batch_size,
                shuffle=False,
                intra_domain_shuffle=False,
                drop_last=False,
            )

            if config.do_train:
                train(
                    model,
                    loader=train_loader,
                )
            val_loss, val_err = evaluate(
                model,
                loader=valid_loader,
            )
            elapsed = time.time() - epoch_start_time
            logger.info("-" * 89)
            logger.info(
                f"| end of epoch {epoch:3d} | time: {elapsed:5.2f}s | "
                f"valid loss/mse {val_loss:5.4f} | err {val_err:5.4f}"
            )
            logger.info("-" * 89)

            if val_loss < best_val_loss:
                best_val_loss = val_loss
                best_model = copy.deepcopy(model)
                best_model_epoch = epoch
                logger.info(f"Best model with score {best_val_loss:5.4f}")

            scheduler.step()
            if DAB_separate_optim:
                scheduler_dab.step()
            if ADV:
                scheduler_D.step()
                scheduler_E.step()
        pass