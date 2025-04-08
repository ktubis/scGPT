import os
import sys
import numpy as np
from torch.utils.data import Dataset, DataLoader
from torch import nn
from scipy.sparse import issparse
import torch
from typing import Dict
import json
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score

sys.path.insert(0, "../")
from scgpt.tokenizer import tokenize_and_pad_batch, random_mask_value
from scgpt.model import TransformerModel


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

    def __init__(self, pad_value, vocab, config_dict, num_batches, num_celltypes,
                 mask_value=MASK_VALUE, mask_ratio=MASK_RATIO):

        self.model = TransformerModel(ntoken=len(vocab), 
                                      num_batch_labels=num_batches,
                                      n_cls=num_celltypes, 
                                      vocab=vocab, 
                                      **config_dict)
        self.vocab = vocab
        self.mask_ratio = mask_ratio
        self.mask_value = mask_value
        self.pad_value = pad_value
        self.pad_token = config_dict["pad_token"]
        self.criterion = nn.CrossEntropyLoss()

    
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

    def train(self):
        pass
