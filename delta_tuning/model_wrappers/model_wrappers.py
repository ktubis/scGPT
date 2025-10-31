import os
import sys
import numpy as np
import pandas as pd
from torch.utils.data import Dataset, DataLoader
from torch import nn
from scipy.sparse import issparse
import torch
from typing import Dict, Tuple
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
import time
import warnings
from sklearn.model_selection import train_test_split
#import wandb
from collections import namedtuple
from transformers import get_linear_schedule_with_warmup
import anndata as ad
from abc import ABC, abstractmethod
import scanpy as sc
import traceback
import matplotlib.pyplot as plt
from enum import Enum
import json
from copy import deepcopy
from gears import PertData
import load_ds

from arcitecture.transformer_wrapper import copy_original_model

sys.path.insert(0, "../")
from scgpt.tokenizer import tokenize_and_pad_batch, random_mask_value
from scgpt.model import TransformerModel
import scgpt
from scgpt import SubsetsBatchSampler
from scgpt.utils import eval_scib_metrics
from scgpt.preprocess import Preprocessor
from scgpt.loss import (
    masked_mse_loss,
    masked_relative_error,
    criterion_neg_log_bernoulli,
)
from scgpt.utils import map_raw_id_to_vocab_id, compute_perturbation_metrics

INPUT_LAYER = "X_binned"
LOG_INTERVAL = 100
RETRAINED_MODELS_DIR = "retrained_models/"
EARLY_STOPPING_EPOCHS_AVG = 5
EARLY_STOPPING_PATIENCE = 5
LR_FINDER_LOG_DIR = "cell_annotation_logs/lr_finder/"
FIND_LR_PERIOD = 600
FIND_LR_GAMMA = 3
TASKS_MODELS_CONFIG_DIR = "tasks_model_configs/"
MODEL_CONFIG_DIR = "model_configs/"
PREDICTIONS_DIR = "predictions/"
INTERMEDIATE_EMBEDDINGS_DIR = "intermediate_embeddings/"
CELL_EMBEDDINGS_DIR = "cell_embeddings/"

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

def add_tokens_to_vocab(pad_token, vocab):
        special_tokens = [pad_token, "<cls>", "<eoc>"]
        for s in special_tokens:
            if s not in vocab:
                vocab.append_token(s)


# TODO: move to outer file
class SeqDataset(Dataset):
    def __init__(self, data: Dict[str, torch.Tensor]):
        self.data = data

    def __len__(self):
        return self.data["gene_ids"].shape[0]

    def __getitem__(self, idx):
        return {k: v[idx] for k, v in self.data.items()}


class ModelDataloader():

    def __init__(self, data_name, vocab, pad_token, pad_value, batch_size, max_seq_len,
                 mask_ratio, mask_value):
        self.vocab = vocab
        self.data_name = data_name
        self.adata_train = None
        self.adata_test = None
        add_tokens_to_vocab(pad_token, self.vocab)
        self.vocab.set_default_index(vocab[pad_token])
        self.batch_size = batch_size
        self.pad_token = pad_token
        self.pad_value = pad_value
        self.max_seq_len = max_seq_len
        self.mask_ratio = mask_ratio
        self.mask_value = mask_value

    @abstractmethod
    def preprocess(self):
        pass

    def get_dataloader(self):
        pass

    @abstractmethod
    def get_num_training_steps(self):
        pass




# Used by celltype annotation and batch correction
class SimpleDataloader(ModelDataloader):

    def __init__(self, data_name, vocab, n_input_bins, n_hvg, pad_token,
                 pad_value, batch_size, max_seq_len, model_task, mask_ratio,
                 mask_value):
        super().__init__(data_name, vocab, pad_token, pad_value, batch_size, max_seq_len,
                         mask_ratio, mask_value)
        self.ds_loader = load_ds.get_data_loader(self.data_name)
        self._process_train_test(n_input_bins, n_hvg, model_task)

    def get_num_celltypes(self):
        return self.ds_loader.get_num_celltypes()
    
    def get_num_batches(self):
        return self.ds_loader.get_num_batches()

    def _preprocess_data(self, adata, preprocessor, vocab, batch_key=None):
        """
        Filter genes in the AnnData object based on their presence in the vocabulary.
        
        Args:
            adata (AnnData): The AnnData object containing gene expression data.
        
        Returns:
            adata (AnnData): The filtered AnnData object.
        """
        preprocessor(adata, batch_key=batch_key)
        adata.var["id_in_vocab"] = [
            1 if gene in vocab else -1 for gene in adata.var.index
        ]
        adata = adata[:, adata.var["id_in_vocab"] >= 0]


    def _process_train_test(self, n_input_bins, n_hvg, model_task):
        self.adata_train, self.adata_test = self.ds_loader.get_train_test()

        # set up the preprocessor, use the args to config the workflow
        preprocessor = Preprocessor(
            use_key="X",
            normalize_total=0.0,
            binning=n_input_bins,
            result_binned_key=INPUT_LAYER,
            subset_hvg=n_hvg if n_hvg < self.adata_train.shape[1] else False,
            hvg_flavor="seurat"    # Assumes data is not raw
        )

        batch_key = "str_batch" if model_task == "batch_correction" else None
        self._preprocess_data(self.adata_train, preprocessor, self.vocab, batch_key)
        self._preprocess_data(self.adata_test, preprocessor, self.vocab, batch_key)

    # Needs to be called before training, once
    def _get_split_data(self):
        all_counts = (
            self.adata_train.layers[INPUT_LAYER].A
            if issparse(self.adata_train.layers[INPUT_LAYER])
            else self.adata_train.layers[INPUT_LAYER]
        )

        celltypes_labels = self.adata_train.obs["celltype_id"].tolist()  # make sure count from 0
        celltypes_labels = np.array(celltypes_labels)

        batch_ids = self.adata_train.obs["batch_id"].tolist()
        batch_ids = np.array(batch_ids)

        genes = self.adata_train.var.index.tolist()
        gene_ids = np.array(self.vocab(genes), dtype=int)

        split_data = TrainTestSplitResults(*train_test_split(
            all_counts, celltypes_labels, batch_ids, test_size=0.1, shuffle=True
        ))
        return split_data, gene_ids

    
    def _prepare_data_for_train(self, sort_seq_batch=False) -> Tuple[Dict[str, torch.Tensor]]:

        input_gene_ids_train, input_gene_ids_valid = (
            self.tokenized_train["genes"],
            self.tokenized_valid["genes"],
        )

        masked_values_train = random_mask_value(
            self.tokenized_train["values"],
            mask_ratio=self.mask_ratio,
            mask_value=self.mask_value,
            pad_value=self.pad_value,
        )

        masked_values_valid = random_mask_value(
            self.tokenized_valid["values"],
            mask_ratio=self.mask_ratio,
            mask_value=self.mask_value,
            pad_value=self.pad_value,
        )

        input_values_train, input_values_valid = masked_values_train, masked_values_valid
        target_values_train, target_values_valid = (
            self.tokenized_train["values"],
            self.tokenized_valid["values"],
        )

        tensor_batch_labels_train = torch.from_numpy(self.split_data.train_batch_labels).long()
        tensor_batch_labels_valid = torch.from_numpy(self.split_data.valid_batch_labels).long()

        tensor_celltype_labels_train = torch.from_numpy(self.split_data.train_celltype_labels).long()
        tensor_celltype_labels_valid = torch.from_numpy(self.split_data.valid_celltype_labels).long()

        if sort_seq_batch:
            train_sort_ids = np.argsort(self.split_data.train_batch_labels)
            input_gene_ids_train = input_gene_ids_train[train_sort_ids]
            input_values_train = input_values_train[train_sort_ids]
            target_values_train = target_values_train[train_sort_ids]
            tensor_batch_labels_train = tensor_batch_labels_train[train_sort_ids]

            valid_sort_ids = np.argsort(self.split_data.valid_batch_labels)
            input_gene_ids_valid = input_gene_ids_valid[valid_sort_ids]
            input_values_valid = input_values_valid[valid_sort_ids]
            target_values_valid = target_values_valid[valid_sort_ids]
            tensor_batch_labels_valid = tensor_batch_labels_valid[valid_sort_ids]


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
    
    def get_num_training_steps(self):
        return np.ceil(len(self.split_data.train_data) / self.batch_size)

    def prepare_train_and_valid(self, include_zero_gene):
        self.split_data, gene_ids = self._get_split_data()
        self.tokenized_train = tokenize_and_pad_batch(
            self.split_data.train_data,
            gene_ids,
            max_len=self.max_seq_len,
            vocab=self.vocab,
            pad_token=self.pad_token,
            pad_value=self.pad_value,
            append_cls=True,  # append <cls> token at the beginning
            include_zero_gene=include_zero_gene
        )
    
        self.tokenized_valid = tokenize_and_pad_batch(
            self.split_data.valid_data,
            gene_ids,
            max_len=self.max_seq_len,
            vocab=self.vocab,
            pad_token=self.pad_token,
            pad_value=self.pad_value,
            append_cls=True,  # append <cls> token at the beginning
            include_zero_gene=include_zero_gene
        )
        
    def _prepare_dataloader(
            self,
            data_pt: Dict[str, torch.Tensor],
            shuffle: bool = False,
            drop_last: bool = False,
            num_workers: int = 0,
            intra_domain_shuffle: bool = False,
            per_seq_batch_sample: bool = False,
    ) -> DataLoader:
        if num_workers == 0:
            num_workers = min(len(os.sched_getaffinity(0)), self.batch_size // 2)

        dataset = SeqDataset(data_pt)

        if per_seq_batch_sample:
            # find the indices of samples in each seq batch
            subsets = []
            batch_labels_array = data_pt["batch_labels"].numpy()
            for batch_label in np.unique(batch_labels_array):
                batch_indices = np.where(batch_labels_array == batch_label)[0].tolist()
                subsets.append(batch_indices)
            data_loader = DataLoader(
                dataset=dataset,
                batch_sampler=SubsetsBatchSampler(
                    subsets,
                    self.batch_size,
                    intra_subset_shuffle=intra_domain_shuffle,
                    inter_subset_shuffle=shuffle,
                    drop_last=drop_last,
                ),
                num_workers=num_workers,
                pin_memory=True,
            )
            return data_loader

        data_loader = DataLoader(
            dataset=dataset,
            batch_size=self.batch_size,
            shuffle=shuffle,
            drop_last=drop_last,
            num_workers=num_workers,
            pin_memory=True,
        )
        return data_loader
    
    # Should be called every epoch, because it randomly selects the genes for each cell
    def get_train_valid_data_per_epoch(self,
                                       sort_batches,
                                       train_dataloader_params_dict,
                                       test_dataloader_params_dict):
        # mask_ratio, mask_value, pad_value, sort_seq_batch=False
        train_data_pt, valid_data_pt = self._prepare_data_for_train(sort_batches)

        train_loader = self._prepare_dataloader(
            train_data_pt,
            **train_dataloader_params_dict 
        )

        valid_loader = self._prepare_dataloader(
            valid_data_pt,
            **test_dataloader_params_dict 
        )
        return train_loader, valid_loader
    

    def get_test_loader(self, include_zero_gene, celltype_list = []):
        if celltype_list:
            self.adata_test = self.adata_test[self.adata_test.obs["celltype"].isin(celltype_list)]
        all_counts = (
            self.adata_test.layers[INPUT_LAYER].A
            if issparse(self.adata_test.layers[INPUT_LAYER])
            else self.adata_test.layers[INPUT_LAYER]
        )

        celltypes_labels = self.adata_test.obs["celltype_id"].tolist()
        celltypes_labels = np.array(celltypes_labels)

        batch_ids = self.adata_test.obs["batch_id"].tolist()
        batch_ids = np.array(batch_ids)

        genes = self.adata_test.var.index.tolist()
        gene_ids = np.array(self.vocab(genes), dtype=int)

        # Evaluate cls cell embeddings
        tokenized = tokenize_and_pad_batch(
            all_counts,
            gene_ids,
            max_len=self.max_seq_len,
            vocab=self.vocab,
            pad_token=self.pad_token,
            pad_value=self.pad_value,
            append_cls=True,  # append <cls> token at the beginning
            include_zero_gene=include_zero_gene
        )

    
        input_values_test = random_mask_value(
            tokenized["values"],
            mask_ratio=self.mask_ratio,
            mask_value=self.mask_value,
            pad_value=self.pad_value,
        )

        test_data_pt = {
            "gene_ids": tokenized["genes"],
            "values": input_values_test,
            "target_values": tokenized["values"],
            "batch_labels": torch.from_numpy(batch_ids).long(),
            "celltype_labels": torch.from_numpy(celltypes_labels).long(),
        }

        test_loader = DataLoader(
            dataset=SeqDataset(test_data_pt),
            batch_size=self.batch_size,
            shuffle=False,
            drop_last=False,
            num_workers=min(len(os.sched_getaffinity(0)), self.batch_size // 2),
            pin_memory=True,
        )

        return test_loader

        

class PertDataloader(ModelDataloader):

    #TODO: get gene_ids and stuff. See how it's done in the other classes.
    def __init__(self, data_name, vocab, n_input_bins, n_hvg, pad_token, pad_value, batch_size, max_seq_len,
                         mask_ratio, mask_value):
        super().__init__(data_name, vocab, pad_token, pad_value, batch_size, max_seq_len,
                         mask_ratio, mask_value)
        self.ds_loader = load_ds.get_data_loader(self.data_name)
        data_path = self.ds_loader.get_data_path()
        self.pert_data = PertData("./data")
        self.pert_data.load(data_path=data_path)
        self.pert_data.prepare_split(split="simulation", seed=1)
        self.pert_data.get_dataloader(batch_size=batch_size, test_batch_size=batch_size)
        self.ds_loader.load_train_test(self.pert_data)
    
    def get_train_valid_data_per_epoch(self,
                                       sort_batches,
                                       train_dataloader_params_dict,
                                       test_dataloader_params_dict):
        return (
            self.pert_data.dataloader["train_loader"],
            self.pert_data.dataloader["val_loader"]
            )
    
    def get_test_loader(self, include_zero_gene):
        return self.pert_data.dataloader["test_loader"]
    
    def prepare_train_and_valid(self, include_zero_gene):
        return
    
    def get_num_genes(self):
        n_genes = len(self.pert_data.adata.var["gene_name"].tolist())
        return n_genes
    
    def get_num_celltypes(self):
        return self.ds_loader.get_num_celltypes()
    
    def get_num_batches(self):
        return self.ds_loader.get_num_batches()
    
    def get_gene_ids(self):
        genes = self.pert_data.adata.var["gene_name"].tolist()
        gene_ids = np.array(
            [self.vocab[gene] if gene in self.vocab else self.vocab[self.pad_token] for gene in genes],
            dtype=int
        )
        return gene_ids
    
    def get_ctrl_condition(self):
        return self.pert_data.adata[self.pert_data.adata.obs["condition"] == "ctrl"]
    
    def get_num_training_steps(self):
        return len(self.pert_data.dataloader["train_loader"])

def move_optimizer_params_to_cuda(optimizer):
        for group in optimizer.param_groups:
            for p in group['params']:
                if not p.is_cuda:
                    print(f"Found optimizer parameter on CPU: {p.shape}")
                    p.data = p.data.cuda()
                    if p.grad is not None:
                        p.grad.data = p.grad.data.cuda()

class ModelLoader():
    class SupportedTasks(Enum):
        CELLTYPE_ANNOTATION = "celltype_annotation"
        BATCH_CORRECTION = "batch_correction"
        GE_PREDICTION = "ge_prediction"
        PERTURBATION = "perturbation"

    def __init__(self, model_task, small_model=False):
        supported = [task.value for task in ModelLoader.SupportedTasks]
        if model_task not in supported:
            raise ValueError(f"Unsupported task '{model_task}'. \
                             Supported downstream tasks are: {supported}")
        self.model_task = model_task
        with open(f"{TASKS_MODELS_CONFIG_DIR}{model_task}.json") as f:
            self.task_train_dict = json.load(f)
        if small_model:
            model_config_file = f"{MODEL_CONFIG_DIR}cpu_model.json"
        else:
            model_config_file = f"{MODEL_CONFIG_DIR}{model_task}.json"
        with open(model_config_file) as f:
            self.model_dict = json.load(f)

    def get_hvg(self):
        if self.model_task == ModelLoader.SupportedTasks.BATCH_CORRECTION.value:
            return self.model_dict['seq_len'] - 1
        else:
            return self.model_dict['seq_len']
    
    def get_input_bins(self):
        return self.task_train_dict['n_input_bins']
    
    def get_pad_token(self):
        return self.task_train_dict['pad_token']

    def get_model(self, model_config_dict):
        if self.model_task == ModelLoader.SupportedTasks.CELLTYPE_ANNOTATION.value:
            return CellAnnotation(task_training_dict=self.task_train_dict.copy(),
                                  config_dict=self.model_dict.copy(), **model_config_dict,
                                  n_hvg=self.get_hvg(), n_input_bins=self.task_train_dict['n_input_bins'])
        
        if self.model_task == ModelLoader.SupportedTasks.BATCH_CORRECTION.value:
            return BatchCorrection(task_training_dict=self.task_train_dict.copy(),
                                   config_dict=self.model_dict.copy(), **model_config_dict,
                                   n_hvg=self.get_hvg(), n_input_bins=self.task_train_dict['n_input_bins'])
        
        if self.model_task == ModelLoader.SupportedTasks.PERTURBATION.value:
            return Perturbation(task_training_dict=self.task_train_dict.copy(),
                                config_dict=self.model_dict.copy(), **model_config_dict,
                                n_hvg=self.get_hvg(), n_input_bins=self.task_train_dict['n_input_bins'])
        
        if self.model_task == ModelLoader.SupportedTasks.GE_PREDICTION.value:
            return GEPrediction(task_training_dict=self.task_train_dict.copy(),
                                config_dict=self.model_dict.copy(), **model_config_dict,
                                n_hvg=self.get_hvg(), n_input_bins=self.task_train_dict['n_input_bins'])
        
    
                        

class ScGPTModelWrapper(ABC):

    def __init__(self,
                 task_training_dict, config_dict, model_path,
                 vocab, num_batches, delta_config, batch_size,
                 model_name="awesome_model", log_dir="cell_annotation_logs/", wandb_config=None,):
        self.logger = scgpt.logger

        self.mask_ratio = task_training_dict["mask_ratio"]
        self.mask_value = task_training_dict["mask_value"]
        task_training_dict.pop("mask_ratio")
        task_training_dict.pop("mask_value")

        self.model = TransformerModel(ntoken=len(vocab), 
                            num_batch_labels=num_batches,
                            vocab=vocab,
                            train_batch_size=batch_size,
                            **config_dict,
                            **task_training_dict)
        self.add_delta_method(delta_config, wandb_config)
        if model_path is not None:
            self.load_model(model_path)

        self.vocab = vocab
        self.pad_token = task_training_dict["pad_token"]
        #self.criterion = nn.CrossEntropyLoss()
        #self.eps = EPS
        self.batch_size = batch_size
        self.log_interval = LOG_INTERVAL
        scgpt.utils.add_file_handler(self.logger, log_dir + f"{model_name}.log")
        #self.max_seq_len = seq_len
        self.model_name = model_name
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.model.to(self.device)
        self.include_zero_gene = False

        # That's if input_embed_style in the config is "continuous"
        self.pad_value = task_training_dict['pad_value']

        self.criterion = None
        self.eps = None
        self.training_config_dict = None
        self.need_predictions = False
        self.need_embeddings = False

    @abstractmethod
    def get_forward_params_for_training(self, batch_data):
        pass

    @abstractmethod
    def get_forward_params_for_evaluation(self, batch_data):
        pass

    @abstractmethod
    def init_loss(self):
        pass            

    @abstractmethod
    def get_loss_for_batch(self, output_dict: dict, batch_data: dict) -> float:
        pass

    @ abstractmethod
    def update_loss(self, loss):
        pass

    @abstractmethod
    def log_training(self, lr, epoch, batch, num_batches):
        pass

    @abstractmethod
    def calc_eval_metrics(self, output, batch_data) -> tuple[float, float]:
        pass

    @abstractmethod
    def get_predictions_from_model_output(self, output):
        pass

    @abstractmethod
    def calc_test_metrics(self, test_results_dict):
        pass

    @abstractmethod
    def init_eval_metrics(self):
        pass

    @abstractmethod
    def evaluate_epoch(self):
        pass

    @abstractmethod
    def get_val_loss_for_comparison(self, resuls_dict):
        pass

    def test_mode(self):
        self.mask_ratio = 0.0

    def load_model(self, model_path):
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

    
    def _evaluate(self, loader: DataLoader, return_raw: bool = False, return_embs=False,
                  get_intermediate_outputs=False, get_gene_embs=False, get_attn_maps=False) -> float:
        """
        Evaluate the model on the evaluation data.
        """
        self.model.eval()
        total_loss = 0.0
        total_error = 0.0
        predictions = []
        intermediate_embeddings = [[] for _ in range(self.model.nlayers + 1)]
        # avg_attention_maps holds the average attention for each layer across the entire
        # test dataset.
        avg_attention_maps = [[] for _ in range(self.model.nlayers)]
        # cell_attention_maps holds the attention between the cell and all the genes,
        # for each cell in the dataset individually.
        cell_attention_maps = [[] for _ in range(self.model.nlayers)]
        embeddings = []
        self.init_eval_metrics()
        total_num = 0
        with torch.no_grad():
            for batch_data in loader:
                with torch.cuda.amp.autocast(enabled=True):
                    forward_input_dict = self.get_forward_params_for_evaluation(batch_data)
                    forward_input_dict["get_intermediate_outputs"] = get_intermediate_outputs
                    forward_input_dict["get_gene_embs"] = get_gene_embs
                    forward_input_dict["get_attention_maps"] = get_attn_maps
                    output = self.model(
                        **forward_input_dict
                    )

                if get_intermediate_outputs:
                    if get_gene_embs:
                        # For the first batch
                        if len(intermediate_embeddings[0]) == 0:
                            # output is a list of length num_layers + 1,
                            # and each element stores a tensor of size (num_genes * 512),
                            # of SUMS of the genes per batch.
                            for i in range(self.model.nlayers + 1):
                                intermediate_embeddings[i] = output[i].cpu().numpy()
                        else:
                            for i in range(self.model.nlayers + 1):
                                intermediate_embeddings[i] += output[i].cpu().numpy()

                    # Get the cell embeddings
                    else:
                        # Add 1 for the embeddings after model._encoder
                        for i in range(self.model.nlayers + 1):
                            intermediate_embeddings[i].append(output[i].cpu().numpy())
                elif get_attn_maps:
                    for i in range(self.model.nlayers):
                        if get_gene_embs:
                            if len(avg_attention_maps[i]) == 0:
                                    avg_attention_maps[i] = output[i].cpu().numpy()
                            else:
                                avg_attention_maps[i] += output[i].cpu().numpy()
                        else:
                            # get the attention of the <cls>
                            cell_attention_maps[i].append(output[i].cpu().numpy())
                else:
                    self.calc_eval_metrics(output, batch_data)
                    #total_loss += total_batch_loss
                    #total_error += total_batch_error
                    if return_raw:
                        predictions.append(self.get_predictions_from_model_output(output))
                    if return_embs:
                        embeddings.append(output["cell_emb"].cpu().numpy())

                total_num += len(batch_data["gene_ids"].to(self.device))

        if get_intermediate_outputs:
            for i in range(len(intermediate_embeddings)):
                if get_gene_embs:
                    # Gene embeddings were summed across the batches, now to get the mean we need
                    # to divide them.
                    intermediate_embeddings[i] /= total_num
                else:
                    intermediate_embeddings[i] = np.concatenate(intermediate_embeddings[i], axis=0)
            return intermediate_embeddings
        
        if get_attn_maps:
            for i in range(self.model.nlayers):
                if get_gene_embs:
                    avg_attention_maps[i] /= total_num
                else:
                    cell_attention_maps[i] = np.concatenate(cell_attention_maps[i], axis=0)
            if get_gene_embs:
                return avg_attention_maps
            else:
                return cell_attention_maps

        eval_dict = {}
        if return_raw:
            eval_dict["predictions"] = np.concatenate(predictions, axis=0)

        if return_embs:
            eval_dict["embeddings"] = np.concatenate(embeddings, axis=0)

        eval_results = self.evaluate_epoch()
        eval_dict.update(eval_results)
        return eval_dict


    def test(self, intermediate_embeddings_file=None,
             predictions_file=None, embeddings_file=None, get_gene_embs=False,
             attention_maps_file=None, celltype_list=[]):
        self.include_zero_gene = True
        test_loader = self.data_loader.get_test_loader(self.include_zero_gene, celltype_list)
        adata_test = self.data_loader.adata_test

        if intermediate_embeddings_file:
            embeddings = self._evaluate(
                loader=test_loader,
                get_intermediate_outputs=True,
                get_gene_embs=get_gene_embs,
            )
            if get_gene_embs:
                embeddings_adata = ad.AnnData(obs=adata_test.var.copy())
            else:
                embeddings_adata = ad.AnnData(obs=adata_test.obs.copy())
            for i in range(len(embeddings)):
                embeddings_adata.obsm[f"transformer_layer_{i}"] = embeddings[i]
            embeddings_adata.write_h5ad(intermediate_embeddings_file)
            return
        
        if attention_maps_file:
            if get_gene_embs:
                avg_attn_maps = self._evaluate(
                    loader=test_loader,
                    get_intermediate_outputs=False,
                    get_gene_embs=True,
                    get_attn_maps=True
                )
            else:
                cell_attn_maps = self._evaluate(
                loader=test_loader,
                get_intermediate_outputs=False,
                get_attn_maps=True
            )
            var_with_cls = adata_test.var.copy()
            if "highly_variable" in var_with_cls:
                cols_to_drop = ["highly_variable",
                                "highly_variable_intersection",
                                "feature_is_filtered"]
                var_with_cls.drop(columns=[col for col in cols_to_drop if col in var_with_cls.columns], 
                                  inplace=True)
            cls_row = pd.DataFrame(index=['cls'])
            var_with_cls = pd.concat([cls_row, var_with_cls])

            if get_gene_embs:
                avg_attn_adata = ad.AnnData(obs=var_with_cls)
                for i in range(self.model.nlayers):
                    avg_attn_adata.obsm[f"transformer_layer_{i}"] = avg_attn_maps[i]
                avg_attn_file = attention_maps_file + "_avg.h5ad"
                avg_attn_adata.write_h5ad(avg_attn_file)
            else:
                cell_attn_adata = ad.AnnData(obs=adata_test.obs.copy(), var=var_with_cls)
                for i in range(self.model.nlayers):
                    cell_attn_adata.obsm[f"transformer_layer_{i}"] = cell_attn_maps[i]
                cell_attn_file = attention_maps_file + "_cell.h5ad"
                cell_attn_adata.write_h5ad(cell_attn_file)
            return

        save_predictions = predictions_file and self.need_predictions
        save_embeddings = embeddings_file and self.need_embeddings

        test_eval_dict = self._evaluate(loader=test_loader,
                                        return_raw=self.need_predictions,
                                        return_embs=self.need_embeddings)
        
        if self.need_predictions:
            adata_test.obs["predictions"] = test_eval_dict["predictions"]
        if self.need_embeddings:
            adata_test.obsm["X_scGPT"] = test_eval_dict["embeddings"]

        if embeddings_file:
            adata_test.write_h5ad(embeddings_file)

        #TODO: clean up this code so it will do something useful
        
        #test_eval_dict["celltype"] = adata.obs["celltype"]
        #test_eval_dict["str_batch"] = adata.obs["str_batch"]

        test_metrics = self.calc_test_metrics(adata_test)
        #wandb.log(test_metrics)

        if save_predictions or save_embeddings:
            predictions_df = adata_test.obs.copy()       
            if save_predictions:
                predictions_df["predicted_celltype_id"] = test_eval_dict["predictions"]   
            """      
            if save_embeddings:
                predictions_adata = ad.AnnData(obs=predictions_df)
                predictions_adata.obsm["X_emb"] = embeddings
                predictions_adata.write_h5ad(embeddings_file)"""
                    
        return test_metrics


    def _train_step(self, data_loader: DataLoader, optimizer, scheduler, scaler, epoch):

        self.model.train()
        self.init_loss()

        num_batches = len(data_loader)
        for batch, batch_data in enumerate(data_loader):
            with torch.cuda.amp.autocast(enabled=True):
                output_dict = self.model(
                    **self.get_forward_params_for_training(batch_data)
                )
                loss = self.get_loss_for_batch(output_dict, batch_data)

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

            total_loss_in_epoch = self.update_loss(loss.item())

            if batch % self.log_interval == 0 and batch > 0:
                self.log_training(scheduler.get_last_lr()[0], epoch, batch, num_batches)
            scheduler.step()
        return total_loss_in_epoch / num_batches


    def train(self, lr, num_epochs, find_lr=False, warm_up_percentage=0, early_stop=False):

        self.data_loader.prepare_train_and_valid(self.include_zero_gene)
        best_val_loss = float("inf")

        optimizer = torch.optim.Adam(
            self.model.parameters(), lr=lr, eps=self.eps
        )
        steps_in_epoch = self.data_loader.get_num_training_steps()
        if find_lr:
            assert warm_up_percentage == 0, "Warm up epochs are not supported in lr finding mode."
            # If in the lr finding mode.
            scheduler = torch.optim.lr_scheduler.StepLR(
                optimizer, steps_in_epoch, gamma=FIND_LR_GAMMA
            )
        else:
            num_training_steps = num_epochs * steps_in_epoch
            warm_up_steps = (warm_up_percentage / 100.) * num_training_steps
            scheduler = get_linear_schedule_with_warmup(
                optimizer,
                num_warmup_steps=warm_up_steps,
                num_training_steps=num_epochs * steps_in_epoch,
            )
        move_optimizer_params_to_cuda(optimizer)
        scaler = torch.cuda.amp.GradScaler(enabled=True)

        # for early stopping
        last_val_losses = np.zeros(EARLY_STOPPING_EPOCHS_AVG)
        early_stopping_counter = 0
        prev_epoch_loss = np.inf
        curr_val_loss_avg = np.inf
        for epoch in range(1, num_epochs + 1):
            epoch_start_time = time.time()
            train_loader, valid_loader = self.data_loader.get_train_valid_data_per_epoch(
                self.sort_batches,
                self.train_dataloader_params_dict,
                self.test_dataloader_params_dict
            )
            epoch_loss = self._train_step(train_loader, optimizer, scheduler, scaler, epoch)
            last_lr = float(scheduler.get_last_lr()[0])
            if find_lr:
                with open(LR_FINDER_LOG_DIR + self.model_name, 'a') as f:
                    f.write(f"{last_lr} {epoch_loss}\n")
            eval_results = self._evaluate(loader=valid_loader)

            if not find_lr:
                test_results = self.test()

                self.logger.info(
                    f"Test results: {test_results}"
                )

            elapsed = time.time() - epoch_start_time
            self.log_epoch(epoch, elapsed, eval_results)

            # Early stopping mechanism
            # get the loss used for comparing performance of model epochs.
            val_loss = self.get_val_loss_for_comparison(eval_results)
            prev_val_loss_avg = curr_val_loss_avg
            last_val_losses[(epoch - 1) % EARLY_STOPPING_EPOCHS_AVG] = val_loss
            curr_val_loss_avg = np.sum(last_val_losses) / np.minimum(EARLY_STOPPING_EPOCHS_AVG, epoch)
            if curr_val_loss_avg >= prev_val_loss_avg:
                early_stopping_counter += 1
            else:
                early_stopping_counter = 0
            if early_stop and early_stopping_counter >= EARLY_STOPPING_PATIENCE:
                # return the value for optuna hyperparameter search
                self.logger.info(
                    f"Early stopping at epoch {epoch} with score {curr_val_loss_avg:5.4f}"
                )
                return

            if find_lr and (val_loss > prev_val_loss_avg):
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
            """
            wandb.log(
                {
                    "train_loss": epoch_loss,
                    "best_val_loss": best_val_loss,
                    "val_loss": val_loss,
                    "epoch": epoch,
                    "lr": last_lr,
                }
            )"""


    ####################################################################################################################
    ############################################# Delta Modules ########################################################
    ####################################################################################################################


    def freeze_all_modules(self):
        for param in self.model.parameters():
            param.requires_grad = False

    def finetune_custom_module(self, delta_model_config, wandb_config):
        self.freeze_all_modules()
        module_name = delta_model_config["module"]
        print(module_name)
        for param in getattr(self.model, module_name).parameters():
            param.requires_grad = True
        
        wandb_config["trainable module"] = module_name


    def train_transformer_modules(self, delta_model_config: dict, wandb_config: dict):
        """
        Freeze all of the modules except the specified transformer modules in the config.
        """
        # train_modules: list, sublayers: list = None
        if not delta_model_config.get("train_modules", None):
            raise ValueError("train_modules must be specified in the delta model config for 'specify_transformer_layers' delta type.")
        
        self.freeze_all_modules()
        train_modules = delta_model_config["train_modules"]
        sublayers = delta_model_config.get("sublayers", None)
        wandb_config["train_modules"] = train_modules
        wandb_config["sublayers"] = sublayers if sublayers is not None else []

        modules_to_train = ["transformer_encoder.layers." + str(i) for i in train_modules]
        if sublayers is not None:
            modules_to_train = [f"{module}.{sublayer}" for module in modules_to_train for sublayer in sublayers]

        for name, module in self.model.named_modules():
            if name.startswith(tuple(modules_to_train)):
                for param in module.parameters():
                    param.requires_grad = True
                self.logger.info(f"Un-freezing module {name}")


    def add_delta_method(self, delta_model_config: dict, wandb_config: dict):
        """
        Add the delta method to the model.
        :param cam: CellAnnotationModelWrapper
        :param model_config: model configuration
        :param wandb_config: wandb configuration
        """
        SUPPORTED_DELTA_TYPES = ["all_weights", "specify_transformer_layers"]
        delta_method = delta_model_config["delta_method"]
        
        if delta_method == "all_weights":
            # nothing to do, all weights are trainable
            print("Delta type is 'all_weights', no additional configuration needed.")
            pass
        elif delta_method == "specify_transformer_layers":
            print("Delta type is 'specify_transformer_layers', training specified transformer modules.")
            self.train_transformer_modules(delta_model_config, wandb_config)
        elif delta_method == "custom_module":
            print("Delta type: custom_module, for module: ", delta_model_config["module"])
            self.finetune_custom_module(delta_model_config, wandb_config)
        else:
            raise ValueError(f"Unsupported delta type: {delta_method}. Supported types are: {SUPPORTED_DELTA_TYPES}")
        
        trainable_params = self.model.get_trainable_parameters()
        wandb_config["delta_method"] = delta_method
        wandb_config["trainable_params"] = trainable_params
        #print(self.model)


class CellAnnotation(ScGPTModelWrapper):
    
    def __init__(self,
                 task_training_dict,
                 config_dict,
                 model_path,
                 vocab,
                 delta_config, batch_size,
                 n_input_bins, n_hvg, data_name,
                 model_name="awesome_model",
                 log_dir="cell_annotation_logs/", wandb_config=None):
        
        seq_len = config_dict["seq_len"]
        self.data_loader = SimpleDataloader(data_name, vocab, n_input_bins, n_hvg,
                                            pad_token=task_training_dict["pad_token"],
                                            pad_value=task_training_dict["pad_value"],
                                            batch_size=batch_size,
                                            max_seq_len=seq_len,
                                            model_task="celltype_annotation",
                                            mask_ratio=task_training_dict["mask_ratio"],
                                            mask_value=task_training_dict["mask_value"])
        config_dict["n_cls"] = self.data_loader.get_num_celltypes()
        num_batches = self.data_loader.get_num_batches()
        
        super().__init__(
                         task_training_dict=task_training_dict,
                         config_dict=config_dict,
                         model_path=model_path,
                         vocab=vocab, 
                         num_batches=num_batches,
                         delta_config=delta_config, 
                         log_dir=log_dir,
                         model_name=model_name,
                         wandb_config=wandb_config,
                         batch_size=batch_size
        )

        self.criterion = nn.CrossEntropyLoss()
        self.eps = 1e-8
        self.need_predictions = True
        self.training_config_dict = {
            'batch_labels': None,
            'CLS': True
        }
        self.sort_batches = False
        self.train_dataloader_params_dict = {
            'shuffle': False,
            'drop_last': False
        }
        self.test_dataloader_params_dict = {
            'shuffle': False,
            'drop_last': False
        }

    def get_forward_params_for_training(self, batch_data):
        input_gene_ids = batch_data["gene_ids"].to(self.device)
        input_values = batch_data["values"].to(self.device)
        src_key_padding_mask = input_gene_ids.eq(self.vocab[self.pad_token])

        forward_input_dict = {
            "input_ids": input_gene_ids,
            "inputs_embeds": input_values,
            "src_key_padding_mask": src_key_padding_mask,
            "CLS": True,
            "do_sample": False,
            "batch_labels": None
        }

        return forward_input_dict
    
    def get_forward_params_for_evaluation(self, batch_data):
        return self.get_forward_params_for_training(batch_data)
        

    def unfreeze_classifier(self):
        for param in self.model.cls_decoder.parameters():
            param.requires_grad = True

    def add_delta_method(self, delta_model_config: dict, wandb_config: dict):
        if delta_model_config["delta_method"] == "classifier":
            print("Delta type is 'classifier', finetuning the classifier decoder.")
            wandb_config["delta_method"] = "classifier"
            self.freeze_all_modules()
        else:
            super().add_delta_method(delta_model_config, wandb_config)
        self.unfreeze_classifier()
        self.model.get_trainable_parameters()

    def init_loss(self):
        self.total_loss = 0.0
        self.total_error = 0.0
        self.start_time = time.time()
        self.total_loss_in_epoch = 0.0
        
    def init_eval_metrics(self):
        self.eval_loss = 0.0
        self.eval_err = 0.0
        self.total_num_in_eval = 0

    def get_loss_for_batch(self, output_dict: dict, batch_data: dict) -> float:
        metrics_to_log = {}
        celltype_labels = batch_data["celltype_labels"].to(self.device)
        self.loss_cls = self.criterion(output_dict["cls_output"], celltype_labels)
        metrics_to_log.update({"train/cls": self.loss_cls.item()})

        self.error_rate = 1 - (
            (output_dict["cls_output"].argmax(1) == celltype_labels)
            .sum()
            .item()
        ) / celltype_labels.size(0)

        return self.loss_cls

    def update_loss(self, loss) -> float:
        #wandb.log({"train/loss": loss, "train/err": self.error_rate})

        self.total_loss += loss
        self.total_loss_in_epoch += loss
        self.total_error += self.error_rate

        return self.total_loss_in_epoch

    def log_training(self, lr, epoch, batch, num_batches):
        ms_per_batch = (time.time() - self.start_time) * 1000 / self.log_interval
        cur_loss = self.total_loss / self.log_interval
        
        cur_error = self.total_error / self.log_interval
        self.logger.info(
            f"| epoch {epoch:3d} | {batch:3d}/{num_batches:3d} batches | "
            f"lr {lr:05.8f} | ms/batch {ms_per_batch:5.2f} | "
            f"loss {cur_loss:5.2f} | "
            + (f"err {cur_error:5.2f} | ")
        )     
        self.total_loss = 0
        self.total_error = 0
        self.start_time = time.time()
    
    def calc_eval_metrics(self, output, batch_data):
        output_values = output["cls_output"]
        celltype_labels = batch_data["celltype_labels"].to(self.device)
        input_gene_ids = batch_data["gene_ids"].to(self.device)
        loss = self.criterion(output_values, celltype_labels)
        self.eval_loss += loss.item() * len(input_gene_ids)
        accuracy = (output_values.argmax(1) == celltype_labels).sum().item()
        self.eval_err += (1 - accuracy / len(input_gene_ids)) * len(input_gene_ids)
        self.total_num_in_eval += len(input_gene_ids)

    def evaluate_epoch(self):
        results = {}
        results["total_loss"] = self.eval_loss / self.total_num_in_eval
        results["total_err"] = self.eval_err / self.total_num_in_eval
        # the values should reset with a call to init_eval_metrics()
        return results

    def log_epoch(self, epoch, elapsed, eval_results):
        self.logger.info("-" * 89)
        self.logger.info(
            f"| end of epoch {epoch:3d} | time: {elapsed:5.2f}s | "
            f"valid loss/mse {eval_results['total_loss']:5.4f} | err {eval_results['total_err']:5.4f}"
        )
        self.logger.info("-" * 89)

    def get_val_loss_for_comparison(self, resuls_dict):
        return resuls_dict["total_loss"]

    
    def get_predictions_from_model_output(self, output):
        output_values = output["cls_output"]
        return output_values.argmax(1).cpu().numpy()
    
    def calc_test_metrics(self, adata):

        celltypes_labels = np.array(adata.obs["celltype_id"].tolist())
        predictions = np.array(adata.obs["predictions"].tolist())


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

        #wandb.log(results)

        return results


class BatchCorrection(ScGPTModelWrapper):

    def __init__(self, task_training_dict, config_dict, model_path, vocab, 
                 delta_config, batch_size,
                 n_hvg, n_input_bins, data_name,
                 model_name="awesome_model",
                 log_dir="cell_annotation_logs/", wandb_config=None):
        print("data name:", data_name)        
        seq_len = config_dict["seq_len"]
        self.data_loader = SimpleDataloader(data_name, vocab, n_input_bins, n_hvg,
                                            pad_token=task_training_dict["pad_token"],
                                            pad_value=task_training_dict["pad_value"],
                                            batch_size=batch_size,
                                            max_seq_len=seq_len,
                                            model_task="batch_correction",
                                            mask_ratio=task_training_dict["mask_ratio"],
                                            mask_value=task_training_dict["mask_value"])
        config_dict["n_cls"] = self.data_loader.get_num_celltypes()
        num_batches = self.data_loader.get_num_batches()
   
        super().__init__(
                         task_training_dict=task_training_dict,
                         config_dict=config_dict,
                         model_path=model_path,
                         vocab=vocab, 
                         num_batches=num_batches,
                         delta_config=delta_config,
                         batch_size=batch_size,
                         log_dir=log_dir,
                         model_name=model_name, 
                         wandb_config=wandb_config,
                         )

        # Different functions for logging purposes
        self.criterion = masked_mse_loss
        self.criterion_dab = nn.CrossEntropyLoss()
        self.criterion_neg_log_bernoulli = criterion_neg_log_bernoulli

        self.eps = 1e-4
        self.include_zero_gene = True
        self.sort_batches = True
        self.need_embeddings = True

        self.train_dataloader_params_dict = {
            'shuffle': False,
            'drop_last': False,
            'intra_domain_shuffle': True,
            'per_seq_batch_sample': True
        }
        self.test_dataloader_params_dict = {
            'shuffle': False,
            'drop_last': False,
            'intra_domain_shuffle': False,
            'per_seq_batch_sample': True
        }

        self.training_config_dict = {
            'batch_labels': None,
            'CLS': False
        }

    def unfreeze_batch_correction_modules(self):
        for param in self.model.batch_encoder.parameters():
            param.requires_grad = True
        # Adversarial descriminator to predict which batch a cell is coming from.
        for param in self.model.grad_reverse_discriminator.parameters():
            param.requires_grad = True
        # Decoder for the masked value prediction for cell embeddings.
        for param in self.model.mvc_decoder.parameters():
            param.requires_grad = True
        # Domain specific batch norm
        for param in self.model.dsbn.parameters():
            param.requires_grad = True

    def add_delta_method(self, delta_model_config: dict, wandb_config: dict):
        if delta_model_config["delta_method"] == "classifier":
            print("Delta type is 'classifier', finetuning the classifier decoder.")
            wandb_config["delta_method"] = "classifier"
            self.freeze_all_modules()
        else:
            super().add_delta_method(delta_model_config, wandb_config)
        self.unfreeze_batch_correction_modules()
        self.model.get_trainable_parameters()

    def get_forward_params_for_training(self, batch_data):
        input_gene_ids = batch_data["gene_ids"].to(self.device)
        input_values = batch_data["values"].to(self.device)
        batch_labels = batch_data["batch_labels"].to(self.device)

        src_key_padding_mask = input_gene_ids.eq(self.vocab[self.pad_token])

        forward_input_dict = {
            "input_ids": input_gene_ids,
            "inputs_embeds": input_values,
            "src_key_padding_mask": src_key_padding_mask,
            "batch_labels": batch_labels,
            "MVC": True,
            "ECS": True,
        }

        return forward_input_dict
    
    def get_forward_params_for_evaluation(self, batch_data):
        eval_input_dict = self.get_forward_params_for_training(batch_data)
        eval_input_dict["MVC"] = False
        eval_input_dict["ECS"] = False
        return eval_input_dict

    def init_loss(self):
        self.start_time = time.time()            
        self.total_loss_in_epoch = 0.0
        self.metrics_lo_log = {}
        self.losses_dict = {
            "mse": 0.0,
            "gepc": 0.0,
            "error": 0.0,
            "total": 0.0
        }

    def get_loss_for_batch(self, output_dict: dict, batch_data: dict) -> float:
        input_values = batch_data["values"].to(self.device)
        target_values = batch_data["target_values"].to(self.device)
        batch_labels = batch_data["batch_labels"].to(self.device)

        masked_positions = input_values.eq(self.mask_value)  # the postions to predict

        loss = loss_mse = self.criterion(
            output_dict["mlm_output"], target_values, masked_positions
        )
        self.losses_dict["mse"] += loss.item()
        metrics_to_log = {"train/mse": loss_mse.item()}
        loss_zero_log_prob = self.criterion_neg_log_bernoulli(
            output_dict["mlm_zero_probs"], target_values, masked_positions
        )
        loss = loss + loss_zero_log_prob
        metrics_to_log.update({"train/nzlp": loss_zero_log_prob.item()})
        loss_gepc = self.criterion(
            output_dict["mvc_output"], target_values, masked_positions
        )
        self.losses_dict["gepc"] += loss_gepc.item()
        loss = loss + loss_gepc
        metrics_to_log.update({"train/mvc": loss_gepc.item()})

        loss_gepc_zero_log_prob = self.criterion_neg_log_bernoulli(
            output_dict["mvc_zero_probs"], target_values, masked_positions
        )
        loss = loss + loss_gepc_zero_log_prob
        metrics_to_log.update(
            {"train/mvc_nzlp": loss_gepc_zero_log_prob.item()}
        )
        loss_ecs = 10 * output_dict["loss_ecs"]
        loss = loss + loss_ecs
        metrics_to_log.update({"train/ecs": loss_ecs.item()})
        loss_dab = self.criterion_dab(output_dict["dab_output"], batch_labels)
        #TODO: set dab_weight = 1.0 as a parameter in a config
        loss = loss + 1.0 * loss_dab
        self.losses_dict["total"] += loss

        with torch.no_grad():
            mre = masked_relative_error(
                output_dict["mlm_output"], target_values, masked_positions
            )
        self.losses_dict["error"] += mre.item()
        
        return loss


    def update_loss(self, loss) -> float:
        self.total_loss_in_epoch += loss
        return self.total_loss_in_epoch

    def log_training(self, lr, epoch, batch, num_batches):
        ms_per_batch = (time.time() - self.start_time) * 1000 / self.log_interval
        cur_loss = self.losses_dict["total"] / self.log_interval
        cur_mse = self.losses_dict["mse"] / self.log_interval
        cur_gepc = self.losses_dict["gepc"] / self.log_interval
        cur_error = self.losses_dict["error"] / self.log_interval
        self.logger.info(
            f"| epoch {epoch:3d} | {batch:3d}/{num_batches:3d} batches | "
            f"lr {lr:05.8f} | ms/batch {ms_per_batch:5.2f} | "
            f"loss {cur_loss:5.2f} | mse {cur_mse:5.2f} | mre {cur_error:5.2f} |"
            + (f"gepc {cur_gepc:5.2f} |")
        )
        self.start_time = time.time()
        for key in self.losses_dict.keys():
            self.losses_dict[key] = 0.0
        

    def calc_test_metrics(self, adata):
        adata_t = adata.copy()
        
        results = {}
        try:
            results = eval_scib_metrics(adata_t)
        except Exception as e:
            traceback.print_exc()
            self.logger.error(e)

        return results
    
    def init_eval_metrics(self):
        self.eval_loss = 0.0
        self.eval_err = 0.0
        self.total_num_in_eval = 0
    
    def calc_eval_metrics(self, output, batch_data) -> tuple[float, float]:
        input_gene_ids = batch_data["gene_ids"].to(self.device)
        target_values = batch_data["target_values"].to(self.device)
        input_values = batch_data["values"].to(self.device)
        masked_positions = input_values.eq(self.mask_value)

        output_values = output["mlm_output"]
        loss = self.criterion(
            output_values, target_values, masked_positions
            ).item() * len(input_gene_ids)
        error = masked_relative_error(
                output_values, target_values, masked_positions
            ).item() * len(input_gene_ids)
        self.eval_loss += loss
        self.eval_err += error
        self.total_num_in_eval += len(input_gene_ids)

    def log_epoch(self, epoch, elapsed, eval_results):
        self.logger.info("-" * 89)
        self.logger.info(
            f"| end of epoch {epoch:3d} | time: {elapsed:5.2f}s | "
            f"valid loss/mse {eval_results['total_loss']:5.4f} | err {eval_results['total_err']:5.4f}"
        )
        self.logger.info("-" * 89)

    def evaluate_epoch(self):
        results = {}
        results["total_loss"] = self.eval_loss / self.total_num_in_eval
        results["total_err"] = self.eval_err / self.total_num_in_eval
        # the values should reset with a call to init_eval_metrics()
        return results

    def get_val_loss_for_comparison(self, resuls_dict):
        return resuls_dict["total_loss"]

    def get_predictions_from_model_output(self, output):
        return np.arange(10)
    

class GEPrediction(ScGPTModelWrapper):
    def __init__(self, task_training_dict, config_dict, model_path, vocab, 
                 delta_config, batch_size,
                 n_hvg, n_input_bins, data_name,
                 model_name="awesome_model",
                 log_dir="cell_annotation_logs/", wandb_config=None):
        print("data name:", data_name)        
        seq_len = config_dict["seq_len"]
        self.data_loader = SimpleDataloader(data_name, vocab, n_input_bins, n_hvg,
                                            pad_token=task_training_dict["pad_token"],
                                            pad_value=task_training_dict["pad_value"],
                                            batch_size=batch_size,
                                            max_seq_len=seq_len,
                                            model_task="ge_prediction",
                                            mask_ratio=task_training_dict["mask_ratio"],
                                            mask_value=task_training_dict["mask_value"])
        config_dict["n_cls"] = self.data_loader.get_num_celltypes()
        num_batches = self.data_loader.get_num_batches()
   
        super().__init__(
                         task_training_dict=task_training_dict,
                         config_dict=config_dict,
                         model_path=model_path,
                         vocab=vocab, 
                         num_batches=num_batches,
                         delta_config=delta_config,
                         batch_size=batch_size,
                         log_dir=log_dir,
                         model_name=model_name, 
                         wandb_config=wandb_config,
                         )

        # Different functions for logging purposes
        self.criterion = masked_mse_loss
        self.criterion_neg_log_bernoulli = criterion_neg_log_bernoulli

        self.eps = 1e-4
        self.include_zero_gene = True
        self.sort_batches = False
        self.need_embeddings = True

        self.train_dataloader_params_dict = {
            'shuffle': False,
            'drop_last': False,
            'intra_domain_shuffle': True,
            'per_seq_batch_sample': False
        }
        self.test_dataloader_params_dict = {
            'shuffle': False,
            'drop_last': False,
            'intra_domain_shuffle': False,
            'per_seq_batch_sample': False
        }

        self.training_config_dict = {
            'batch_labels': None,
            'CLS': False
        }

    def unfreeze_prediction_heads(self):
        # Decoder for the masked value prediction for cell embeddings.
        for param in self.model.mvc_decoder.parameters():
            param.requires_grad = True

    def add_delta_method(self, delta_model_config: dict, wandb_config: dict):
        if delta_model_config["delta_method"] == "classifier":
            print("Delta type is 'classifier', finetuning the classifier decoder.")
            wandb_config["delta_method"] = "classifier"
            self.freeze_all_modules()
        else:
            super().add_delta_method(delta_model_config, wandb_config)
        self.unfreeze_prediction_heads()
        self.model.get_trainable_parameters()

    def get_forward_params_for_training(self, batch_data):
        input_gene_ids = batch_data["gene_ids"].to(self.device)
        input_values = batch_data["values"].to(self.device)

        src_key_padding_mask = input_gene_ids.eq(self.vocab[self.pad_token])

        forward_input_dict = {
            "input_ids": input_gene_ids,
            "inputs_embeds": input_values,
            "src_key_padding_mask": src_key_padding_mask,
            "MVC": True,
            "ECS": True,
        }

        return forward_input_dict
    
    def get_forward_params_for_evaluation(self, batch_data):
        eval_input_dict = self.get_forward_params_for_training(batch_data)
        eval_input_dict["MVC"] = False
        eval_input_dict["ECS"] = False
        return eval_input_dict

    def init_loss(self):
        self.start_time = time.time()            
        self.total_loss_in_epoch = 0.0
        self.metrics_lo_log = {}
        self.losses_dict = {
            "mse": 0.0,
            "gepc": 0.0,
            "error": 0.0,
            "total": 0.0
        }

    def get_loss_for_batch(self, output_dict: dict, batch_data: dict) -> float:
        input_values = batch_data["values"].to(self.device)
        target_values = batch_data["target_values"].to(self.device)
        masked_positions = input_values.eq(self.mask_value)  # the postions to predict

        loss = loss_mse = self.criterion(
            output_dict["mlm_output"], target_values, masked_positions
        )
        self.losses_dict["mse"] += loss.item()
        metrics_to_log = {"train/mse": loss_mse.item()}
        loss_zero_log_prob = self.criterion_neg_log_bernoulli(
            output_dict["mlm_zero_probs"], target_values, masked_positions
        )
        loss = loss + loss_zero_log_prob
        metrics_to_log.update({"train/nzlp": loss_zero_log_prob.item()})
        loss_gepc = self.criterion(
            output_dict["mvc_output"], target_values, masked_positions
        )
        self.losses_dict["gepc"] += loss_gepc.item()
        loss = loss + loss_gepc
        metrics_to_log.update({"train/mvc": loss_gepc.item()})

        loss_gepc_zero_log_prob = self.criterion_neg_log_bernoulli(
            output_dict["mvc_zero_probs"], target_values, masked_positions
        )
        loss = loss + loss_gepc_zero_log_prob
        metrics_to_log.update(
            {"train/mvc_nzlp": loss_gepc_zero_log_prob.item()}
        )
        loss_ecs = 10 * output_dict["loss_ecs"]
        loss = loss + loss_ecs
        metrics_to_log.update({"train/ecs": loss_ecs.item()})
        self.losses_dict["total"] += loss

        with torch.no_grad():
            mre = masked_relative_error(
                output_dict["mlm_output"], target_values, masked_positions
            )
        self.losses_dict["error"] += mre.item()
        
        return loss


    def update_loss(self, loss) -> float:
        self.total_loss_in_epoch += loss
        return self.total_loss_in_epoch

    def log_training(self, lr, epoch, batch, num_batches):
        ms_per_batch = (time.time() - self.start_time) * 1000 / self.log_interval
        cur_loss = self.losses_dict["total"] / self.log_interval
        cur_mse = self.losses_dict["mse"] / self.log_interval
        cur_gepc = self.losses_dict["gepc"] / self.log_interval
        cur_error = self.losses_dict["error"] / self.log_interval
        self.logger.info(
            f"| epoch {epoch:3d} | {batch:3d}/{num_batches:3d} batches | "
            f"lr {lr:05.8f} | ms/batch {ms_per_batch:5.2f} | "
            f"loss {cur_loss:5.2f} | mse {cur_mse:5.2f} | mre {cur_error:5.2f} |"
            + (f"gepc {cur_gepc:5.2f} |")
        )
        self.start_time = time.time()
        for key in self.losses_dict.keys():
            self.losses_dict[key] = 0.0
        

    def calc_test_metrics(self, adata):
        return "All Good!"
    
    def init_eval_metrics(self):
        self.eval_loss = 0.0
        self.eval_err = 0.0
        self.total_num_in_eval = 0
    
    def calc_eval_metrics(self, output, batch_data) -> tuple[float, float]:
        input_gene_ids = batch_data["gene_ids"].to(self.device)
        target_values = batch_data["target_values"].to(self.device)
        input_values = batch_data["values"].to(self.device)
        masked_positions = input_values.eq(self.mask_value)

        output_values = output["mlm_output"]
        loss = self.criterion(
            output_values, target_values, masked_positions
            ).item() * len(input_gene_ids)
        error = masked_relative_error(
                output_values, target_values, masked_positions
            ).item() * len(input_gene_ids)
        self.eval_loss += loss
        self.eval_err += error
        self.total_num_in_eval += len(input_gene_ids)

    def log_epoch(self, epoch, elapsed, eval_results):
        self.logger.info("-" * 89)
        self.logger.info(
            f"| end of epoch {epoch:3d} | time: {elapsed:5.2f}s | "
            f"valid loss/mse {eval_results['total_loss']:5.4f} | err {eval_results['total_err']:5.4f}"
        )
        self.logger.info("-" * 89)

    def evaluate_epoch(self):
        results = {}
        results["total_loss"] = self.eval_loss / self.total_num_in_eval
        results["total_err"] = self.eval_err / self.total_num_in_eval
        # the values should reset with a call to init_eval_metrics()
        return results

    def get_val_loss_for_comparison(self, resuls_dict):
        return resuls_dict["total_loss"]

    def get_predictions_from_model_output(self, output):
        return np.arange(10)
    

class Perturbation(ScGPTModelWrapper):

    def __init__(self, task_training_dict, config_dict, model_path, vocab, 
                 delta_config, batch_size,
                 n_hvg, n_input_bins, data_name,
                 model_name="awesome_model",
                 log_dir="cell_annotation_logs/", wandb_config=None):
        print("data name:", data_name)        
        seq_len = config_dict["seq_len"]
        self.data_loader = PertDataloader(data_name, vocab, n_input_bins, n_hvg,
                                            pad_token=task_training_dict["pad_token"],
                                            pad_value=task_training_dict["pad_value"],
                                            batch_size=batch_size,
                                            max_seq_len=seq_len,
                                            mask_ratio=task_training_dict["mask_ratio"],
                                            mask_value=task_training_dict["mask_value"])
        config_dict["n_cls"] = self.data_loader.get_num_celltypes()
        num_batches = self.data_loader.get_num_batches()
        self.n_genes = self.data_loader.get_num_genes()
        self.sort_batches = False
   
        super().__init__(
                         task_training_dict=task_training_dict,
                         config_dict=config_dict,
                         model_path=model_path,
                         vocab=vocab, 
                         num_batches=num_batches,
                         delta_config=delta_config,
                         batch_size=batch_size,
                         log_dir=log_dir,
                         model_name=model_name, 
                         wandb_config=wandb_config
                         )
        self.eps = 1e-8
        self.train_dataloader_params_dict = {}
        self.test_dataloader_params_dict = {}
        self.max_seq_len = seq_len
        self.criterion = masked_mse_loss
        self.input_gene_ids_for_batch = None
        self.input_gene_ids_for_eval = None
        
    def init_loss(self):
        self.total_loss = 0.0
        self.start_time = time.time()
        self.total_loss_in_epoch = 0.0

    def get_forward_params_for_training(self, batch_data):
        batch_data.to(self.device)
        x: torch.Tensor = batch_data.x  # (batch_size * n_genes, 2)
        ori_gene_values = x[:, 0].view(self.batch_size, self.n_genes)
        pert_flags = x[:, 1].long().view(self.batch_size, self.n_genes)

        self.input_gene_ids_for_batch = torch.arange(self.n_genes, device=self.device, dtype=torch.long)
        if len(self.input_gene_ids_for_batch) > self.max_seq_len:
            self.input_gene_ids_for_batch = torch.randperm(len(self.input_gene_ids_for_batch), device=self.device)[
                :self.max_seq_len
            ]
        input_values = ori_gene_values[:, self.input_gene_ids_for_batch]
        input_pert_flags = pert_flags[:, self.input_gene_ids_for_batch]

        gene_ids = self.data_loader.get_gene_ids()
        mapped_input_gene_ids = map_raw_id_to_vocab_id(self.input_gene_ids_for_batch, gene_ids)
        mapped_input_gene_ids = mapped_input_gene_ids.repeat(self.batch_size, 1)

        # src_key_padding_mask = mapped_input_gene_ids.eq(vocab[pad_token])
        src_key_padding_mask = torch.zeros_like(
            input_values, dtype=torch.bool, device=self.device
        )

        training_dict = {
            "input_ids": mapped_input_gene_ids,
            "inputs_embeds": input_values,
            "input_pert_flags": input_pert_flags,
            "src_key_padding_mask": src_key_padding_mask,
        }

        return training_dict
    
    def get_loss_for_batch(self, output_dict, batch_data):
        target_gene_values = batch_data.y  # (batch_size, n_genes)
        target_values = target_gene_values[:, self.input_gene_ids_for_batch]
        #input_gene_ids = torch.arange(self.n_genes, device=self.device, dtype=torch.long)
        output_values = output_dict["mlm_output"]
        x: torch.Tensor = batch_data.x  # (batch_size * n_genes, 2)
        ori_gene_values = x[:, 0].view(self.batch_size, self.n_genes)
        input_values = ori_gene_values[:, self.input_gene_ids_for_batch]
        masked_positions = torch.ones_like(
                input_values, dtype=torch.bool, device=self.device
            )  # Use all
        loss = self.criterion(output_values, target_values, masked_positions)

        return loss
    
    def update_loss(self, loss) -> float:
        self.total_loss += loss
        self.total_loss_in_epoch += loss

        return self.total_loss_in_epoch
    
    def log_training(self, lr, epoch, batch, num_batches):
        ms_per_batch = (time.time() - self.start_time) * 1000 / self.log_interval
        cur_loss = self.total_loss / self.log_interval
        
        self.logger.info(
            f"| epoch {epoch:3d} | {batch:3d}/{num_batches:3d} batches | "
            f"lr {lr:05.8f} | ms/batch {ms_per_batch:5.2f} | "
            f"loss {cur_loss:5.2f} | "
        )     
        self.total_loss = 0
        self.start_time = time.time()

    def get_forward_params_for_evaluation(self, batch):
        batch.to(self.device)
        batch_size = len(batch.pert)
        x: torch.Tensor = batch.x
        ori_gene_values = x[:, 0].view(batch_size, -1)  # (batch_size, n_genes)
        pert_flags = x[:, 1].long().view(batch_size, -1)
        gene_ids = self.data_loader.get_gene_ids()
        assert gene_ids is not None
        input_gene_ids = torch.arange(ori_gene_values.size(1), device=self.device)
        self.input_gene_ids_for_eval = torch.arange(self.n_genes, device=self.device, dtype=torch.long)
        if len(self.input_gene_ids_for_eval) > self.max_seq_len:
            self.input_gene_ids_for_eval = torch.randperm(len(self.input_gene_ids_for_eval), device=self.device)[
                :self.max_seq_len
            ]
        input_values = ori_gene_values[:, self.input_gene_ids_for_eval]
        input_pert_flags = pert_flags[:, self.input_gene_ids_for_eval]

        mapped_input_gene_ids = map_raw_id_to_vocab_id(self.input_gene_ids_for_eval, gene_ids)
        #mapped_input_gene_ids = self.input_gene_ids_for_eval.repeat(batch_size, 1)
        mapped_input_gene_ids = mapped_input_gene_ids.repeat(batch_size, 1)

        src_key_padding_mask = torch.zeros_like(
            input_values, dtype=torch.bool, device=self.device
        )

        eval_params_dict = {
            "input_ids": mapped_input_gene_ids,
            "inputs_embeds": input_values,
            "input_pert_flags": input_pert_flags,
            "src_key_padding_mask": src_key_padding_mask,
            "do_sample": True
        }

        return eval_params_dict
    
    def init_eval_metrics(self):
        self.pert_cat = []
        self.pred = []
        self.truth = []
        self.pred_de = []
        self.truth_de = []
        self.total_num_in_eval = 0
    
    # This is called at the end of batch
    def calc_eval_metrics(self, output, batch) -> tuple[float, float]:
        output_values = output["mlm_output"].float()
        batch_size = len(batch.pert)
        ori_gene_values = batch.x[:, 0].view(batch_size, -1)
        pred_gene_values = torch.zeros_like(ori_gene_values)
        pred_gene_values[:, self.input_gene_ids_for_eval] = output_values
        t = batch.y
        self.pred.extend(pred_gene_values.cpu())
        self.truth.extend(t.cpu())
        self.pert_cat.extend(batch.pert)

        # Differentially expressed genes
        for itr, de_idx in enumerate(batch.de_idx):
            self.pred_de.append(pred_gene_values[itr, de_idx])
            self.truth_de.append(t[itr, de_idx])

        self.total_num_in_eval += batch.x.shape[0]
    
    # this is called at the end of every epoch
    def evaluate_epoch(self):
        results = {}
        results["pert_cat"] = np.array(self.pert_cat)
        self.pred = torch.stack(self.pred)
        self.truth = torch.stack(self.truth)
        results["pred"] = self.pred.detach().cpu().numpy().astype(float)
        results["truth"] = self.truth.detach().cpu().numpy().astype(float)

        self.pred_de = torch.stack(self.pred_de)
        self.truth_de = torch.stack(self.truth_de)
        results["pred_de"] = self.pred_de.detach().cpu().numpy().astype(float)
        results["truth_de"] = self.truth_de.detach().cpu().numpy().astype(float)

        val_metrics = compute_perturbation_metrics(
            results, self.data_loader.get_ctrl_condition()
        )

        truth_tensor = torch.as_tensor(results["truth"], device=self.device)
        pred_tensor = torch.as_tensor(results["pred"], device=self.device)
        masked_positions = torch.ones_like(truth_tensor, dtype=torch.bool)
        loss = self.criterion(pred_tensor, truth_tensor, masked_positions)
        val_metrics["mse"] = loss

        return val_metrics

    def log_epoch(self, epoch, elapsed, eval_results):
        self.logger.info("-" * 89)
        self.logger.info(
            f"| end of epoch {epoch:3d} | time: {elapsed:5.2f}s | "
        )
        self.logger.info(eval_results)
        self.logger.info("-" * 89)

    def calc_test_metrics(self, adata):
        return "All Good!"

    def get_val_loss_for_comparison(self, resuls_dict):
        # The metric taken to evaluate the model is pearson correlation - 
        # the higher the better, so we return the negative value as loss.
        return -resuls_dict["pearson"]
    
    def get_predictions_from_model_output(self, output):
        return np.arange(10)