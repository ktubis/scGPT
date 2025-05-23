# %%
from abc import ABC, abstractmethod
import pandas as pd
import anndata as ad
import numpy as np
from pathlib import Path
from enum import Enum
import warnings


class SupportedDatasets(Enum):
    FILTERED_PANCREAS = "filtered_pancreas"
    MS = "ms"
    MYE = "mye"
    PANCREAS = "pancreas"
    
def get_data_loader(ds_name):
    if ds_name == SupportedDatasets.FILTERED_PANCREAS.value:
        return FilteredPancreas()
    if ds_name == SupportedDatasets.MS.value:
        return MS()
    if ds_name == SupportedDatasets.MYE.value:
        return Myeloid()
    if ds_name == SupportedDatasets.PANCREAS.value:
        return PancreaticDataset()
    raise ValueError("Invalid dataset name. Supported names are: ",
                     [ds.value for ds in SupportedDatasets])

class DataLoader(ABC):
    @abstractmethod
    def get_num_celltypes(self):
        pass

    @abstractmethod
    def get_train_test(self, test_ds_name=None):
        pass


class FilteredPancreas(DataLoader):
    DATASET_PATH = "data/datasets/annotation_pancreas/"
    TRAIN_DATA = "demo_train.h5ad"
    TEST_DATA = "demo_test.h5ad"

    def __init__(self):
        self.train_data = ad.read_h5ad(FilteredPancreas.DATASET_PATH + FilteredPancreas.TRAIN_DATA)
        self.train_data.obs.rename(columns={"Celltype": "celltype"}, inplace=True)
        self.test_data = ad.read_h5ad(FilteredPancreas.DATASET_PATH + FilteredPancreas.TEST_DATA)
        self.test_data.obs.rename(columns={"Celltype": "celltype"}, inplace=True)

    def get_num_celltypes(self):
        return len(self.train_data.obs["celltype"].unique())
    
    def get_train_test(self, test_ds_name=None):
        return self.train_data, self.test_data
    

class MS(DataLoader):
    DATASET_PATH = "data/datasets/ms/"
    TRAIN_DATA = "filtered_ms_adata.h5ad"
    TEST_DATA = "c_data.h5ad"

    def __init__(self):
        self.train_data = ad.read_h5ad(MS.DATASET_PATH + MS.TRAIN_DATA)
        self.test_data = ad.read_h5ad(MS.DATASET_PATH + MS.TEST_DATA)

    def get_num_celltypes(self):
        return len(self.train_data.obs["celltype"].unique())
    
    def get_train_test(self, test_ds_name=None):
        return self.train_data, self.test_data

class Myeloid(DataLoader):
    DATASET_PATH = "data/datasets/mye/"
    TRAIN_DATA = "reference_adata.h5ad"
    TEST_DATA = "query_adata.h5ad"

    def __init__(self):
        self.train_data = ad.read_h5ad(Myeloid.DATASET_PATH + Myeloid.TRAIN_DATA)
        self.train_data.obs.rename(columns={"cell_type": "celltype"}, inplace=True)
        self.test_data = ad.read_h5ad(Myeloid.DATASET_PATH + Myeloid.TEST_DATA)
        self.test_data.obs.rename(columns={"cell_type": "celltype"}, inplace=True)

    def get_num_celltypes(self):
        return len(self.train_data.obs["celltype"].unique())

    def get_train_test(self, test_ds_name=None):
        return self.train_data, self.test_data


class PancreaticDataset(DataLoader):

    DATASET_PATH = "../delta_tuning/data/scRNAseq_Benchmark_datasets/Intra-dataset/Pancreatic_data/"

    class SupportedDatasets(Enum):
        BARON = "Baron_Human"
        MURARO = "Muraro"
        SEGERSTOLPE = "Segerstolpe"
        XIN = "Xin"

    _DATA_FILE = "data/datasets/pancreas_all_genes/all_data.h5ad"

    def __init__(self):
        """
        Initialize the PancreaticDataset with the directory and filtering option.
        """

        self.adata = ad.read_h5ad(PancreaticDataset._DATA_FILE)
        self.adata_train = None
        self.adata_test = None

    
    def get_train_test(self, test_ds_name=None):
        """
        Split the dataset into training and testing sets based on the dataset name.

        Args:
            ds_name (str): The name of the dataset to split.

        Returns:
            adata_train (AnnData): The training set AnnData object.
            adata_test (AnnData): The testing set AnnData object.
        """

        available_names = [ds.value for ds in PancreaticDataset.SupportedDatasets]
        if test_ds_name not in available_names:
            raise ValueError("Name of the test dataset is not supported. Should be one of: ",
                             available_names)
                
        self.adata_test = self.adata[self.adata.obs["dataset"] == test_ds_name].copy()
        self.adata_train = self.adata[self.adata.obs["dataset"] != test_ds_name].copy()
        
        return self.adata_train, self.adata_test
    
    def get_num_celltypes(self):
        """
        Get the number of unique cell types in the dataset.

        Returns:
            int: The number of unique cell types.
        """
        if self.adata_train is None:
            warnings.warn("The train dataset is not specified. Returning the celltypes for the entire dataset, including test.")
            return len(self.adata.obs["celltype"].unique())
        return len(self.adata_train.obs["celltype"].unique())


## %%
if __name__ == "__main__":
    # Example usage
    dataloader = PancreaticDataset()
    adata, id2type = dataloader.load_all_datasets()

# %%
