# %%
from abc import ABC, abstractmethod
import pandas as pd
import anndata as ad
import numpy as np
from pathlib import Path
from enum import Enum
import warnings
from scipy.sparse import issparse


class SupportedDatasets(Enum):
    FILTERED_PANCREAS = "pc"
    MS = "ms"
    MYE = "mye"
    PANCREAS = "pancreas"
    PANCREAS_OLD = "pancreas_old"  # For backward compatibility, if needed
    SWAPPED_PANCREAS = "swapped_pc"
    COLORECTAL_CANCER = "cc"
    HIGH_CORR_CC = "high_cc"
    WEAK_CORR_CC = "weak_cc"
    WEAK_CORR_CC2 = "weak_cc2"
    INTESTINE = "int"
    DOWNSAMPLED_INT = "downsampled_int"
    DOWNSAMPLED_CC = "downsampled_cc"
    PBMC = "pbmc"
    COVID = "covid"

    
def get_data_loader(ds_name):
    print("DS NAME", ds_name)
    if ds_name == SupportedDatasets.FILTERED_PANCREAS.value:
        return FilteredPancreas()
    if ds_name == SupportedDatasets.MS.value:
        return MS()
    if ds_name == SupportedDatasets.MYE.value:
        return Myeloid()
    if ds_name == SupportedDatasets.PANCREAS.value:
        return PancreaticDataset()
    if ds_name == SupportedDatasets.PANCREAS_OLD.value:
        return PancreaticDatasetOLD()
    if ds_name == SupportedDatasets.SWAPPED_PANCREAS.value:
        return SwappedPancreas()
    if ds_name == SupportedDatasets.COLORECTAL_CANCER.value:
        return ColorectalCancer()
    if ds_name == SupportedDatasets.HIGH_CORR_CC.value:
        return HighCorrColorectalCancer()
    if ds_name == SupportedDatasets.WEAK_CORR_CC.value:
        return WeakCorrColorectalCancer()
    if ds_name == SupportedDatasets.WEAK_CORR_CC2.value:
        return WeakCorrColorectalCancer2()
    if ds_name == SupportedDatasets.INTESTINE.value:
        return Intestine()
    if ds_name == SupportedDatasets.PBMC.value:
        return PBMC()
    if ds_name == SupportedDatasets.COVID.value:
        return Covid()
    if ds_name == SupportedDatasets.DOWNSAMPLED_INT.value:
        return DownsampledInt()
    if ds_name == SupportedDatasets.DOWNSAMPLED_CC.value:
        return DownsampledCC()
    raise ValueError("Invalid dataset name. Supported names are: ",
                     [ds.value for ds in SupportedDatasets])

class DataLoader(ABC):

    def __init__(self):
        self.train_data = None
        self.test_data = None
        self.num_celltypes = 0
        self.train_path = ""
        self.test_path = ""

    def read_data(self):
        self.train_data = ad.read_h5ad(self.train_path)
        self.test_data = ad.read_h5ad(self.test_path)

    def get_num_batches(self):
        """
        Get the number of unique batches in the dataset.
        Returns:
            int: The number of unique batches.
        """
        if "batch_id" in self.train_data.obs:
            return len(self.train_data.obs["batch_id"].unique())
        else:
            raise ValueError("Batch information is not available in the dataset. "
                             "Ensure that the dataset contains a 'batch_id' column in obs.")
        
    def add_celltype_id(self):
        """
        Add a celltype_id column to the AnnData object.
        The celltype_id is an integer code representing the cell type.
        """
        if "celltype" not in self.train_data.obs or "celltype" not in self.test_data.obs:
            raise ValueError("Cell type information is not available in the dataset. "
                             "Ensure that the dataset contains a 'celltype' column in obs.")
        
        celltypes = set(self.train_data.obs["celltype"].values) | set(self.test_data.obs["celltype"].values)
        self.num_celltypes = len(celltypes)
        print("CELLTYPES:", celltypes)
        celltype_to_id = {ct: i for i, ct in enumerate(sorted(celltypes))}
        self.train_data.obs["celltype_id"] = self.train_data.obs["celltype"].map(celltype_to_id)
        print("Celltype to ID mapping:", self.train_data.obs["celltype_id"].unique())
        self.test_data.obs["celltype_id"] = self.test_data.obs["celltype"].map(celltype_to_id)
        print("Celltype to ID mapping test:", self.test_data.obs["celltype_id"].unique())

    def get_num_celltypes(self):
        return self.num_celltypes
    
    def get_train_test(self, test_ds_name=None):
        return self.train_data, self.test_data
        

class FilteredPancreas(DataLoader):
    def __init__(self):
        self.train_path = "data/datasets/annotation_pancreas/demo_train.h5ad"
        self.test_path = "data/datasets/annotation_pancreas/demo_test.h5ad"
        self.read_data()
        self.train_data.obs.rename(columns={"Celltype": "celltype"}, inplace=True)
        self.test_data.obs.rename(columns={"Celltype": "celltype"}, inplace=True)
        self.add_celltype_id()

class SwappedPancreas(DataLoader):
    def __init__(self):
        self.train_path = "data/datasets/pancreas_swapped/train.h5ad"
        self.test_path = "data/datasets/pancreas_swapped/test.h5ad"
        self.read_data()
        self.train_data.obs.rename(columns={"Celltype": "celltype"}, inplace=True)
        self.test_data.obs.rename(columns={"Celltype": "celltype"}, inplace=True)
        self.add_celltype_id()

    

class MS(DataLoader):

    @staticmethod
    def preprocess_data(adata):
        """
        Preprocess the AnnData object by renaming columns and ensuring categorical types.
        """
        adata.obs.rename(columns={"str_batch": "batch_id"}, inplace=True)
        adata.obs["batch_id"] = adata.obs["batch_id"].astype("int")

    def __init__(self):
        self.train_path = "data/datasets/ms/filtered_ms_adata.h5ad"
        self.test_path = "data/datasets/ms/c_data.h5ad"
        self.read_data()
        self.__class__.preprocess_data(self.train_data)
        self.__class__.preprocess_data(self.test_data)
        self.add_celltype_id()


class Myeloid(DataLoader):
    @staticmethod
    def preprocess_data(adata):
        """
        Preprocess the AnnData object by renaming columns and ensuring categorical types.
        """
        adata.obs.rename(columns={"cell_type": "celltype"}, inplace=True)
        adata.obs.rename(columns={"batch": "batch_id"}, inplace=True)
        adata.obs["batch_id"] = adata.obs["batch_id"].astype("int")

    def __init__(self):
        self.train_path = "data/datasets/mye/reference_adata.h5ad"
        self.test_path = "data/datasets/mye/query_adata.h5ad"
        self.read_data()
        self.__class__.preprocess_data(self.train_data)
        self.__class__.preprocess_data(self.test_data)
        self.add_celltype_id()

class ColorectalCancer(DataLoader):
    @staticmethod
    def preprocess_data(adata):
        """
        Preprocess the AnnData object by renaming columns and ensuring categorical types.
        """
        adata.obs["batch_id"] = adata.obs["batchID"].astype("category").cat.codes
        adata.obs.rename(columns={"celltype": "celltype_small"}, inplace=True)
        adata.obs.rename(columns={"celltype_mid": "celltype"}, inplace=True)

    def __init__(self):
        self.train_path = "data/datasets/10x_v2_v3_immune_hubs/train_all_cells_3k.h5ad"
        self.test_path = "data/datasets/10x_v2_v3_immune_hubs/test_all_cells_3k.h5ad"
        self.read_data()
        self.__class__.preprocess_data(self.train_data)
        self.__class__.preprocess_data(self.test_data)
        self.add_celltype_id()

class DownsampledCC(DataLoader):
    @staticmethod
    def preprocess_data(adata):
        """
        Preprocess the AnnData object by renaming columns and ensuring categorical types.
        """
        adata.obs["batch_id"] = adata.obs["batchID"].astype("category").cat.codes
        adata.obs.rename(columns={"celltype": "celltype_small"}, inplace=True)
        adata.obs.rename(columns={"celltype_mid": "celltype"}, inplace=True)

    def __init__(self):
        self.train_path = "data/datasets/10x_v2_v3_immune_hubs/downsampled_train.h5ad"
        self.test_path = "data/datasets/10x_v2_v3_immune_hubs/test_all_cells_3k.h5ad"
        self.read_data()
        self.__class__.preprocess_data(self.train_data)
        self.__class__.preprocess_data(self.test_data)
        self.add_celltype_id()


class WeakCorrColorectalCancer(DataLoader):
    @staticmethod
    def preprocess_data(adata):
        """
        Preprocess the AnnData object by renaming columns and ensuring categorical types.
        """
        adata.obs["batch_id"] = adata.obs["batch_id"].astype("int")

    def __init__(self):
        self.train_path = "data/datasets/10x_v2_v3_immune_hubs/weak_corr_train.h5ad"
        self.test_path = "data/datasets/10x_v2_v3_immune_hubs/weak_corr_test.h5ad"
        self.read_data()
        self.__class__.preprocess_data(self.train_data)
        self.__class__.preprocess_data(self.test_data)
        self.add_celltype_id()

class WeakCorrColorectalCancer2(DataLoader):
    @staticmethod
    def preprocess_data(adata):
        """
        Preprocess the AnnData object by renaming columns and ensuring categorical types.
        """
        adata.obs["batch_id"] = adata.obs["batch_id"].astype("int")

    def __init__(self):
        self.train_path = "data/datasets/10x_v2_v3_immune_hubs/weak_corr2_train.h5ad"
        self.test_path = "data/datasets/10x_v2_v3_immune_hubs/weak_corr2_test.h5ad"
        self.read_data()
        self.__class__.preprocess_data(self.train_data)
        self.__class__.preprocess_data(self.test_data)
        self.add_celltype_id()

class HighCorrColorectalCancer(DataLoader):
    @staticmethod
    def preprocess_data(adata):
        """
        Preprocess the AnnData object by renaming columns and ensuring categorical types.
        """
        adata.obs["batch_id"] = adata.obs["batchID"].astype("category").cat.codes
        adata.obs.rename(columns={"celltype": "celltype_small"}, inplace=True)
        adata.obs.rename(columns={"celltype_mid": "celltype"}, inplace=True)

    def __init__(self):
        self.train_path = "data/datasets/10x_v2_v3_immune_hubs/hard_train.h5ad"
        self.test_path = "data/datasets/10x_v2_v3_immune_hubs/hard_test.h5ad"
        self.read_data()
        self.__class__.preprocess_data(self.train_data)
        self.__class__.preprocess_data(self.test_data)
        self.add_celltype_id()


class Intestine(DataLoader):
    @staticmethod
    def preprocess_data(adata):
        """
        Preprocess the AnnData object by renaming columns and ensuring categorical types.
        """
        adata.obs["batch_id"] = adata.obs["batch"].astype("int")
        adata.obs.rename(columns={"cell_type": "celltype"}, inplace=True)

    def __init__(self):
        self.train_path = "data/datasets/intestine/train_not_beumer_not_holloway.h5ad"
        self.test_path = "data/datasets/intestine/test_beumer.h5ad"
        self.read_data()
        self.__class__.preprocess_data(self.train_data)
        self.__class__.preprocess_data(self.test_data)
        self.add_celltype_id()

class DownsampledInt(DataLoader):
    @staticmethod
    def preprocess_data(adata):
        """
        Preprocess the AnnData object by renaming columns and ensuring categorical types.
        """
        adata.obs["batch_id"] = adata.obs["batch"].astype("int")
        adata.obs.rename(columns={"cell_type": "celltype"}, inplace=True)

    def __init__(self):
        self.train_path = "data/datasets/intestine/downsampled_train.h5ad"
        self.test_path = "data/datasets/intestine/test_beumer.h5ad"
        self.read_data()
        self.__class__.preprocess_data(self.train_data)
        self.__class__.preprocess_data(self.test_data)
        self.add_celltype_id()

class PBMC(DataLoader):
    @staticmethod
    def preprocess_data(adata):
        """
        Preprocess the AnnData object by renaming columns and ensuring categorical types.
        """
        adata.obs["batch_id"] = adata.obs["batch"].astype("category").cat.codes
        adata.obs["str_batch"] = adata.obs["batch"].astype(str)
        adata.obs.rename(columns={"str_labels": "celltype"}, inplace=True)

    def __init__(self):
        self.train_path = "data/batch_correction/pbmc.h5ad"
        self.train_data = ad.read_h5ad(self.train_path)
        self.__class__.preprocess_data(self.train_data)
        self.test_data = self.train_data[self.train_data.obs["batch_id"].argsort()].copy()
        self.add_celltype_id()

class Covid(DataLoader):
    @staticmethod
    def preprocess_data(adata):
        """
        Preprocess the AnnData object by renaming columns and ensuring categorical types.
        """
        adata.obs["batch_id"] = adata.obs["study2"].astype("category").cat.codes
        adata.obs["str_batch"] = adata.obs["study2"].astype(str)
        adata.obs.rename(columns={"str_labels": "celltype"}, inplace=True)

    def __init__(self):
        self.train_path = "data/batch_correction/covid19.h5ad"
        self.train_data = ad.read_h5ad(self.train_path)
        self.__class__.preprocess_data(self.train_data)
        self.test_data = self.train_data[self.train_data.obs["batch_id"].argsort()].copy()
        self.add_celltype_id()


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
        self.adata.obs.rename(columns={"ds_on_concat": "batch_id"}, inplace=True)
        self.adata.obs["batch_id"] = self.adata.obs["batch_id"].astype("int")
        self.train_data = None
        self.test_data = None


    def get_train_test(self, test_ds_name=None):
        """
        Split the dataset into training and testing sets based on the dataset name.

        Args:
            ds_name (str): The name of the dataset to split.

        Returns:
            train_data (AnnData): The training set AnnData object.
            test_data (AnnData): The testing set AnnData object.
        """

        available_names = [ds.value for ds in PancreaticDataset.SupportedDatasets]
        print("test_ds_name:", test_ds_name)
        if test_ds_name not in available_names:
            raise ValueError("Name of the test dataset is not supported. Should be one of: ",
                             available_names)
                
        self.test_data = self.adata[self.adata.obs["dataset"] == test_ds_name].copy()
        self.train_data = self.adata[self.adata.obs["dataset"] != test_ds_name].copy()
        self.add_celltype_id()
        
        return self.train_data, self.test_data
    

class PancreaticDatasetOLD():

    DATASET_PATH = "../delta_tuning/data/scRNAseq_Benchmark_datasets/Intra-dataset/Pancreatic_data/"

    class SupportedDatasets(Enum):
        BARON = "Baron_Human"
        MURARO = "Muraro"
        SEGERSTOLPE = "Segerstolpe"
        XIN = "Xin"

    _DATA_FILE = "Filtered_data.csv"
    _LABEL_FILE = "Labels.csv"

    def __init__(self):
        """
        Initialize the PancreaticDataset with the directory and filtering option.
        """
        self.data_dir = Path(PancreaticDataset.DATASET_PATH)
        self.curr_batch_id = 0
        self.dataset_batch_dict = {}    # a dict that maps dataset name to batch id
        self.adata = None
        self.id2type = None
        self._load_all_datasets()

    def _get_curr_batch_id(self):
        """
        Get the current batch ID, and increase the counter by 1.
        """
        batch_id = self.curr_batch_id
        self.curr_batch_id += 1
        return batch_id


    def _preprocess_data(self, dataset):
        """
        Preprocess data from a folder containing 'Filtered_data.csv' and 'Labels.csv'.
        Args:
            folder_path (str or Path): Path to the folder containing the data files.
        Returns:
            adata (AnnData): The preprocessed AnnData object.
            id2type (dict): Mapping from cell type IDs to cell type names.
        """
        dataset_name = dataset.value
        folder_path = self.data_dir / dataset_name

        # Load gene counts and labels
        gene_counts = pd.read_csv(folder_path / self._DATA_FILE, index_col=0)
        if "Unnamed: 0" in gene_counts.columns:
            gene_counts = gene_counts.rename(columns={"Unnamed: 0": "gene_name"})
        labels = pd.read_csv(folder_path / self._LABEL_FILE, index_col=0)

        if dataset == self.SupportedDatasets.MURARO:
            genes = [gene.split('__')[0] for gene in gene_counts.columns]
            gene_counts.columns = genes
            labels.rename(index={"duct": "ductal"}, inplace=True)

        # Filter duplicate genes
        gene_counts = gene_counts.loc[:, ~gene_counts.columns.duplicated()]

        # Create AnnData object
        adata = ad.AnnData(X=gene_counts)

        # Add cell annotations
        adata.obs["celltype"] = labels.index

        # Ensure celltype is categorical
        if "celltype" in adata.obs:
            adata.obs["celltype"] = adata.obs["celltype"].astype("category")

        # Add batch information
        adata.obs["batch_id"] = self._get_curr_batch_id()  # Default batch ID

        # Set gene names as index
        adata.var["gene_name"] = gene_counts.columns
        adata.var.set_index("gene_name", inplace=True)


        return adata


    def _load_all_datasets(self):
        """
        Load and preprocess all supported datasets, and combine them into a single AnnData object.
        The datasets are filtered to keep only common genes across all datasets.
        Returns:
            combined_adata (AnnData): The combined AnnData object with common genes.
            combined_id2type (dict): Mapping from cell type IDs to cell type names for the combined dataset.
        """
        adata_list = []

        for dataset in self.SupportedDatasets:
            self.dataset_batch_dict[dataset.value] = self.curr_batch_id
            adata  = self._preprocess_data(dataset)
            adata_list.append(adata)

        self.adata = ad.concat(adata_list, join="inner")

        # Add celltype_id column
        celltype_id_labels = self.adata.obs["celltype"].astype("category").cat.codes.values
        self.adata.obs["celltype_id"] = celltype_id_labels

        # Create id2type mapping for the combined dataset
        self.id2type = dict(enumerate(self.adata.obs["celltype"].astype("category").cat.categories))

    def get_train_test(self, test_ds_name="Xin"):
        """
        Split the dataset into training and testing sets based on the dataset name.
        Args:
            ds_name (str): The name of the dataset to split.
        Returns:
            adata_train (AnnData): The training set AnnData object.
            adata_test (AnnData): The testing set AnnData object.
        """

        if test_ds_name not in self.dataset_batch_dict:
            raise ValueError(f"Dataset {test_ds_name} not found in the dataset batch dictionary.")
        
        batch_id = self.dataset_batch_dict[test_ds_name]
        # Split the dataset into training and testing sets
        adata_test = self.adata[self.adata.obs["batch_id"] == batch_id].copy()
        adata_train = self.adata[self.adata.obs["batch_id"] != batch_id].copy()

        return adata_train, adata_test
    
    def get_num_celltypes(self):
        """
        Get the number of unique cell types in the dataset.
        Returns:
            int: The number of unique cell types.
        """
        return len(self.adata.obs["celltype"].cat.categories)


    def get_num_batches(self):
        """
        Get the number of unique batches in the dataset.
        Returns:
            int: The number of unique batches.
        """
        return len(self.adata.obs["batch_id"].unique())


## %%
if __name__ == "__main__":
    # Example usage
    dataloader = PancreaticDataset()
    adata, id2type = dataloader.load_all_datasets()

# %%
