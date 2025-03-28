# %%

import pandas as pd
import anndata as ad
import numpy as np
from pathlib import Path
from enum import Enum

class Katya:
    def __init__(self):
        print("Katya class initialized")

class PancreaticDataset:

    DATASET_PATH = "/local/home/ktubis/PycharmProjects/scGPT/delta_tuning/data/scRNAseq_Benchmark_datasets/Intra-dataset/Pancreatic_data/"

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

    def _get_curr_batch_id(self):
        """
        Get the current batch ID, and increase the counter by 1.
        """
        batch_id = str(self.curr_batch_id)
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
        
        # Add celltype_id column
        celltype_id_labels = adata.obs["celltype"].astype("category").cat.codes.values
        adata.obs["celltype_id"] = celltype_id_labels
        
        return adata
    

    def load_all_datasets(self):
        """
        Load and preprocess all supported datasets, and combine them into a single AnnData object.
        The datasets are filtered to keep only common genes across all datasets.
        Returns:
            combined_adata (AnnData): The combined AnnData object with common genes.
            combined_id2type (dict): Mapping from cell type IDs to cell type names for the combined dataset.
        """
        combined_adata = None
        adata_list = []

        for dataset in self.SupportedDatasets:
            # Load and preprocess the dataset
            adata  = self._preprocess_data(dataset)
            adata_list.append(adata)
        
        combined_adata = ad.concat(adata_list, join="inner")
        """
        # Filter combined_adata to keep only common genes
        if combined_adata is not None and common_genes is not None:
            combined_adata = combined_adata[:, list(common_genes)]"""

        # Create id2type mapping for the combined dataset
        combined_id2type = dict(enumerate(combined_adata.obs["celltype"].astype("category").cat.categories))

        return combined_adata, combined_id2type

## %%
if __name__ == "__main__":
    # Example usage
    dataloader = PancreaticDataset()
    adata, id2type = dataloader.load_all_datasets()

# %%
