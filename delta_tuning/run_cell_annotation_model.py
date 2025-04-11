#%%


#%load_ext autoreload
#%autoreload 2


import json
import sys
import load_ds
from model_wrappers.cell_annotation import CellAnnotationModelWrapper
import random
import torch
import numpy as np
import os
import argparse
from pathlib import Path

sys.path.insert(0, "../")
from scgpt.model import TransformerModel
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.preprocess import Preprocessor

#TODO: make a model attributes file for the cell annotation model
INPUT_LAYER = "X_binned"
VOCAB_PATH = "../pretrained_models/best_model/vocab.json"
SEED = 42


def set_seed(seed):
    os.environ['PYTHONHASHSEED'] = str(seed)
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)  # if using multiple GPUs
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False


def preprocess_data(adata, preprocessor, vocab):
    """
    Filter genes in the AnnData object based on their presence in the vocabulary.
    
    Args:
        adata (AnnData): The AnnData object containing gene expression data.
    
    Returns:
        adata (AnnData): The filtered AnnData object.
    """
    preprocessor(adata, batch_key=None)
    adata.var["id_in_vocab"] = [
        1 if gene in vocab else -1 for gene in adata.var.index
    ]
    return adata[:, adata.var["id_in_vocab"] >= 0]


def main():
    parser = argparse.ArgumentParser(description="Cell Annotation Model")
    parser.add_argument("--model_config_path", type=str, required=True, help="Path to the model config file")
    parser.add_argument("--model", type=str, default=None, help="Path to the pretrained model to load. Must match the model config file. If None, will initialize a new model.")
    parser.add_argument("--results_file", type=str, required=True, help="Path to the file in which to save the results")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    parser.add_argument("--test_data", type=str, default="Muraro", help="Which dataset to use as the test data")
    parser.add_argument("--test_batch_size", type=int, default=32, help="Batch size for testing")
    parser.add_argument("--max_seq_len", type=int, default=301, help="Maximum sequence length")
    parser.add_argument("--log_file", type=str, default="log.txt", help="Name of the log file")
    parser.add_argument("--train", type=bool, default=False, help="Whether to train the model or not")
    parser.add_argument("--epochs", type=int, default=10, help="Number of epochs to train the model")
    parser.add_argument("--model_name", type=str, default="awesome_model", help="The name of the model to be saved")
    args = parser.parse_args()

    set_seed(args.seed)

    with open(args.model_config_path, "r") as f:
        config_dict = json.load(f)

    vocab = GeneVocab.from_file(VOCAB_PATH)
    special_tokens = [config_dict["pad_token"], "<cls>", "<eoc>"]
    for s in special_tokens:
        if s not in vocab:
            vocab.append_token(s)

    ds_loader = load_ds.PancreaticDataset()
    adata_train, adata_test = ds_loader.get_train_test(args.test_data)

    # set up the preprocessor, use the args to config the workflow
    preprocessor = Preprocessor(
        use_key="X",
        normalize_total=0.0,
        binning=config_dict["n_input_bins"],
        result_binned_key=INPUT_LAYER,
    )

    adata_train = preprocess_data(adata_train, preprocessor, vocab)
    adata_test = preprocess_data(adata_test, preprocessor, vocab)

    num_celltypes = ds_loader.get_num_celltypes()
    num_batches = ds_loader.get_num_batches()

    cam = CellAnnotationModelWrapper(
        model_path=args.model,
        max_seq_len=args.max_seq_len,
        pad_value=config_dict["pad_value"],
        vocab=vocab,
        config_dict=config_dict,
        num_batches=num_batches,
        num_celltypes=num_celltypes,
        model_name=args.model_name,
    )

    print(cam.model)

    if args.train:
        cam.train(args.epochs, adata_train)
    
    predictions, celltypes_labels, results = cam.test(adata_test, eval_batch_size=args.test_batch_size)
    
    results_file = Path(args.results_file)
    if not results_file.exists():
        results_file.touch()
        print("File created successfully")        
    
    with open(args.results_file, "w") as f:
        json.dump(results, f)

    log_file = Path(args.log_file)
    if not log_file.exists():
        log_file.touch()
        print("Log file created successfully")

    with open(args.log_file, "w") as f:
        f.write(f"Predictions:\n {list(predictions)}\n")
        f.write(f"Celltype Labels:\n {list(celltypes_labels)}\n")


if __name__ == "__main__":
    main()