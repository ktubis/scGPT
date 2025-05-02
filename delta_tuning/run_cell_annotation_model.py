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
from datetime import datetime
from opendelta import AdapterModel, AutoDeltaConfig, AutoDeltaModel

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


def add_delta_model(model_config, cam_model):
        if model_config["delta_type"] == "adapter":
            modified_modules = []
            for i in range(cam_model.nlayers):
                modified_modules.append(f"transformer_encoder.layers.{i}.linear2")
                model_config["modified_modules"] = modified_modules
                delta_model = AdapterModel(backbone_model=cam_model, bottleneck_dim=model_config["bottleneck_dim"], modified_modules=modified_modules)
                cam_model.to("cuda")
                delta_model.freeze_module(exclude=['deltas', 'cls_decoder'])
        #delta_config = AutoDeltaConfig.from_dict(model_config)
        #delta_model = AutoDeltaModel.from_config(delta_config, backbone_model=cam_model)
        #cam_model.to("cuda")
        #delta_model.freeze_module(exclude=['deltas', 'cls_decoder'])


def main():
    parser = argparse.ArgumentParser(description="Cell Annotation Model")
    parser.add_argument("--model_config_path", type=str, default="model_configs/scgpt_pretrained_model.json", help="Path to the model config file")
    parser.add_argument("--model", type=str, default="../pretrained_models/best_model/best_model.pt", help="Path to the pretrained model to load. Must match the model config file. If None, will initialize a new model.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    parser.add_argument("--test_data", type=str, default="Muraro", help="Which dataset to use as the test data")
    parser.add_argument("--max_seq_len", type=int, default=301, help="Maximum sequence length")
    parser.add_argument("--epochs", type=int, default=1000, help="Number of epochs to train the model")
    parser.add_argument("--model_name", type=str, default="awesome_model", help="The name of the model to be saved")
    parser.add_argument("--finetune_decoder", action='store_true', help="Whether to finetune only the decoder.")
    parser.add_argument("--lr", type=float, default=1e-4)
    parser.add_argument("--wandb", action='store_true', help="Whether to use wandb for logging")
    parser.add_argument("--find_lr", action="store_true", help="Turns on the learning rate finder. The lr argument is the starting lr in that case, and it should be a negative power of 10.")
    parser.add_argument("--delta_configs_file", default=None, help="The file that stores the configs of the delta models. Must be a dict of dicts. If None, doesn't add a delta model.")
    parser.add_argument("--finetune_all_weights", action='store_true', help="Whether to finetune all the weights.")
    parser.add_argument("--inference", action='store_true', help="Whether to run inference on the given model.")
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
    adata_test2 = ds_loader.adata[ds_loader.adata.obs["batch_id"] == ds_loader.dataset_batch_dict["Xin"]].copy()

    # set up the preprocessor, use the args to config the workflow
    preprocessor = Preprocessor(
        use_key="X",
        normalize_total=0.0,
        binning=config_dict["n_input_bins"],
        result_binned_key=INPUT_LAYER,
    )

    adata_train = preprocess_data(adata_train, preprocessor, vocab)
    adata_test = preprocess_data(adata_test, preprocessor, vocab)
    adata_test2 = preprocess_data(adata_test2, preprocessor, vocab)

    num_celltypes = ds_loader.get_num_celltypes()
    num_batches = ds_loader.get_num_batches()

    # Get the current date and time, format as 'yymmddhhmm'
    formatted_time = datetime.now().strftime('%y%m%d%H%M')

    # model_name = args.model.split('/')[-1].split('.')[0]
    if not args.inference:
        model_name = args.model_name + f"_{formatted_time}"

    model_init_params = {
        "model_path": args.model,
        "max_seq_len": args.max_seq_len,
        "pad_value": config_dict["pad_value"],
        "vocab": vocab,
        "config_dict": config_dict,
        "num_batches": num_batches,
        "num_celltypes": num_celltypes,
        "model_name": model_name,
        "lr": args.lr,
        "wandb": args.wandb
    }

    if not args.inference:
        if args.delta_configs_file:
            with open(args.delta_configs_file, 'r') as f:
                configs = json.load(f)
            for model_name, config in configs.items():
                model_init_params["model_name"] = model_name
                cam = CellAnnotationModelWrapper(**model_init_params)
                add_delta_model(config, cam.model)
                print(cam.model)
                cam.train(args.epochs, adata_train, args.seed, adata_test1=adata_test, adata_test2=adata_test2, find_lr=args.find_lr)
        elif args.finetune_decoder:
            cam = CellAnnotationModelWrapper(**model_init_params)
            cam.finetune_cls_decoder()
            cam.train(args.epochs, adata_train, args.seed, adata_test1=adata_test, adata_test2=adata_test2, find_lr=args.find_lr)
        elif args.finetune_all_weights:
            cam = CellAnnotationModelWrapper(**model_init_params)
            cam.train(args.epochs, adata_train, args.seed, adata_test1=adata_test, adata_test2=adata_test2, find_lr=args.find_lr)
    
    _, _, results = cam.test(adata_test)
    
    results_file = Path("results/" + model_name)
    if not results_file.exists():
        results_file.touch()
    
    with open(results_file, "w") as f:
        json.dump(results, f)


if __name__ == "__main__":
    main()
# %%
