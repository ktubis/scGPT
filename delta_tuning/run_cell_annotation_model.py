#%%


#%load_ext autoreload
#%autoreload 2


import json
import sys
import load_ds
from model_wrappers.cell_annotation import CellAnnotationModelWrapper, move_optimizer_params_to_cuda
import random
import torch
import numpy as np
import os
import argparse
from pathlib import Path
from datetime import datetime
import optuna
import wandb
from OpenDelta.opendelta import AdapterModel, SoftPromptModel, LoraModel
from transformers import get_linear_schedule_with_warmup
import logging
import add_delta_method

sys.path.insert(0, "../")
from scgpt.model import TransformerModel
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.preprocess import Preprocessor


#TODO: make a model attributes file for the cell annotation model
INPUT_LAYER = "X_binned"
VOCAB_PATH = "../pretrained_models/best_model/vocab.json"
SEED = 42
FINETUNED_MOODELS_PATH = "retrained_models/"
OPTUNA_PRUNING_EPOCHS = 5
OPTUNA_PRUNING_PERCENTAGE = 0.05
HYPERPARAMS_SEARCH_DIR = "hyperparam_search_results/"
OPTUNA_LOGS_DIR = HYPERPARAMS_SEARCH_DIR + "optuna_logs/"


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


def init_wandb(wandb_config):
    wandb.init(
        config=wandb_config,
        project="cell_annotation",
        reinit=True,
        settings=wandb.Settings(start_method="fork"),
        mode='offline',
        name=wandb_config["model_name"],
    )

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


def get_data_loaders(ds_loader, vocab, config_dict, test_data):
    adata_train, adata_test = ds_loader.get_train_test(test_data)

    # set up the preprocessor, use the args to config the workflow
    preprocessor = Preprocessor(
        use_key="X",
        normalize_total=0.0,
        binning=config_dict["n_input_bins"],
        result_binned_key=INPUT_LAYER,
    )

    adata_train = preprocess_data(adata_train, preprocessor, vocab)
    adata_test = preprocess_data(adata_test, preprocessor, vocab)

    return adata_train, adata_test

def add_tokens_to_vocab(pad_token, vocab):
    special_tokens = [pad_token, "<cls>", "<eoc>"]
    for s in special_tokens:
        if s not in vocab:
            vocab.append_token(s)


def hyperparameter_search(cam, num_epochs, adata_train, adata_test, trial, warm_up_percentage, batch_size, logger):
    split_data, gene_ids = cam._get_split_data(adata_train)
    optimizer = torch.optim.Adam(
        cam.model.parameters(), lr=cam.lr, eps=cam.eps
    )
    move_optimizer_params_to_cuda(optimizer)
    num_training_steps = num_epochs * np.ceil(len(split_data.train_data) / batch_size)
    warm_up_steps = (warm_up_percentage / 100.) * num_training_steps
    scheduler = get_linear_schedule_with_warmup(
                optimizer,
                num_warmup_steps=warm_up_steps,
                num_training_steps=num_training_steps,
            )
    scaler = torch.cuda.amp.GradScaler(enabled=True)
    bad_epochs_counter = 0
    assert num_epochs > 0, "num_epochs must be greater than 0"
    best_eval_loss = np.inf
    prev_eval_loss = np.inf

    for epoch in range(1, num_epochs + 1):
        train_loader, valid_loader = cam._get_train_valid_data_per_epoch(split_data, gene_ids)
        epoch_loss = cam._train_step(train_loader, optimizer, scheduler, scaler, epoch)
        eval_loss, _ = cam._evaluate(valid_loader)
        if eval_loss < best_eval_loss:
            best_eval_loss = eval_loss
        if eval_loss > prev_eval_loss:
            bad_epochs_counter += 1
        prev_eval_loss = eval_loss
        _, _, test_results = cam.test(adata_test)
        f1_score = test_results["test/macro_f1"]
        trial.report(eval_loss, step=epoch)
        logger.info(
            f"Epoch {epoch:03d} | "
            f"train_loss: {epoch_loss:.4f} | "
            f"eval_loss: {eval_loss:.4f} | "
            f"best_eval: {best_eval_loss:.4f} | "
            f"f1: {f1_score:.4f} | "
            f"lr: {scheduler.get_last_lr()[0]:.4e}"
        )

        if cam.log_wandb:
            wandb.log(
                {
                    "epoch": epoch,
                    "train_loss": epoch_loss,
                    "eval_loss": eval_loss,
                    "best_eval_loss": best_eval_loss,
                    "f1_score": f1_score,
                    "lr": scheduler.get_last_lr()[0],
                }
            )
                
        # Return if the model doesn't improve for a certain number of epochs
        if bad_epochs_counter >= OPTUNA_PRUNING_EPOCHS:
            logger.info(f"Early stopping at epoch {epoch} due to no improvement in eval loss, with best eval loss: {best_eval_loss}")
            return best_eval_loss

        is_in_warmup_phase = epoch <= warm_up_percentage * num_epochs / 100
        if not is_in_warmup_phase and trial.should_prune():
            logger.info(f"Trial {trial.number} pruned at epoch {epoch} with eval loss: {eval_loss:.4f}")
            raise optuna.exceptions.TrialPruned()

    # return the eval loss of the best epoch
    return best_eval_loss
        

def update_model_config(model_config, hyperparams):
    model_config["lr"] = hyperparams["lr"]
    model_config["schedule_interval"] = hyperparams["schedule_interval"]
    model_config["schedule_ratio"] = hyperparams["schedule_ratio"]
    return model_config

def find_hyperparams(model_init_params, adata_train, adata_test, delta_config, hyperparam_search_config, wandb_config, logger):
    def optuna_objective(trial):
        min_lr = hyperparam_search_config["min_lr"]
        max_lr = hyperparam_search_config["max_lr"]
        lr = trial.suggest_float("lr", min_lr, max_lr, log=True)
        num_epochs = trial.suggest_categorical("epochs", hyperparam_search_config["epochs"])
        warm_up_percentage = trial.suggest_categorical("warm_up_percentage", hyperparam_search_config["warm_up_percentages"])
        if hyperparam_search_config.get("delta_params", None) is not None:
            delta_param = trial.suggest_categorical("delta_param", hyperparam_search_config["delta_params"])
            delta_config["delta_param"] = delta_param
        else:
            delta_param = -1

        delta_method = delta_config["delta_method"]
        logger.info(f"Trial {trial.number}: lr: {lr}, delta_param: {delta_param}, num_epochs: {num_epochs}, warm_up_percentage: {warm_up_percentage}, delta_method: {delta_method}")

        wandb_config["lr"] = lr
        wandb_config["num_epochs"] = num_epochs
        wandb_config["warm_up_percentage"] = warm_up_percentage
        wandb_config["model_name"] = f"{wandb_config['model_name']}_{trial.number}"
        init_wandb(wandb_config)

        # Update the model config with the hyperparameters and create model
        model_init_params["lr"] = lr
        cam = CellAnnotationModelWrapper(**model_init_params)

        add_delta_method.add_delta_method(cam, delta_config, wandb_config)

        wandb.watch(cam.model, log="all", log_graph=True)
        return hyperparameter_search(cam, num_epochs, adata_train, adata_test, trial,
                                     warm_up_percentage=warm_up_percentage,
                                     batch_size=model_init_params["batch_size"],
                                     logger=logger)
    return optuna_objective

def run_optuna_hyperparam_search(model_init_params, adata_train, adata_test, delta_config,
                                 wandb_config, hyperparam_search_config_file, logger):
        os.makedirs(OPTUNA_LOGS_DIR, exist_ok=True)
        name_log_file = f"optuna_trials_{wandb_config['model_name']}_delta_method_{delta_config['delta_method']}.log"
        log_path = os.path.join(OPTUNA_LOGS_DIR, name_log_file)
        # Set the Optuna logging level to INFO
        optuna.logging.set_verbosity(optuna.logging.INFO)
        logger = optuna.logging.get_logger("optuna")
        optuna.logging.disable_default_handler()
        logger.setLevel(logging.INFO)

        # Prevent double logging
        #if not logger.hasHandlers():
        file_handler = logging.FileHandler(log_path)
        formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        
        print("Running hyperparameter search...")
        study = optuna.create_study(direction="minimize")
        with open(hyperparam_search_config_file, 'r') as f:
            hyperparam_search_config = json.load(f)
        objective = find_hyperparams(model_init_params, adata_train, adata_test, delta_config,
                                     hyperparam_search_config, wandb_config, logger)
        n_trials = hyperparam_search_config["n_trials"]
        study.optimize(objective, n_trials=n_trials)
        best_trial = study.best_trial
        logger.info("Best trial:")
        logger.info("  Value: {}".format(best_trial.value))
        logger.info("  Params: ")
        for key, value in best_trial.params.items():
            logger.info("    {}: {}".format(key, value))
        df = study.trials_dataframe()

        os.makedirs(HYPERPARAMS_SEARCH_DIR, exist_ok=True)
        df.to_csv(f"{HYPERPARAMS_SEARCH_DIR}/{wandb_config['model_name']}.csv", index=False)


def train_model(args, adata_train, adata_test, cam, delta_config, wandb_config, warm_up_percentage=0):
    print(delta_config)
    add_delta_method.add_delta_method(cam, delta_config, wandb_config)
    if args.wandb:
        init_wandb(wandb_config)
        wandb.watch(cam.model, log="all", log_graph=True)
    cam.train(args.epochs, adata_train, args.seed, adata_test=adata_test, find_lr=args.find_lr,
              warm_up_percentage=warm_up_percentage, early_stop=args.early_stop)


def main():
    parser = argparse.ArgumentParser(description="Cell Annotation Model")
    parser.add_argument("--model_config_path", type=str, default="model_configs/scgpt_pretrained_model.json", help="Path to the model config file")
    parser.add_argument("--model", type=str, default=None, help="Path to the pretrained model to load. Must match the model config file. If None, will initialize a new model.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    parser.add_argument("--test_data", type=str, default="Muraro", help="Which dataset to use as the test data, only for the Pancreatic dataset. Must be one of: Baron_Human, Muraro, Segerstolpe, Xin")
    parser.add_argument("--train_data", type=str, default="pancreas", help="Which dataset to use for the train data.")
    parser.add_argument("--max_seq_len", type=int, default=3001, help="Maximum sequence length")
    parser.add_argument("--epochs", type=int, default=1000, help="Number of epochs to train the model")
    parser.add_argument("--model_name", type=str, default="awesome_model", help="The name of the model to be saved")
    parser.add_argument("--lr", type=float, default=1e-4)
    parser.add_argument("--wandb", action='store_true', help="Whether to use wandb for logging")
    parser.add_argument("--find_lr", action="store_true", help="Turns on the learning rate finder. The lr argument is the starting lr in that case, and it should be a negative power of 10.")
    parser.add_argument("--delta_config_file", help="The file that stores the configs of the delta models. Must be a dict of dicts. If None, doesn't add a delta model.")
    parser.add_argument("--inference", action='store_true', help="Whether to run inference on the given model.")
    parser.add_argument("--schedule_interval", type=int, default=20, help="The interval at which to schedule the learning rate.")
    parser.add_argument("--schedule_ratio", type=float, default=0.9, help="The ratio of the learning rate to schedule.")
    parser.add_argument("--hyperparam_search_config", default=None, help="The configuration from which to do hyperparameter search.")
    parser.add_argument("--batch_size", type=int, default=32, help="The batch size to use for training.")
    parser.add_argument("--warm_up_percentage", type=int, default=0, help="The percentage of epochs to dedicate to the warm up.")
    parser.add_argument("--early_stop", action='store_true', help="Whether to use early stopping.")
    args = parser.parse_args()

    supported_datasets = [ds.value for ds in load_ds.SupportedDatasets]
    if args.train_data not in supported_datasets:
        raise ValueError("Name of the train dataset is not supported. Should be one of: ",
                         supported_datasets)
    set_seed(args.seed)
    with open(args.model_config_path, "r") as f:
        config_dict = json.load(f)

    vocab = GeneVocab.from_file(VOCAB_PATH)
    vocab = GeneVocab.from_file(VOCAB_PATH)
    add_tokens_to_vocab(config_dict["pad_token"], vocab)
    print("train data:", args.train_data in supported_datasets)
    ds_loader = load_ds.get_data_loader(args.train_data)
    adata_train, adata_test = get_data_loaders(ds_loader, vocab, config_dict, args.test_data)

    num_celltypes = ds_loader.get_num_celltypes()
    num_batches = ds_loader.get_num_batches()

    # Get the current date and time, format as 'yymmddhhmm'
    formatted_time = datetime.now().strftime('%y%m%d%H%M')

    # model_name = args.model.split('/')[-1].split('.')[0]
    if args.inference:
        model_name = args.model.split('/')[-1].split('.')[0]
    else:
        model_name = args.model_name + f"_{formatted_time}"

    if not args.inference:
        with open(args.delta_config_file, 'r') as f:
            delta_config = json.load(f)

    wandb_config = {
        "learning_rate": args.lr,
        "model_name": model_name,
        "epochs": args.epochs,
        "seed": args.seed,
        "log_wandb": args.wandb,
        "dataset": args.train_data,
    }

    model_init_params = {
        "model_path": args.model,
        "max_seq_len": args.max_seq_len,
        "pad_value": config_dict["pad_value"],
        "vocab": vocab,
        "config_dict": config_dict,
        "num_celltypes": num_celltypes,
        "num_batches": num_batches,
        "model_name": model_name,
        "lr": args.lr,
        "log_wandb": args.wandb,
        "schedule_interval": args.schedule_interval,
        "schedule_ratio": args.schedule_ratio,
        "batch_size": args.batch_size,
        "eval_batch_size": args.batch_size,
    }

    if args.hyperparam_search_config:
        assert not args.inference, "Hyperparameter search cannot be run in inference mode."
        run_optuna_hyperparam_search(model_init_params, adata_train, adata_test, delta_config,
                                     wandb_config, args.hyperparam_search_config, logging.getLogger())
    else:
        cam = CellAnnotationModelWrapper(**model_init_params)
        if not args.inference:
            train_model(args, adata_train, adata_test, cam, delta_config, wandb_config, warm_up_percentage=args.warm_up_percentage)
            test_kwargs = {}
        else:
            predictions_file = f"predictions/{model_name}_{args.train_data}_{args.seed}"
            test_kwargs = {
                "predictions_file": predictions_file if args.inference else None,
                "save_embeddings": True,
            }

        _, _, results = cam.test(adata_test, **test_kwargs)
        print("Results:", results)
        
        results_file = Path("results/" + model_name)
        if not results_file.exists():
            results_file.touch()
        
        with open(results_file, "w") as f:
            json.dump(results, f)


if __name__ == "__main__":
    main()
# %%
