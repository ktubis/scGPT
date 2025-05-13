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

def get_delta_config(method_name, delta_param):
    if method_name == "adapter":
        return {
            "delta_type": "adapter",
            "bottleneck_dim": delta_param,
        }
    if method_name == "soft_prompt":
        return {
            "delta_type": "soft_prompt",
            "soft_token_num": delta_param,
        }
    if method_name == "lora":
        return {
            "delta_type": "lora",
            "lora_r": delta_param,
        }

def add_delta_model(model_config, cam_model, wandb_config):
    if model_config["delta_type"] == "adapter":
        if "modified_modules" in model_config:
            modified_modules = model_config["modified_modules"]
        else:
            modified_modules = []
            for i in range(cam_model.nlayers):
                modified_modules.append(f"transformer_encoder.layers.{i}.linear2")
        delta_model = AdapterModel(backbone_model=cam_model, bottleneck_dim=model_config["bottleneck_dim"],
                                   modified_modules=modified_modules)
        cam_model.to("cuda")
        delta_model.freeze_module(exclude=['deltas', 'cls_decoder'])
        wandb_config["bottleneck_dim"] = model_config["bottleneck_dim"]
    if model_config["delta_type"] == "soft_prompt":
        delta_model = SoftPromptModel(backbone_model=cam_model, soft_token_num=model_config["soft_token_num"])
        wandb_config["soft_token_num"] = model_config["soft_token_num"]
    if model_config["delta_type"] == "lora":
        modified_modules = []
        for i in range(cam_model.nlayers):
            modified_modules.append(f"transformer_encoder.layers.{i}.self_attn.q_proj_weight")
            modified_modules.append(f"transformer_encoder.layers.{i}.self_attn.v_proj_weight")
        print("LEN MODIFIED MODULES:", len(modified_modules))
        delta_model = LoraModel(backbone_model=cam_model, lora_r=model_config["lora_r"], modified_modules=modified_modules)
        cam_model.to("cuda")
        delta_model.freeze_module(exclude=['deltas', 'cls_decoder'])
        wandb_config["lora_r"] = model_config["lora_r"]

    trainable_params = sum(p.numel() for p in cam_model.parameters() if p.requires_grad) / sum(p.numel() for p in cam_model.parameters())
    print(f"Trainable parameters: {trainable_params}")
    #delta_config = AutoDeltaConfig.from_dict(model_config)
    #delta_model = AutoDeltaModel.from_config(delta_config, backbone_model=cam_model)
    #cam_model.to("cuda")
    #delta_model.freeze_module(exclude=['deltas', 'cls_decoder'])


def get_data_loaders(ds_loader, vocab, config_dict, test_data):
    adata_train, adata_test = ds_loader.get_train_test(test_data)
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

    return adata_train, adata_test, adata_test2

def add_tokens_to_vocab(pad_token, vocab):
    special_tokens = [pad_token, "<cls>", "<eoc>"]
    for s in special_tokens:
        if s not in vocab:
            vocab.append_token(s)


def hyperparameter_search(cam, num_epochs, adata_train, adata_test, trial):
    split_data, gene_ids = cam._get_split_data(adata_train)
    optimizer = torch.optim.Adam(
        cam.model.parameters(), lr=cam.lr, eps=cam.eps
    )
    move_optimizer_params_to_cuda(optimizer)
    scheduler = torch.optim.lr_scheduler.StepLR(
            optimizer, cam.schedule_interval, gamma=cam.schedule_ratio
        )
    scaler = torch.cuda.amp.GradScaler(enabled=True)
    prev_epoch_loss = np.inf
    bad_epochs_counter = 0
    assert num_epochs > 0, "num_epochs must be greater than 0"
    for epoch in range(1, num_epochs + 1):
        train_loader, valid_loader = cam._get_train_valid_data_per_epoch(split_data, gene_ids)
        epoch_loss = cam._train_step(train_loader, optimizer, scheduler, scaler, epoch)
        print("epoch loss:", epoch_loss)
        _, _, test_results = cam.test(adata_test)
        f1_score = test_results["test/macro_f1"]
        print("f1 score:", f1_score)
        trial.report(f1_score, step=epoch)

        # prune if doesn't improve on train in the last past epochs
        if epoch_loss >= prev_epoch_loss - (OPTUNA_PRUNING_PERCENTAGE * prev_epoch_loss) / 100:
            bad_epochs_counter += 1
        else:
            bad_epochs_counter = 0
        if bad_epochs_counter >= OPTUNA_PRUNING_EPOCHS or trial.should_prune():
            raise optuna.exceptions.TrialPruned()
        scheduler.step()
    # return the loss of the last epoch
    return f1_score
        

def update_model_config(model_config, hyperparams):
    model_config["lr"] = hyperparams["lr"]
    model_config["schedule_interval"] = hyperparams["schedule_interval"]
    model_config["schedule_ratio"] = hyperparams["schedule_ratio"]
    return model_config

def find_hyperparams(model_init_params, grid_search_config, num_epochs, adata_train, adata_test, delta_method, wandb_config, log_wandb=False):
    def optuna_objective(trial):
        print("trial number:", trial.number)
        with open(grid_search_config, "r") as f:
            grid_search_dict = json.load(f)
        min_lr = grid_search_dict["min_lr"]
        max_lr = grid_search_dict["max_lr"]
        min_delta_param = grid_search_dict["min_delta_param"]
        max_delta_param = grid_search_dict["max_delta_param"]
        min_schaduler_interval = grid_search_dict["min_scheduler_interval"]
        max_scheduler_interval = grid_search_dict["max_scheduler_interval"]
        min_scheduler_ratio = grid_search_dict["min_scheduler_ratio"]
        max_scheduler_ratio = grid_search_dict["max_scheduler_ratio"]

        lr = trial.suggest_loguniform("lr", min_lr, max_lr)
        schedule_interval = trial.suggest_int("schedule_interval", min_schaduler_interval, max_scheduler_interval)
        schedule_ratio = trial.suggest_float("schedule_ratio", min_scheduler_ratio, max_scheduler_ratio)
        delta_param = trial.suggest_int("delta_param", min_delta_param, max_delta_param)
        hyperparams = {
            "lr": lr,
            "schedule_interval": schedule_interval,
            "schedule_ratio": schedule_ratio,
        }

        print(hyperparams)
        print("delta_param:", delta_param)

        wandb_config["lr"] = lr
        wandb_config["schedule_interval"] = schedule_interval
        wandb_config["schedule_ratio"] = schedule_ratio
        wandb_config["delta_param"] = delta_param
        wandb_config["model_name"] = f"{wandb_config['model_name']}_{trial.number}"
        init_wandb(wandb_config)
        log_wandb = True

        # Update the model config with the hyperparameters and create model
        update_model_config(model_init_params, hyperparams)
        cam = CellAnnotationModelWrapper(**model_init_params)

        if delta_method == "all_weights":
            pass
        elif delta_method == "finetune_classifier":
            cam.finetune_cls_decoder()
        else:
            delta_config = get_delta_config(delta_method, delta_param)
            add_delta_model(delta_config, cam.model, wandb_config)

        if log_wandb:
            wandb.watch(cam.model, log="all", log_graph=True)

        return hyperparameter_search(cam, num_epochs, adata_train, adata_test, trial)
    return optuna_objective


def train_model(args, adata_train, adata_test, cam, wandb_config, adata_test2=None):
    if args.delta_configs_file:
        with open(args.delta_configs_file, 'r') as f:
            delta_config = json.load(f)
        print(delta_config)
        add_delta_model(delta_config, cam.model, wandb_config)
        print(cam.model)
    elif args.finetune_decoder:
        cam.finetune_cls_decoder()
        wandb_config["delta_method"] = "finetune_classifier"
    else:
        assert (args.finetune_all_weights == True), "If you want to finetune all weights, please set the finetune_decoder flag to True."
        wandb_config["delta_method"] = "all_weights"
    if args.wandb:
        init_wandb(wandb_config)
        wandb.watch(cam.model, log="all", log_graph=True)
    cam.train(args.epochs, adata_train, args.seed, adata_test1=adata_test, adata_test2=adata_test2, find_lr=args.find_lr)


def main():
    parser = argparse.ArgumentParser(description="Cell Annotation Model")
    parser.add_argument("--model_config_path", type=str, default="model_configs/scgpt_pretrained_model.json", help="Path to the model config file")
    parser.add_argument("--model", type=str, default=None, help="Path to the pretrained model to load. Must match the model config file. If None, will initialize a new model.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    parser.add_argument("--test_data", type=str, default="Muraro", help="Which dataset to use as the test data")
    parser.add_argument("--max_seq_len", type=int, default=3001, help="Maximum sequence length")
    parser.add_argument("--epochs", type=int, default=1000, help="Number of epochs to train the model")
    parser.add_argument("--model_name", type=str, default="awesome_model", help="The name of the model to be saved")
    parser.add_argument("--finetune_decoder", action='store_true', help="Whether to finetune only the decoder.")
    parser.add_argument("--lr", type=float, default=1e-4)
    parser.add_argument("--wandb", action='store_true', help="Whether to use wandb for logging")
    parser.add_argument("--find_lr", action="store_true", help="Turns on the learning rate finder. The lr argument is the starting lr in that case, and it should be a negative power of 10.")
    parser.add_argument("--delta_configs_file", default=None, help="The file that stores the configs of the delta models. Must be a dict of dicts. If None, doesn't add a delta model.")
    parser.add_argument("--finetune_all_weights", action='store_true', help="Whether to finetune all the weights.")
    parser.add_argument("--inference", action='store_true', help="Whether to run inference on the given model.")
    parser.add_argument("--schedule_interval", type=int, default=20, help="The interval at which to schedule the learning rate.")
    parser.add_argument("--schedule_ratio", type=float, default=0.9, help="The ratio of the learning rate to schedule.")
    parser.add_argument("--hyperparam_search_config", default=None, help="The configuration from which to do hyperparameter search.")
    args = parser.parse_args()

    set_seed(args.seed)
    with open(args.model_config_path, "r") as f:
        config_dict = json.load(f)

    vocab = GeneVocab.from_file(VOCAB_PATH)
    vocab = GeneVocab.from_file(VOCAB_PATH)
    add_tokens_to_vocab(config_dict["pad_token"], vocab)

    ds_loader = load_ds.PancreaticDataset()
    adata_train, adata_test, adata_test2 = get_data_loaders(ds_loader, vocab, config_dict, args.test_data)

    num_celltypes = ds_loader.get_num_celltypes()
    num_batches = ds_loader.get_num_batches()

    # Get the current date and time, format as 'yymmddhhmm'
    formatted_time = datetime.now().strftime('%y%m%d%H%M')

    # model_name = args.model.split('/')[-1].split('.')[0]
    if args.inference:
        model_name = args.model.split('/')[-1].split('.')[0]
    else:
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
        "wandb": args.wandb,
        "schedule_interval": args.schedule_interval,
        "schedule_ratio": args.schedule_ratio
    }

    wandb_config = {
        "learning_rate": args.lr,
        "model_name": model_name,
        "epochs": args.epochs,
        "schedule_ratio": args.schedule_ratio,
        "schedule_interval": args.schedule_interval,
        "seed": args.seed,
    }

    if args.hyperparam_search_config:
        if args.finetune_all_weights:
            delta_method = "all_weights"
        elif args.finetune_decoder:
            delta_method = "finetune_classifier"
        elif args.delta_configs_file:
            with open(args.delta_configs_file, 'r') as f:
                delta_config = json.load(f)
            delta_method = delta_config["delta_type"]
        else:
            raise ValueError("Please provide a delta method for hyperparameter search.")
        
        print("Running hyperparameter search...")
        study = optuna.create_study(direction="maximize")
        objective = find_hyperparams(model_init_params, args.hyperparam_search_config, args.epochs, adata_train, adata_test,
                                        delta_method, wandb_config, log_wandb=args.wandb)
        study.optimize(objective, n_trials=50)
        best_trial = study.best_trial
        print("Best trial:")
        print("  Value: {}".format(best_trial.value))
        print("  Params: ")
        for key, value in best_trial.params.items():
            print("    {}: {}".format(key, value))

    else:
        cam = CellAnnotationModelWrapper(**model_init_params)
        if not args.inference:
            train_model(args, adata_train, adata_test, cam, wandb_config, adata_test2)

        _, _, results = cam.test(adata_test)
        
        results_file = Path("results/" + model_name)
        if not results_file.exists():
            results_file.touch()
        
        with open(results_file, "w") as f:
            json.dump(results, f)


if __name__ == "__main__":
    main()
# %%
