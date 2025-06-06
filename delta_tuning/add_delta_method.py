from OpenDelta.opendelta import AdapterModel, LoraModel
from model_wrappers.cell_annotation import CellAnnotationModelWrapper

SUPPORTED_DELTA_TYPES = ["all_weights", "classifier", "specify_transformer_layers", "adapter", "lora"]

def add_open_delta_model(cam_model, delta_model_config, wandb_config):
    if delta_model_config["delta_method"] == "adapter":
        if "modified_modules" in delta_model_config:
            modified_modules = delta_model_config["modified_modules"]
        else:
            modified_modules = []
            for i in range(cam_model.nlayers):
                modified_modules.append(f"transformer_encoder.layers.{i}.linear2")
        delta_model = AdapterModel(backbone_model=cam_model, bottleneck_dim=delta_model_config["delta_param"],
                                   modified_modules=modified_modules)
        wandb_config["delta_param"] = delta_model_config["delta_param"]
    if delta_model_config["delta_method"] == "lora":
        modified_modules = []
        for i in range(cam_model.nlayers):
            modified_modules.append(f"transformer_encoder.layers.{i}.self_attn.Q")
            modified_modules.append(f"transformer_encoder.layers.{i}.self_attn.V")
        delta_model = LoraModel(backbone_model=cam_model, lora_r=delta_model_config["delta_param"], modified_modules=modified_modules)
        wandb_config["delta_param"] = delta_model_config["delta_param"]
    cam_model.to("cuda")
    delta_model.freeze_module(exclude=['deltas', 'cls_decoder'], set_state_dict=False)


def finetune_cls_decoder(cam: CellAnnotationModelWrapper):
    for param in cam.model.parameters():
        param.requires_grad = False
    
    for param in cam.model.cls_decoder.parameters():
        param.requires_grad = True

    trainable_params = sum(p.numel() for p in cam.model.parameters() if p.requires_grad) / sum(p.numel() for p in cam.model.parameters())
    cam.logger.info(f"Trainable parameters: {trainable_params}")


def train_transformer_modules(cam: CellAnnotationModelWrapper, delta_model_config: dict, wandb_config: dict):
    """
    Freeze all of the modules except the specified transformer modules in the config.
    """
    # train_modules: list, sublayers: list = None
    finetune_cls_decoder(cam)  # make sure cls decoder is trainable
    if not delta_model_config.get("train_modules", None):
        raise ValueError("train_modules must be specified in the delta model config for 'specify_transformer_layers' delta type.")
    
    train_modules = delta_model_config["train_modules"]
    sublayers = delta_model_config.get("sublayers", None)
    wandb_config["train_modules"] = train_modules
    wandb_config["sublayers"] = sublayers if sublayers is not None else []

    modules_to_train = ["transformer_encoder.layers." + str(i) for i in train_modules]
    if sublayers is not None:
        modules_to_train = [f"{module}.{sublayer}" for module in modules_to_train for sublayer in sublayers]

    for name, module in cam.model.named_modules():
        if name.startswith(tuple(modules_to_train)):
            for param in module.parameters():
                param.requires_grad = True
            cam.logger.info(f"Un-freezing module {name}")

    trainable_params = sum(p.numel() for p in cam.model.parameters() if p.requires_grad) / sum(p.numel() for p in cam.model.parameters())
    cam.logger.info(f"Trainable parameters after unfreezing transformers: {trainable_params}")


def add_delta_method(cam: CellAnnotationModelWrapper, delta_model_config: dict, wandb_config: dict):
    """
    Add the delta method to the model.
    :param cam: CellAnnotationModelWrapper
    :param model_config: model configuration
    :param wandb_config: wandb configuration
    """
    delta_method = delta_model_config["delta_method"]
    
    if delta_method == "all_weights":
        # nothing to do, all weights are trainable
        print("Delta type is 'all_weights', no additional configuration needed.")
        pass
    elif delta_method == "classifier":
        print("Delta type is 'classifier', finetuning the classifier decoder.")
        finetune_cls_decoder(cam)
    elif delta_method == "specify_transformer_layers":
        print("Delta type is 'specify_transformer_layers', training specified transformer modules.")
        train_transformer_modules(cam, delta_model_config, wandb_config)
    elif delta_method == "adapter" or delta_method == "lora":
        print(f"Delta type is '{delta_method}', adding OpenDelta model.")
        add_open_delta_model(cam.model, delta_model_config, wandb_config)
    else:
        raise ValueError(f"Unsupported delta type: {delta_method}. Supported types are: {SUPPORTED_DELTA_TYPES}")
    
    wandb_config["delta_method"] = delta_method