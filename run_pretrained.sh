#!/bin/bash

python delta_tuning/run_cell_annotation_model.py --model_config_path="delta_tuning/model_configs/scgpt_pretrained_model.json" --model="pretrained_models/best_model/best_model.pt" --results_file="delta_tuning/results/pretrained_model_results.json" --log_file="delta_tuning/logs/pretrained_model.txt"

