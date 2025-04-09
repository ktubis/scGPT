#!/bin/bash

python run_cell_annotation_model.py --model_config_path="model_configs/scgpt_pretrained_model.json" --model="../pretrained_models/best_model/best_model.pt" --results_file="results/pretrained_model_results.json" --log_file="logs/pretrained_model.txt" --test_max_seq_len=3001

