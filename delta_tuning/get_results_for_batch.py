import numpy as np
import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from sklearn.decomposition import PCA
import umap
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from sklearn.metrics import confusion_matrix, classification_report
import distinctipy
from scipy.stats import pearsonr
import sys
from load_ds import get_data_loader
import json

sys.path.insert(0, "../")
from scgpt.utils import eval_scib_metrics


def get_results(model_name, ds, test_adata, layer = 11):
    results = {}
    seeds = (101, 102, 103, 104, 105)
    for seed in seeds:
        output = ad.read_h5ad(f"predictions/batch/{model_name}_{ds}_{ds}_{seed}")
        print(f"predictions/batch/{model_name}_{ds}_{ds}_{seed}")
        adata = ad.AnnData(test_adata.X, obs=output.obs)
        adata.obsm['X_scGPT'] = output.obsm[f"transformer_layer_${layer}"]
        del output
        results[seed] = eval_scib_metrics(adata)
    arr_results = {}
    for key in results[seeds[0]]:
        arr_results[key] = [results[seed][key] for seed in seeds]
    return arr_results


def plot_umaps(model_name, ds, layer=11):
    output = ad.read_h5ad(f"predictions/batch/{model_name}_{ds}_{ds}_101")
    print(f"predictions/batch/{model_name}_{ds}_{ds}_101")
    adata = ad.AnnData(obs=output.obs)
    adata.obsm['X_scGPT'] = output.obsm[f"transformer_layer_${layer}"]
    del output

    sc.pp.neighbors(adata, use_rep="X_scGPT")
    sc.tl.umap(adata, min_dist=0.3)
    fig = sc.pl.umap(
        adata,
        color=["str_batch"],
        title= f"{model_name} {ds} batches",
        frameon=False,
        return_fig=True,
    )
    fig.savefig(f"plots/batch/umaps/{model_name}_{ds}_layer_{layer}_batch.png")

    fig = sc.pl.umap(
        adata,
        color=["celltype"],
        title=f"{model_name} {ds} celltype",
        frameon=False,
        return_fig=True,
    )
    fig.savefig(f"plots/batch/umaps/{model_name}_{ds}_layer_{layer}_celltype.png")


def main():
    datasets = ["covid", "cc_batch"]
    models = ["all_weights", "lora", "first_layer", "last_layer"]

    test_adatas = {}
    for ds in datasets:
        test_adatas[ds] = get_data_loader(ds).get_test_data()

    for ds in datasets:
        sc.pp.pca(test_adatas[ds], n_comps=512)
        sc.pp.neighbors(test_adatas[ds])
        sc.tl.umap(test_adatas[ds])
        fig = sc.pl.umap(
            test_adatas[ds],
            color=["celltype"],
            title=f"Original {ds} celltype",
            frameon=False,
            return_fig=True,
        )
        fig.savefig(f"plots/batch/umaps/Original_{ds}_celltype.png")

        fig = sc.pl.umap(
            test_adatas[ds],
            color=["str_batch"],
            title=f"Original {ds} batch",
            frameon=False,
            return_fig=True,
        )
        fig.savefig(f"plots/batch/umaps/Original_{ds}_batch.png")

    results = {}
    for ds in datasets:
        results[ds] = {}
        for model in models:
            results[ds][model] = get_results(model, ds, test_adatas[ds])
            plot_umaps(model, ds)

    with open("results/batch/covid_cc.json", "w") as f:
        json.dump(results, f)


    


if __name__ == "__main__":
    main()