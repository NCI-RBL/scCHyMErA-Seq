#importing libraries
import sys
import random
import os
import pandas as pd
import numpy as np
import torch
import scanpy as sc
import anndata as ad
from scipy.stats import median_abs_deviation
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import rcParams
FIGSIZE = (4, 4)
rcParams["figure.figsize"] = FIGSIZE

#setting seeds for reproducibility
random.seed(0)
np.random.seed(0)
torch.manual_seed(0)

#loading file
matrix_input = sys.argv[1]
##example path to matrix file
#matrix_input = "/data/kumarib2/all_single_cell/scanpy_seurat/AGG_main/outs/count"
adata = sc.read_10x_h5(os.path.join(matrix_input,"filtered_feature_bc_matrix.h5"),gex_only=True)
adata.var_names_make_unique()
##mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")
print(sum(adata.var['mt']))
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
adata.obs.rename({"n_genes_by_counts": "Genes Detected"}, axis=1, inplace=True)
adata.obs.rename({"total_counts": "UMI Counts"}, axis=1, inplace=True)
adata.obs.rename({"pct_counts_mt": "Mitochondrial Gene Percentage"}, axis=1, inplace=True)
##Plot genes detected
def plot_genes_detected(adata, output_dir):
    fig, ax = plt.subplots(dpi=300)
    ax = sc.pl.violin(adata, ['Genes Detected'], jitter=0.4, show=False)
    ax.tick_params(labelsize=12)
    ax.set_xticklabels([])
    plt.xticks([])
    ax.set_xlabel("")
    ax.set_ylabel("Number of genes", fontsize=14, family='Arial')
    ax.set_title("Genes expressed per cell", fontsize=14, family='Arial')
    plt.savefig(os.path.join(output_dir, "gene_count.pdf"), transparent=True, dpi=300, format='pdf')
    plt.close(fig)

#Plot UMI counts
def plot_umi_counts(adata, output_dir):
    fig, ax = plt.subplots(dpi=300)
    ax = sc.pl.violin(adata, ['UMI Counts'], jitter=0.4, show=False)
    ax.tick_params(labelsize=12)
    ax.set_xticklabels([])
    plt.xticks([])
    ax.set_ylabel("Number of UMIs", fontsize=14, family='Arial')
    ax.set_title("UMI counts per cell", fontsize=14, family='Arial')
    plt.savefig(os.path.join(output_dir, "umi_count.pdf"), transparent=True, dpi=300, format='pdf')
    plt.close(fig)

#Plot mitochondrial gene percentage
def plot_mito_pct(adata, output_dir):
    fig, ax = plt.subplots(dpi=300)
    ax = sc.pl.violin(adata, ['Mitochondrial Gene Percentage'], jitter=0.4, show=False)
   #ax.set_ylim(0,50)
    ax.tick_params(labelsize=12)
    ax.set_xticklabels([])
    plt.xticks([])
    ax.set_ylabel("% mitochondrial genes", fontsize=14, family='Arial')
    ax.set_title("Mitochondrial gene percentage per cell", fontsize=14, family='Arial')
    fig.set_size_inches(2, 2)
    plt.savefig(os.path.join(output_dir, "mitochondrial.pdf"), transparent=True, dpi=300, format='pdf')
    plt.close(fig)

#Correlation plot for total expressed genes vs UMI
def plot_correlation(adata, output_dir):
    fig, ax = plt.subplots(dpi=300)
    ax = sc.pl.scatter(adata, x='UMI Counts', y='Genes Detected', show=False)
    ax.tick_params(labelsize=12)
    ax.set_xlabel("UMI counts", family='Arial')
    ax.set_ylabel("Number of genes", fontsize=14, family='Arial')
    ax.set_title("Genes vs UMI counts per cell", fontsize=14, family='Arial')
    fig.set_size_inches(2, 2)
    plt.savefig(os.path.join(output_dir, "correlation.pdf"), transparent=True, dpi=300, format='pdf')
    plt.close(fig)
#
