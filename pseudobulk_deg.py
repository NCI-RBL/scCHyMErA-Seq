##importing packages
import random
import os
import decoupler as dc
import pandas as pd
import numpy as np
import torch
import scanpy as sc
import anndata as ad
from pydeseq2.dds import DeseqDataSet, DefaultInference
from pydeseq2.ds import DeseqStats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import median_abs_deviation
##setting seed for reproducibility
random.seed(0)
np.random.seed(0)
torch.manual_seed(0)
##seting figure parameters
sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=300,
    fontsize=7,
    figsize=(8,8),
    facecolor="white",
    frameon=True,
    transparent=True,
    vector_friendly=False
)
##setting input variables
matrix_input = args.matrix_input
annotation_input = args.anno_csv
mixscape_pert = args.mixscape_passed
filt = pd.read_csv(mixscape_pert, sep=",")
filt = filt[filt['mixscape_class_global'] != 'NP']
anno = pd.read_csv(annotation_input, sep=",")
anno = anno.loc[~anno['Cas9_Cas12a_targeted'].str.startswith('Non_Targeting')]
anno = anno[~anno['Cas9_Cas12a_targeted'].str.contains('_KO_')]
#anno = anno[anno['Cas9_Cas12a_targeted'].str.contains('_KO_|_intergenic_')]
split_interval = anno["Cas9_Cas12a_targeted"].str.split("_", expand=True)
anno["Gene"] = split_interval[0]
anno["Exon_mod"] = anno["Cas9_Cas12a_targeted"].str.extract(r'(.*)_')
anno["Exon"] = anno.loc[:, "Exon_mod"]
anno["Exon"] = anno["Exon"].str.replace('_pos','')
anno["Exon"] = anno["Exon"].str.replace('_neg','')
anno["Replicate"] = anno["Cas9_Cas12a_targeted"].str.replace('.*_','',regex=True)
anno['Perturbation'] = np.where(anno['Gene']=='intergenic', 'intergenic', 'Perturbed')
print(anno)
anno['temp'] = anno.loc[:, 'cell_barcode']
anno = anno.set_index('temp')
anno.index.name = None
adata = sc.read_10x_h5(matrix_input),gex_only=True)
adata.var_names_make_unique()
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
##filter cells with less than 200 expressed genes and genes which were found in less than 3 cells.
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.obs['cell_barcode'] = adata.obs.index
bdata = adata[adata.obs['cell_barcode'].isin(anno.cell_barcode)]
bdata.obs = bdata.obs.merge(anno[['Cas9_Cas12a_targeted','num_features','Gene','Exon','Replicate','Perturbation','cell_type']], how='left',left_index=True,right_index=True).copy()
#print(bdata.obs)
#print(bdata.var)
##save count data
bdata.X = np.round(bdata.X)
bdata.layers['counts'] = bdata.X.copy()
##normalize to 10,000 reads per cell
sc.pp.normalize_total(bdata, target_sum=1e4)
##logarithmize the data
sc.pp.log1p(bdata)
#bdata.layers['normalized'] = bdata.X
bdata.raw = bdata
sc.pp.scale(bdata, max_value=10)
##dimensionality reduction
sc.tl.pca(bdata, svd_solver="arpack",random_state=0)
sc.pl.pca_variance_ratio(bdata, log=True)
##Compute distances in the PCA space, and find cell neighbors
sc.pp.neighbors(bdata, n_neighbors=10, n_pcs=40)
##Generate UMAP features
sc.tl.umap(bdata)
#bdata.write_h5ad("scanpy_all.h5ad")
#bdata.write_csvs("all_cells")
##select perturbed cells and control cells (from the file "obs.csv" after applying mixscape)
filt = pd.read_csv(mixscape_pert, sep=",")
filt = filt[filt['mixscape_class_global'] != 'NP']
pert = bdata[bdata.obs['cell_barcode'].isin(filt.cell_barcode)]
pert.obs["Exon"] = pert.obs["Exon"].str.replace('_','range')
unique_exon = pert.obs[pert.obs['Exon'] != "intergenic"]['Exon'].unique().tolist()
# Do pseudobulk
pdata = dc.get_pseudobulk(
    pert,
    sample_col='Exon',
    groups_col='Replicate',
    layer='counts',
    mode='sum'
)
# Store raw counts in layers
pdata.layers['counts'] = pdata.X.copy()
# Normalize, scale and compute pca
sc.pp.normalize_total(pdata, target_sum=1e4)
sc.pp.log1p(pdata)
sc.pp.scale(pdata, max_value=10)
sc.tl.pca(pdata,svd_solver="arpack",random_state=0)
# Return raw counts to X
dc.swap_layer(pdata, 'counts', X_layer_key=None, inplace=True)
pdata.T.to_df().to_csv("matrix_file_exon")
pdata.write_h5ad("scanpy_pseudobulk.h5ad")
pdata.write_csvs("pseudobulk_cells")
inference = DefaultInference(n_cpus=8)
# Run PyDesq2 to determine differentially expressed genes
dds = DeseqDataSet(
    adata=pdata,
    design_factors='Exon',
    ref_level=['Exon', 'intergenic'],
    refit_cooks=True,
    inference=inference,
)
dds.deseq2()
for unique in unique_exon:
    stat_res= DeseqStats(dds,contrast=["Exon",unique,"intergenic"], inference=inference,)
    stat_res.summary()
    results_df = stat_res.results_df
    results_df.to_csv(f'{unique}.csv')
    dc.plot_volcano_df(results_df,x='log2FoldChange',y='padj',top=5,save=f'{unique}.pdf',figsize=(8, 4))

