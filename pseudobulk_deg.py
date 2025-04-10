# importing libraries
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

######## get_args - start ##########################################################################

def get_args():
    """*get_args* - parses program's arg values.

    Returns
    arguments (*dict*) - user selected arguments parser.parse_args()

    """

    parser = argparse.ArgumentParser(
        prog='pseudobulk_deg.py',
        description="Description: This script generates csv files and plots showing differentially expressed genes for each perturbation.")

    #### Parameters
    parser.add_argument("-v","--version", action='version', version='%(prog)s version: Version 1.0.0 - Feb 2025')
    parser.add_argument("-o", "--out", help="Location of output directory where plots will be written.\nIf not specified, files will be written to the current working directory.", default='./', required=False)
    parser.add_argument("--analysis", help="KO or Exon analysis.\nIf not specified, will do KO.", default='KO', required=False)
    parser.add_argument("--control", help="select either intergenic or intergenic+non-targeting as control.\nIf not specified, will consider only intergenic control.", default='intergenic', required=False)

    # https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-h5-matrices

    group_required = parser.add_argument_group("required arguments")
    group_required.add_argument("-m", "--matrix_input", help="Path to matrix input in HDF5 Format (i.e filtered_feature_bc_matrix.h5). ", required=True)
    group_required.add_argument("-a", "--anno_csv", help="Path to guide annotation input.", required=True)
    group_required.add_argument("-p", "--mixscape_passed", help="Path to mixscape generated obs.csv file.", required=True)
    group_required.add_argument("--timestamp", help="TimeStamp ", required=True)

    return parser.parse_args()


######## get_args - end ############################################################################

######## main - start ##############################################################################
# Details on the use of this script can be found in file scanpy_analysis.md. Interpretation of plots can be
# found:
#     https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html
#     https://decoupler-py.readthedocs.io/en/latest/notebooks/pseudobulk.html
####################################################################################################

def main():
    """*main* - main function.

    """
    ######## Parse arguments ########
    args = get_args()

# setting seed for reproducibility
    random.seed(0)
    np.random.seed(0)
    torch.manual_seed(0)

# seting figure parameters
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

# setting input variables
    my_wd=args.out
    matrix_input = args.matrix_input
    annotation_input = args.anno_csv
    mixscape_pert = args.mixscape_passed
    log = open(my_wd+'log.'+args.timestamp+'.txt', 'w+')

    log.write("################################################################################" +'\n')
    log.write("############################# Running pseudobulk_deg.py #############################" +'\n')
    log.write("Output Directory:  " + my_wd +'\n')
    log.write("Matrix Input File: " + matrix_input +'\n')
    log.write("Annotation Input File: " + annotation_input + '\n')
    log.write("Mixscape generated Input file: " + mixscape_pert + '\n')

# loading metadata file with guide annotation for each cell
    anno = pd.read_csv(annotation_input, sep=",")

# selecting Gene knockout / Exons Perturbations
    if analysis_type == 'KO':
        anno = anno[anno['Cas9_Cas12a_targeted'].str.contains('_KO_|intergenic')]
    elif analysis_type == 'Exon':
        anno = anno[~anno['Cas9_Cas12a_targeted'].str.contains('_KO_')]
    else:
        log.write('Error : wrong analysis type' +'\n')

    anno = anno.loc[~anno['Cas9_Cas12a_targeted'].str.startswith('Non_Targeting')]

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

# filter cells with less than 200 expressed genes and genes which were found in less than 3 cells.
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.obs['cell_barcode'] = adata.obs.index
    bdata = adata[adata.obs['cell_barcode'].isin(anno.cell_barcode)]
    bdata.obs = bdata.obs.merge(anno[['Cas9_Cas12a_targeted','num_features','Gene','Exon','Replicate','Perturbation','cell_type']], how='left',left_index=True,right_index=True).copy()
   #print(bdata.obs)
   #print(bdata.var)

# save count data
    bdata.X = np.round(bdata.X)
    bdata.layers['counts'] = bdata.X.copy()

# normalize to 10,000 reads per cell
    sc.pp.normalize_total(bdata, target_sum=1e4)

# logarithmize and scale the data
    sc.pp.log1p(bdata)
    bdata.raw = bdata
    sc.pp.scale(bdata, max_value=10)

#dimensionality reduction
    sc.tl.pca(bdata, svd_solver="arpack",random_state=0)
    sc.pl.pca_variance_ratio(bdata, log=True)

# Compute distances in the PCA space, and find cell neighbors
    sc.pp.neighbors(bdata, n_neighbors=10, n_pcs=40)

##Generate UMAP features
    sc.tl.umap(bdata)
   #bdata.write_h5ad("scanpy_all.h5ad")
   #bdata.write_csvs("all_cells")

# select perturbed and control cells (from the file "obs.csv" after applying mixscape)
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

if __name__ == "__main__":
        main()
    ######## main - end ################################################################################

