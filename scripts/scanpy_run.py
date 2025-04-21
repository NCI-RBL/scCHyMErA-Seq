# importing libraries
import random
import os
import argparse
import sys
import pandas as pd
import numpy as np
from natsort import natsorted
import torch
import subprocess
import scanpy as sc
import pertpy as pt
import anndata as ad
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import rcParams
from scipy.stats import median_abs_deviation

# setting seed for reproducibility
random.seed(0)
np.random.seed(0)
torch.manual_seed(0)
torch.cuda.manual_seed_all(0)
#
######## get_args - start ##########################################################################

def get_args():
    """*get_args* - parses program's arg values.

    Returns
    arguments (*dict*) - user selected arguments parser.parse_args()

    """

    parser = argparse.ArgumentParser(
        prog='scanpy_analysis.py',
        description="Description: This script loads matrix file as scanpy object and performs processing and mixscape analysis.")

    #### Parameters
    parser.add_argument("-v","--version", action='version', version='%(prog)s version: Version 1.0.0 - Feb 2025')
    parser.add_argument("-o", "--out", help="Location of output directory where plots will be written.\nIf not specified, files will be written to the current working directory.", default='./', required=False)
    parser.add_argument("--analysis", choices=['Exon', 'KO'], help="gene or exon perturbation analysis.\nIf not specified, will do gene.", default='KO', required=False)
    parser.add_argument("--control", choices=['all', 'intergenic'], help="select either intergenic or intergenic+non-targeting as control.\nIf not specified, will consider only intergenic control.", default='intergenic', required=False)
    parser.add_argument("--resolution", help="Resolution for leiden clustering.\nIf not specified, will do 0.25.", default='0.25', required=False)
    # https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-h5-matrices

    group_required = parser.add_argument_group("required arguments")
    group_required.add_argument("-m", "--matrix_input", help="Path to matrix input in HDF5 Format (i.e filtered_feature_bc_matrix.h5). ", required=True)
    group_required.add_argument("-a", "--anno_csv", help="Path to guide annotation input.", required=True)
    group_required.add_argument("--timestamp", help="TimeStamp ", required=True)

    return parser.parse_args()

####### get_args - end ############################################################################

######## main - start ##############################################################################
# Details on the use of this script can be found in file scanpy_analysis.md. Interpretation of plots can be
# found:
#     https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html
#     https://pertpy.readthedocs.io/en/latest/tutorials/notebooks/mixscape.html
#     https://satijalab.org/seurat/articles/mixscape_vignette
#     https://www.sc-best-practices.org/cellular_structure/clustering.html
####################################################################################################

def main():
    """*main* - main function.

    """
    ######## Parse arguments ########
    args = get_args()

    #### matplotlib settings
    matplotlib.use('Agg')
    FIGSIZE = (3, 3)
    rcParams["figure.figsize"] = FIGSIZE

    #### Scanpy Settings
    sc.settings.verbosity = 0
    sc.settings.set_figure_params(
        dpi=300,
        fontsize=7,
        facecolor="white",
        frameon=True,
        transparent=True,
        vector_friendly=False
    )

# setting input variables
    my_wd=args.out
    matrix_input = args.matrix_input
    annotation_input = args.anno_csv
    analysis_type = args.analysis
    control_type = args.control
    resolution = float(args.resolution)
    log = open(my_wd+'scanpy_log.'+args.timestamp+'.txt', 'w+')   

    log.write("################################################################################" +'\n')
    log.write("############################# Running scanpy_analysis.py #############################" +'\n')
    log.write("Output Directory:  " + my_wd +'\n')
    log.write("Matrix Input File: " + matrix_input +'\n')
    log.write("Annotation Input File: " + annotation_input + '\n')
    
    ################################################################################
    ## Start analysis
    ################################################################################

    log.write("\n########################### Loading and Formatting Annotation Input File ##########################" +'\n')

    anno = pd.read_csv(annotation_input, sep="\t")

    # Selecting Gene knockout / Exons Perturbations

    if analysis_type == 'KO':
        anno = anno[anno['Cas9_Cas12a_targeted'].str.contains('_KO_|intergenic')]
    elif analysis_type == 'Exon':
        anno = anno[~anno['Cas9_Cas12a_targeted'].str.contains('_KO_')]
    else:
        log.write('Error : wrong analysis type' +'\n')

    mixscape_column = analysis_type

    if control_type == 'intergenic':
        anno = anno.loc[~anno['Cas9_Cas12a_targeted'].str.startswith('Non_Targeting')]
    else:
        anno['Cas9_Cas12a_targeted'] = anno['Cas9_Cas12a_targeted'].str.replace('Non_Targeting','intergenic')

    # Extract name of gene targeted
    split_interval = anno["Cas9_Cas12a_targeted"].str.split("_", expand=True)
    anno["KO"] = split_interval[0]
    anno["Exon_mod"] = anno["Cas9_Cas12a_targeted"].str.extract(r'(.*)_')
    anno["Exon"] = anno.loc[:, "Exon_mod"]
    anno["Exon"] = anno["Exon"].str.replace('_pos','')
    anno["Exon"] = anno["Exon"].str.replace('_neg','')
    anno["Replicate"] = anno["Cas9_Cas12a_targeted"].str.replace('.*_','',regex=True)
    anno['Perturbation'] = np.where(anno['Gene']=='intergenic', 'intergenic', 'Perturbed')
    anno['temp'] = anno.loc[:, 'cell_barcode']
    anno = anno.set_index('temp')
    anno.index.name = None
    
    log.write("\n########################### Loading Matrix Input File ##########################" +'\n')



    # Load matrix
    adata = sc.read_10x_h5(matrix_input, gex_only=True)
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith("MT-")

    # Handle any non-unique
    log.write("Handling: UserWarning by making variable names unique." +'\n')
    adata.var_names_make_unique()


    log.write("\n########################### Filtering and Normalizing Matrix ##########################" +'\n')

    sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    # Filter cells with less than 200 expressed genes and genes which were found in less than 3 cells.
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    # Get barcode from row names
    adata.obs['cell_barcode'] = adata.obs.index


    # Intersect Annotation and Matrix
    bdata = adata[adata.obs['cell_barcode'].isin(anno.cell_barcode)]
    bdata.obs = bdata.obs.merge(anno[['Cas9_Cas12a_targeted','num_features','KO','Exon','Replicate','Perturbation','cell_type']], how='left',left_index=True,right_index=True).copy()

    # save data
    bdata.write_h5ad("bdata." + analysis_type + ".h5ad")
    # bdata = sc.read_h5ad("bdata." + analysis_type + ".h5ad")
    #save raw counts in independent column
    bdata.layers['counts'] = bdata.X.copy()

    # Normalize to 10,000 reads per cell
    sc.pp.normalize_total(bdata, target_sum=1e4)
    # Log transform the data
    sc.pp.log1p(bdata)
    # Feature selection
    sc.pp.highly_variable_genes(bdata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pl.highly_variable_genes(bdata, save=analysis_type + '.genes_variability.pdf')
    bdata.raw = bdata
    bdata = bdata[:, bdata.var.highly_variable]

    # Regress out Total counts per Cell and Mitochondrial gene expression percentage
    sc.pp.regress_out(bdata, ["total_counts", "pct_counts_mt"])
    sc.pp.scale(bdata, max_value=10)
    # Save processed data
    bdata.write_h5ad("bdata." + analysis_type + ".scaled.h5ad")
    bdata.write_csvs("bdata." + analysis_type + ".scaled.cells")

    log.write("\n########################### Reducing Dimensions ##########################" +'\n')
    # Dimensionality reduction
    sc.tl.pca(bdata, svd_solver="arpack",random_state=0)
    sc.pl.pca_variance_ratio(bdata, log=True)
    # Compute distances in the PCA space, and find cell neighbors
    sc.pp.neighbors(bdata, n_neighbors=10, n_pcs=40)

    # Generate UMAP features
    sc.tl.umap(bdata)
    sc.pl.umap(bdata, color=["Perturbation"], use_raw=False, size=2, frameon=True, palette='rainbow', legend_fontsize='xx-small', save=analysis_type + '_origin.pdf')

    # Clustering
    sc.tl.leiden(bdata, resolution=resolution,random_state=0,flavor="igraph",n_iterations=2,directed=False)
    sc.pl.umap(bdata, color=['leiden'], legend_fontsize=8, legend_loc="on data", save=analysis_type + '_leiden.pdf')

    # Save bdata
    bdata.write_csvs(analysis_type + ".pca.cells")
    bdata.write_h5ad(analysis_type + ".pca.h5ad")
    bdata.raw.to_adata().write(analysis_type + ".pca_withoutX.h5ad")

    log.write("\n########################### Calculating local perturbation signatures ##########################" +'\n')

    # Run mixscape
    ms = pt.tl.Mixscape()
    ms.perturbation_signature(adata=bdata, pert_key=mixscape_column, control="intergenic", random_state=0)
    bdata_pert = bdata.copy()
    bdata_pert.X = bdata_pert.layers["X_pert"]

    # PCA
    sc.pp.pca(bdata_pert, svd_solver='arpack', random_state=0)
    sc.pp.neighbors(bdata_pert, metric="cosine", random_state=0)

    # UMAP
    sc.tl.umap(bdata_pert, random_state=0)
    sc.pl.umap(bdata_pert, color=["Perturbation"], palette= 'rainbow', frameon=True, legend_fontsize='xx-small', save=analysis_type + '.pertsig.pdf')

    # Clustering
    sc.tl.leiden(bdata_pert, resolution=resolution,random_state=0,flavor="igraph",n_iterations=2,directed=False)
    sc.pl.umap(bdata_pert, color=['leiden'], legend_fontsize=8, legend_loc="on data", save=analysis_type + '.leiden.pertsig.pdf')

    # Save bdata
    bdata_pert.write_csvs(analysis_type + ".pertsig.cells")
    bdata_pert.write_h5ad(analysis_type + ".pertsig.h5ad")

    log.write("\n########################### Running Mixscape ##########################" +'\n')

    # Run Mixscape
    ms.mixscape(adata=bdata_pert, control="intergenic", labels=mixscape_column)#, layer="X_pert")

    # Run LDA
    ms.lda(bdata, control=mixscape_column, labels="Exon", layer="X_pert")
    ms.plot_lda(bdata, control="intergenic", frameon=True, save=analysis_type + "_lda_plot.pdf", palette='rainbow')

    unique_id = [cls for cls in bdata.obs["mixscape_class"].unique() if cls.endswith("KO") or cls.startswith("intergenic")]
    unique_id_sorted = natsorted(unique_id)
    for feature in unique_id_sorted:
        palette_assign = [
                "red" if cls == feature else ("lightgrey",0.3)
                for cls in unique_id_sorted
                ]
        print(feature,palette_assign)
        fig = plt.figure(figsize=(8,8))
        ax = ms.plot_lda(bdata, control='intergenic', frameon=True, size=15, palette=palette_assign, show=False, legend_loc = None)
        plt.title(feature,fontsize=10)
        plt.savefig(f'{feature}.pdf', format='pdf')
        plt.close(fig)

    bdata_mixscape = bdata_pert.copy()
    
    bdata_mixscape.write_csvs(analysis_type + ".mixscape.cells")
    bdata_mixscape.write_h5ad(analysis_type + ".mixscape.h5ad")

    # Keep perturbed only cells
    bdata_perturb_only = bdata_mixscape[bdata_mixscape.obs["mixscape_class_global"] == "KO"].copy()

    # PCA
    sc.pp.pca(bdata_perturb_only,svd_solver='arpack',random_state=0)
    sc.pp.neighbors(bdata_perturb_only, metric="cosine",random_state=0)

    # UMAP
    sc.tl.umap(bdata_perturb_only,random_state=0)
    sc.pl.umap(bdata_perturb_only, color=["Perturbation"], palette= 'rainbow', frameon=True, legend_fontsize='xx-small', save=analysis_type + '.perturebed_only.pdf')

    # Clustering
    sc.tl.leiden(bdata_perturb_only, resolution=resolution,random_state=0,flavor="igraph",n_iterations=2,directed=False)
    sc.pl.umap(bdata_perturb_only, color=['leiden'], legend_fontsize=8, legend_loc="on data", save=analysis_type + '.leiden.perturebed_only.pdf')

    # Save bdata
    bdata_perturb_only.write_csvs(analysis_type + ".perturb_only_cells")
    bdata_perturb_only.write_h5ad(analysis_type + ".perturb_only.h5ad")
    #bdata_perturb_only = sc.read_h5ad(analysis_type + ".perturb_only.h5ad")


    log.write("\n########################### Computing Correlation Matrices ##########################" +'\n')

    #  Identify genes that are differentially expressed in the clusters 
    sc.tl.rank_genes_groups(bdata_perturb_only, groupby="leiden", method="wilcoxon")


    # Extract log-fold changes
    cluster_list = list(set(bdata_perturb_only.obs['leiden']))
    lfc_df = sc.get.rank_genes_groups_df(bdata_perturb_only, group = cluster_list)
    # Save as CSV
    lfc_df.to_csv(analysis_type + ".logfoldchanges.csv")

    # Visualize marker genes - show only the top 4 scoring gene
    sc.pl.rank_genes_groups_dotplot(
        bdata_perturb_only,
        n_genes=4,
        values_to_plot="logfoldchanges",
        min_logfoldchange=3,
        vmax=7,
        vmin=-7,
        cmap="bwr",
        save=analysis_type + '.rank_genes_groups_dotplot.pdf'
    )

    sc.pl.rank_genes_groups_matrixplot(
        bdata_perturb_only,
        n_genes=4,
        values_to_plot="logfoldchanges",
        vmin=-3,
        vmax=3,
        cmap="bwr",
        save=analysis_type + '.rank_genes_groups_matrixplot.pdf'
    )

    log.write("#Cell Files:\n")
    log.write(os.path.abspath("./*.perturb_only_cells") +'\n')

    log.write("\n# UMAP and Leiden clustering:" +'\n')
    log.write(os.path.abspath("./figures/umap*.pdf") +'\n')

    log.write("\n# Figures after calculating local perturbation signatures:" +'\n')
    log.write(os.path.abspath("./figures/umap*.pertsig.pdf") +'\n')

    log.write("\n# Figures using Perturbed cells only:" +'\n')
    log.write(os.path.abspath("./figures/umap*.perturebed_only.pdf")+'\n')

    log.write("\n# Dot Plots:" +'\n')
    log.write(os.path.abspath("./figures/*.rank_genes_groups_dotplot.pdf") +'\n')

    log.write("\n# Matrix Plots:" +'\n')
    log.write(os.path.abspath("./figures/*.rank_genes_groups_matrixplot.pdf") +'\n')

    log.write("\n# Correlation Matrix:" +'\n')
    log.write(os.path.abspath("./figures/*.correlation_matrix.pdf") +'\n')




if __name__ == "__main__":
        main()
    ######## main - end ################################################################################

