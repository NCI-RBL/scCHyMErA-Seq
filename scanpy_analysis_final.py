## Python Imports for both QC and Analysis
import os
#import decoupler as dc
import pandas as pd
import numpy as np
import anndata as ad
from pydeseq2.dds import DeseqDataSet, DefaultInference
from pydeseq2.ds import DeseqStats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import rcParams
from scipy.stats import median_abs_deviation
import argparse
import anndata as ad
import subprocess
import pertpy as pt
import scanpy as sc
import seaborn as sns


######## get_args - start ##########################################################################

def get_args():
    """*get_args* - parses program's arg values.

    Returns
    arguments (*dict*) - user selected arguments parser.parse_args()

    """

    parser = argparse.ArgumentParser(
        prog='scanpy_analysis.py',
        description="Description: This script 'scanpy_analysis.py' generates runs the scanpy + mixscape analysis.\nSee included documentation file 'scanpy_analysis.md' for details")

    #### Parameters
    parser.add_argument("-v","--version", action='version', version='%(prog)s version: Version 1.0.0 - Feb 2025')
    parser.add_argument("-o", "--out", help="Location of output directory where plots will be written.\nIf not specified, files will be written to the current working directory.", default='./', required=False)
    parser.add_argument("--analysis", help="KO or Exon analysis.\nIf not specified, will do KO.", default='KO', required=False)
    parser.add_argument("--resolution", help="Resolution for leiden clustering.\nIf not specified, will do 0.25.", default='0.25', required=False)

    # https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-h5-matrices

    group_required = parser.add_argument_group("required arguments")
    group_required.add_argument("-m", "--matrix_input", help="Path to matrix input HDF5 Format (i.e filtered_feature_bc_matrix.h5). ", required=True)
    group_required.add_argument("-a", "--anno_csv", help="Path to matrix input directory. ", required=True)
    group_required.add_argument("--timestamp", help="TimeStamp ", required=True)

    return parser.parse_args()


######## get_args - end ############################################################################

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
        facecolor="white",
        frameon=True,
        vector_friendly=False
    )

    # Set seed
    np.random.seed(0)

    # Set Input Variables
    my_wd=args.out
    matrix_input = args.matrix_input
    annotation_input = args.anno_csv
    analysis_type = args.analysis
    resolution = float(args.resolution)
    log = open(my_wd+'scanpy_log.'+args.timestamp+'.txt', 'w+')

    print("################################################################################")
    print("############################# Running scanpy_analysis.py #############################")
    print("################################################################################")
    print("Output Directory:  " + my_wd)
    print("Matrix Input File: " + matrix_input)
    print("Annotation Input File: " + annotation_input)


    ################################################################################
    ## Start analysis
    ################################################################################


    print("\n########################### Loading and Formatting Annotation Input File ##########################")

    #anno = pd.read_csv("select_pairs_1noise.csv", sep="\t")
    anno = pd.read_csv(annotation_input, sep="\t")


    # Separating Genes KO / Exons Perturbations and intergenic
    # Also removing Non_Targeting

    if analysis_type == 'KO':
        anno = anno[anno['Cas9_Cas12a_targeted'].str.contains('_KO_|intergenic')]
    elif analysis_type == 'Exon':
        anno = anno[anno['Cas9_Cas12a_targeted'].str.contains('_pos_|_neg_|intergenic')]
    else:
        print('Error : wrong analysis type')

    # Extract name of gene targeted
    split_interval = anno["Cas9_Cas12a_targeted"].str.split("_", expand=True)
    anno["Gene"] = split_interval[0]

    # Extract Exon / Exon_mod
    if analysis_type == 'KO':
        anno["Exon_mod"] = anno["Cas9_Cas12a_targeted"].str.extract(r'(.*)_')
        anno["Exon"] = anno.loc[:, "Exon_mod"]
        anno["Exon"] = anno["Exon"].str.replace('_KO','')
    elif analysis_type == 'Exon':
        anno["Exon_mod"] = anno["Cas9_Cas12a_targeted"].str.extract(r'([^_]*)_') # extract everything before first underscore
        anno["Exon"] = anno.loc[:, "Exon_mod"]

    # Extract Replicate
    anno["Replicate"] = anno["Cas9_Cas12a_targeted"].str.replace('.*_','',regex=True)

    # Extract Pertubation status 
    anno['Perturbation'] = np.where(anno['Exon']=='intergenic', 'Intergenic', 'Perturbed')

    # Set cell type
    anno['cell_type'] = 'single'

    # Make barcodes row names
    anno['temp'] = anno.loc[:, 'cell_barcode']
    anno = anno.set_index('temp')
    anno.index.name = None

    # Pilot vs Full Project
    anno['Project'] = anno.loc[:, 'cell_barcode']
    split_interval = anno["Project"].str.split("-", expand=True)
    anno["Sample"] = split_interval[1]
    anno["Project"] = np.where(anno["Sample"].isin(['9', '10']), "Pilot", "Full")

    # Keep Full Only
    anno = anno[anno["Project"] == "Full"]


    print("\n########################### Loading Matrix Input File ##########################")

    # Load matrix
    #matrix_input = "/mnt/gridftp/guibletwm/CCBRRBL13/AGG_main_toy/outs/count"
    #adata = sc.read_10x_h5(os.path.join(matrix_input,"filtered_feature_bc_matrix.h5"),gex_only=True)
    adata = sc.read_10x_h5(matrix_input, gex_only=True)
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith("MT-")

    # Handle any non-unique
    print("Handling: UserWarning by making variable names unique.")
    adata.var_names_make_unique()


    print("\n########################### Filtering and Normalizing Matrix ##########################")

    sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    # Filter cells with less than 200 expressed genes and genes which were found in less than 3 cells.
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    # Get barcode from row names
    adata.obs['cell_barcode'] = adata.obs.index


    # Intersect Annotation and Matrix
    bdata = adata[adata.obs['cell_barcode'].isin(anno.cell_barcode)]
    bdata.obs = bdata.obs.merge(anno[['Cas9_Cas12a_targeted','num_features','Gene','Exon','Replicate','Perturbation','cell_type','Project','Sample']], how='left',left_index=True,right_index=True).copy()

    # save data
    bdata.write_h5ad("bdata." + analysis_type + ".h5ad")
    # bdata = sc.read_h5ad("bdata." + analysis_type + ".h5ad")

    #save raw counts in independent column
    bdata.layers['counts'] = bdata.X.copy()

    #normalize to 10,000 reads per cell
    sc.pp.normalize_total(bdata, target_sum=1e4)

    #log transform
    sc.pp.log1p(bdata)

    #feature selection
    sc.pp.highly_variable_genes(bdata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pl.highly_variable_genes(bdata, save=analysis_type + '.genes_variability.pdf')
    bdata.raw = bdata
    bdata = bdata[:, bdata.var.highly_variable]

    # What is regress and scale?
    sc.pp.regress_out(bdata, ["total_counts", "pct_counts_mt"])
    sc.pp.scale(bdata, max_value=10)

    # save scaled data
    bdata.write_h5ad("bdata." + analysis_type + ".scaled.h5ad")
    bdata.write_csvs("bdata." + analysis_type + ".scaled.cells")
    # bdata = sc.read_h5ad(("bdata." + analysis_type + ".scaled.h5ad")


    print("\n########################### Reducing Dimensions ##########################")

    # PCA
    sc.tl.pca(bdata, svd_solver="arpack",random_state=0)
    sc.pl.pca_variance_ratio(bdata, log=True)
    # Compute distances in the PCA space, and find cell neighbors
    sc.pp.neighbors(bdata, n_neighbors=10, n_pcs=40,random_state=0)

    # Generate UMAP features
    sc.tl.umap(bdata,random_state=0)
    sc.pl.umap(bdata, color=["Perturbation","Project","Sample"], use_raw=False, size=2, frameon=True, palette='rainbow', legend_fontsize='xx-small', save=analysis_type + '_origin.pdf')

    # Clustering
    sc.tl.leiden(bdata, resolution=resolution,random_state=0,flavor="igraph",n_iterations=2,directed=False)
    sc.pl.umap(bdata, color=['leiden'], legend_fontsize=8, legend_loc="on data", save=analysis_type + '_leiden.pdf')

    #to_be_removed = {'Intergenic', 'None'}
    #gene_var = list(set(bdata.obs['Exon']))#.remove('intergenic')
    #gene_var = [item for item in gene_var if item not in to_be_removed ]
    #sc.pl.dotplot(
    #    bdata,
    #    groupby="leiden",
    #    var_names=gene_var,
    #    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
    #    save='_base'
    #)

    # Save bdata
    bdata.write_csvs(analysis_type + ".pca.cells")
    bdata.write_h5ad(analysis_type + ".pca.h5ad")
    # bdata = sc.read_h5ad(analysis_type + ".pca.h5ad")


    print("\n########################### Calculating local perturbation signatures ##########################")

    # Run mixscape
    ms = pt.tl.Mixscape()
    ms.perturbation_signature(adata=bdata, pert_key="Exon", control="intergenic",random_state=0)
    bdata_pert = bdata.copy()
    bdata_pert.X = bdata_pert.layers["X_pert"]

    # PCA
    sc.pp.pca(bdata_pert,random_state=0,svd_solver='arpack')
    sc.pp.neighbors(bdata_pert, metric="cosine",random_state=0)

    # UMAP
    sc.tl.umap(bdata_pert,random_state=0)
    sc.pl.umap(bdata_pert, color=["Perturbation","Project","Sample"], palette= 'rainbow', frameon=True, legend_fontsize='xx-small', save=analysis_type + '.pertsig.pdf')

    # Clustering
    sc.tl.leiden(bdata_pert, resolution=resolution,random_state=0,flavor="igraph",n_iterations=2,directed=False)
    sc.pl.umap(bdata_pert, color=['leiden'], legend_fontsize=8, legend_loc="on data", save=analysis_type + '.leiden.pertsig.pdf')
    #sc.pl.dotplot(
    #    bdata_pert,
    #    groupby="leiden",
    #    var_names=gene_var,
    #    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
    #    dendrogram=True,
    #    save='_perturbation_signature'
    #)

    # Save bdata
    bdata_pert.write_csvs(analysis_type + ".pertsig.cells")
    bdata_pert.write_h5ad(analysis_type + ".pertsig.h5ad")
    #bdata_pert = sc.read_h5ad(analysis_type + ".pertsig.h5ad")


    print("\n########################### Running Mixscape ##########################")

    # Run Mixscape
    ms.mixscape(adata=bdata_pert, control="intergenic", labels="Exon")#, layer="X_pert")
    bdata_mixscape = bdata_pert.copy()
    #bdata_mixscape = bdata.copy()
    #bdata_mixscape.X = bdata_mixscape.layers["X_pert"]
    #ms.plot_barplot(adata=bdata, guide_rna_column="Exon") # DOES NOT WORK

    bdata_mixscape.write_csvs(analysis_type + ".mixscape.cells")
    bdata_mixscape.write_h5ad(analysis_type + ".mixscape.h5ad")
    #bdata = sc.read_h5ad(analysis_type + ".mixscape.h5ad")

    # Keep perturbed only cells
    bdata_perturb_only = bdata_mixscape[bdata_mixscape.obs["mixscape_class_global"] == "KO"].copy()

    # PCA
    sc.pp.pca(bdata_perturb_only,random_state=0,svd_solver='arpack')
    sc.pp.neighbors(bdata_perturb_only, metric="cosine",random_state=0)

    # UMAP
    sc.tl.umap(bdata_perturb_only,random_state=0)
    sc.pl.umap(bdata_perturb_only, color=["Perturbation","Project","Sample"], palette= 'rainbow', frameon=True, legend_fontsize='xx-small', save=analysis_type + '.perturebed_only.pdf')

    # Clustering
    sc.tl.leiden(bdata_perturb_only, resolution=resolution,random_state=0,flavor="igraph",n_iterations=2,directed=False)
    sc.pl.umap(bdata_perturb_only, color=['leiden'], legend_fontsize=8, legend_loc="on data", save=analysis_type + '.leiden.perturebed_only.pdf')

    # Save bdata
    bdata_perturb_only.write_csvs(analysis_type + ".perturb_only_cells")
    bdata_perturb_only.write_h5ad(analysis_type + ".perturb_only.h5ad")
    #bdata_perturb_only = sc.read_h5ad(analysis_type + ".perturb_only.h5ad")


    print("\n########################### Computing Correlation Matrices ##########################")


    #  Identify genes that are differentially expressed in the clusters 
    sc.tl.rank_genes_groups(bdata_perturb_only, groupby="leiden", method="wilcoxon")

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

    #sc.pl.rank_genes_groups_heatmap(
    #    bdata_perturb_only,
    #    n_genes=4,
    #    values_to_plot="logfoldchanges",
    #    swap_axes=True,
    #    vmin=-3,
    #    vmax=3,
    #    cmap="bwr",
    #    figsize=(10, 7),
    #    show=False,
    #    save=analysis_type + '.rank_genes_groups_heatmap.pdf'
    #)


    sc.pl.correlation_matrix(bdata_perturb_only, "leiden", figsize=(5, 3.5), save=analysis_type + '.correlation_matrix.pdf')




    print("#Cell Files:\n")
    print(os.path.abspath("./*.perturb_only_cells"))

    print("\n# UMAP and Leiden clustering:")
    print(os.path.abspath("./figures/umap*.pdf"))

    print("\n# Figures after calculating local perturbation signatures:")
    print(os.path.abspath("./figures/umap*.pertsig.pdf"))

    print("\n# Figures using Perturbed cells only:")
    print(os.path.abspath("./figures/umap*.perturebed_only.pdf"))

    print("\n# Dot Plots:")
    print(os.path.abspath("./figures/*.rank_genes_groups_dotplot.pdf"))

    print("\n# Matrix Plots:")
    print(os.path.abspath("./figures/*.rank_genes_groups_matrixplot.pdf"))

    print("\n# Correlation Matrix:")
    print(os.path.abspath("./figures/*.correlation_matrix.pdf"))




if __name__ == "__main__":
        main()
    ######## main - end ################################################################################
