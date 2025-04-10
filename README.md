# scCHyMErA-Seq
Code repository for scCHyMErA-Seq project

scCHyMErA-Seq repository includes codes to efficiently filter and process gRNA targeted cells from large amount of single cell sequencing data using the scverse Packages. scCHyMErA-Seq requires cellranger matrix output "*matrix.h5" and a metadata file with cell barcode and targeting guide information.

## Prepare the matrix file

```
#Cell Ranger's CRISPR Guide Capture Algorithm was used

--------------------------------------------------

#example library file (library.csv)

sample,fastqs,lanes,library_type
GEX,Sample_GEX,Any,Gene Expression
Cas9,Sample_Cas9,Any,CRISPR Guide Capture
Cas12a,Sample_Cas12a,Any,CRISPR Guide Capture

-------------------------------------------------

> Run cellranger count function

module load cellranger
cellranger count --id=s \
       --transcriptome=refdata-gex-GRCh38-2024-A \
       --libraries=library.csv \
       --feature-ref=feature_reference.csv \
       --create-bam=true

```

## Prerequisites
[scanpy](https://github.com/scverse/scanpy)

[anndata](https://github.com/scverse/anndata)

[pertpy](https://github.com/scverse/pertpy) (mixscape analysis)

[DecoupleR](https://decoupler-py.readthedocs.io/en/latest/installation.html) (pseudobulk count matrix calculation)

[PyDESeq2](https://pydeseq2.readthedocs.io/en/stable/usage/installation.html) (determinating differentially expressed genes)

## Usage

### QC plots

python qc_cells.py filtered_feature_bc_matrix.h5

### Matrix preprocessing and mixscape implementation

#### LDA

scanpy_mixscpe.py

Outputs: UMAPs for all processed cells and LDA plots after applying mixscape.

**In addition one LDA plot for each cluster are generated, highlighting the cluster in color while rendering the others in grey to facilitate cluster-specific  analysis.**


#### UMAP and Leiden clustering

Arguments for scanpy_analysis_split.py and scanpy_analysis_combined.py.

- -o, --out : Location of output directory where plots will be written. If not specified, files will be written to the current working directory.
- --analysis : KO or Exon analysis. scanpy_analysis_split.py only.
- --resolution : Resolution for leiden clustering. Value between 0 and 1. Higher value will create more clusters.
- -m, --matrix_input : Path to matrix input HDF5 Format
- -a, --anno_csv : Path to annotation matrix input

**Also generate files from enrichment analysis**

Example slurm run:

```

#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=append
#SBATCH --time=1:00:00
#SBATCH --mem=300g
#SBATCH --job-name=schymeraseq


timestamp=$(date +%Y%m%d_%H%M)

export PYTHONHASHSEED=0
export NUMBA_CPU_NAME=generic

timestamp=$(date +%Y%m%d_%H%M)

python scanpy_analysis_split.py -o ./ --analysis Exon --resolution 0.15 -m filtered_feature_bc_matrix.h5 -a paired_hgRNA_calls_per_cell.csv --timestamp $timestamp
python scanpy_analysis_split.py -o ./ --analysis KO --resolution 0.15 -m filtered_feature_bc_matrix.h5 -a paired_hgRNA_calls_per_cell.csv --timestamp $timestamp
python scanpy_analysis_combined.py -o ./ --resolution 0.15 -m filtered_feature_bc_matrix.h5 -a paired_hgRNA_calls_per_cell.csv --timestamp $timestamp

```

### Determination of differentially expressed genes for each perturbation

```
python pseudobulk_deg.py -m filtered_feature_bc_matrix.h5 -a paired_hgRNA_calls_per_cell.csv -p exon_mxs_obs.csv --timestamp $timestamp
```

