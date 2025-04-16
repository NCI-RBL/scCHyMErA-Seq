# scCHyMErA-Seq
> Code repository for scCHyMErA-Seq project<br/>
scCHyMErA-Seq is a platform that efficiently induces exon perturbations as well as gene knockouts to generate single-cell RNA-sequencing phenotypic readouts. To streamline downstream analysis, we have included a ready-to-use pipeline built with scverse tools.<br/>  
The following files are necessary to run the scCHyMErA-Seq pipeline.
1. A matrix file ("*matrix.h5") produced by cellranger.
2. A metadata file with cell barcode and targeting guide information.

## Preparation of matrix file
Please use Cell Ranger's CRISPR Guide Capture Algorithm.
<ins>Run cellranger count function</ins>
```
module load cellranger
cellranger count --id=s \
       --transcriptome=refdata-gex-GRCh38-2024-A \
       --libraries=library.csv \
       --feature-ref=feature_reference.csv \
       --create-bam=true
Example library file (library.csv)
---------------------------------------------
sample,fastqs,lanes,library_type
GEX,Sample_GEX,Any,Gene Expression
Cas9,Sample_Cas9,Any,CRISPR Guide Capture
Cas12a,Sample_Cas12a,Any,CRISPR Guide Capture
---------------------------------------------



#



```
## Loading of matrix file and downstream analysis
### Prerequisites
[scanpy](https://github.com/scverse/scanpy)<br/>
[anndata](https://github.com/scverse/anndata)<br/>
[pertpy](https://github.com/scverse/pertpy) (mixscape analysis)<br/>
[DecoupleR](https://decoupler-py.readthedocs.io/en/latest/installation.html) (pseudobulk count matrix calculation)<br/>
[PyDESeq2](https://pydeseq2.readthedocs.io/en/stable/usage/installation.html) (determinating differentially expressed genes)

### Usage
#### QC plots
```
$ python qc_cells.py filtered_feature_bc_matrix.h5
```

#### Matrix preprocessing and mixscape implementation

$ python scanpy_analysis_split.py<br/>
$ python scanpy_analysis_combined.py

Outputs: UMAPs for all processed cells and LDA plots after applying mixscape.<br/>
**In addition, one LDA plot for each cluster are generated, highlighting the cluster in color while rendering the others in grey to facilitate cluster-specific analysis.**










> #### UMAP and Leiden clustering

<span style="color:blue">Arguments for scanpy_analysis_split.py and scanpy_analysis_combined.py.</span>

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
#### Determination of differentially expressed genes for each perturbation

```
python pseudobulk_deg.py -m filtered_feature_bc_matrix.h5 -a paired_hgRNA_calls_per_cell.csv -p exon_mxs_obs.csv --timestamp $timestamp
```
