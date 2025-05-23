# **scCHyMErA-Seq**
> Code repository for the **scCHyMErA-Seq** project

**scCHyMErA-Seq** is a platform that enables efficient exon perturbations and gene knockouts, generating single-cell RNA-sequencing phenotypic readouts. To facilitate downstream analysis, this repository includes a ready-to-use pipeline built with `scverse` tools.

![scCHyMErA_github_04292025](https://github.com/user-attachments/assets/d358caf1-802e-491d-887c-703570b1a90c)


---

## Input Files

To run the scCHyMErA-Seq pipeline, you’ll need:

1. A matrix file (`*matrix.h5`) produced by **Cell Ranger**
2. A metadata file containing cell barcodes and guide information

---

## Generating Input Files

Use the **Cell Ranger `count` pipeline** for CRISPR Guide Capture analysis.
 [Cell Ranger documentation](https://www.10xgenomics.com/support/software/cell-ranger/8.0/analysis/running-pipelines/cr-gex-count)

```bash
module load cellranger
cellranger count --id=s \
    --transcriptome=refdata-gex-GRCh38-2024-A \
    --libraries=library.csv \
    --feature-ref=feature_reference.csv \
    --create-bam=true
```
This will create matrix file and protospacer files along with many others


Example: `library.csv`
```csv
sample,fastqs,lanes,library_type
GEX,Sample_GEX,Any,Gene Expression
Cas9,Sample_Cas9,Any,CRISPR Guide Capture
Cas12a,Sample_Cas12a,Any,CRISPR Guide Capture
```
---

## Loading Files and Downstream Analysis

### Prerequisites

Install the following Python packages:

- [`scanpy`](https://github.com/scverse/scanpy)
- [`anndata`](https://github.com/scverse/anndata)
- [`pertpy`](https://github.com/scverse/pertpy) — *for Mixscape analysis*
- [`DecoupleR`](https://decoupler-py.readthedocs.io/en/latest/installation.html) — *for pseudobulk matrix calculation*
- [`PyDESeq2`](https://pydeseq2.readthedocs.io/en/stable/usage/installation.html) — *for differential expression analysys*


---

## Usage

### Quality Control

```bash
python qc_cells.py filtered_feature_bc_matrix.h5
```

---


### Matrix Preprocessing & Mixscape

```bash
python scanpy_analysis_split.py
python scanpy_analysis_combined.py
```

**Outputs**:
- UMAPs of all processed cells
- Cluster-specific LDA plots (highlighted cluster vs grey others)

---

### UMAP + Leiden Clustering

Arguments for `scanpy_analysis_split.py` and `scanpy_analysis_combined.py`:

| Argument               | Description                                                                 |
|------------------------|-----------------------------------------------------------------------------|
| `-o`, `--out`          | Output directory for plots (default: current working directory)             |
| `--analysis`           | Type of analysis: `KO` or `Exon` (used in `scanpy_analysis_split.py` only) |
| `--resolution`         | Leiden clustering resolution (0–1; higher = more clusters)                  |
| `-m`, `--matrix_input` | Path to input matrix file (`.h5`)                                           |
| `-a`, `--anno_csv`     | Path to annotation file (CSV) with cell barcode and guide pairing           |

> These scripts also generate inputs for `chymeraseq.md` and `gprofiler_analysis.md`


### Example SLURM Job

```bash
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

python scanpy_analysis_split.py -o ./ --analysis Exon --resolution 0.15 \
    -m filtered_feature_bc_matrix.h5 -a paired_hgRNA_calls_per_cell.csv \
    --timestamp $timestamp

python scanpy_analysis_split.py -o ./ --analysis KO --resolution 0.15 \
    -m filtered_feature_bc_matrix.h5 -a paired_hgRNA_calls_per_cell.csv \
    --timestamp $timestamp

python scanpy_analysis_combined.py -o ./ --resolution 0.15 \
    -m filtered_feature_bc_matrix.h5 -a paired_hgRNA_calls_per_cell.csv \
    --timestamp $timestamp
```

---

### Bulk Differential Expression Analysis

To identify differentially expressed genes for each perturbation:

```bash
python pseudobulk_deg.py \
    -m filtered_feature_bc_matrix.h5 \
    -a paired_hgRNA_calls_per_cell.csv \
    -p exon_mxs_obs.csv \
    --timestamp $timestamp
```

---

### Functional Enrichment Analysis

To run functional enrichment analysis using [g:Profiler](https://biit.cs.ut.ee/gprofiler/gost):

```bash
# --excel_file - Excel, csv or txt file generated by scanpy.get.rank_genes_groups_df()

python gprofiler_analysis.py \
--excel_file DEG_exons_mod.csv \
--out deg_exons_mod_0.5 \
--lfc_cutoff 0.5 \
--run_gprofiler

```
