# Generalizable Gene-of-Interest scRNA-seq Pipeline

A modular, reusable single-cell RNA-seq analysis pipeline for investigating any gene of interest across cell populations. Demonstrated using **CELSR2** in ER+ breast cancer (GSE161529 — 117,079 cells across 24 patient-derived tumor samples).

## Biological Context
CELSR2 is a planar cell polarity (PCP) GPCR hypothesized to promote collective cancer cell dissemination in ER+ breast cancer. This pipeline investigates its expression, co-expression patterns, and pathway associations across epithelial tumor cell populations.

## Pipeline Overview
1. **QC** — MT/ribosomal filtering, doublet removal (scVI)
2. **Preprocessing** — normalization, log transformation, compression
3. **Clustering** — Harmony batch correction, UMAP, Leiden clustering
4. **Analysis** — gene scoring, exploratory plots, cluster subsetting
5. **Enrichment** — GO enrichment on cluster DEGs (gseapy)

## Usage
Open `notebooks/01_analysis.ipynb` and set your gene of interest. All pipeline functions accept parameters directly — change thresholds, resolutions, and gene lists interactively.

```python
from src.qc import run_qc
from src.analysis import GeneOfInterestAnalysis

adata = run_qc(adata, mt_threshold=20, min_genes=100)
analysis = GeneOfInterestAnalysis(adata, gene="CELSR2")
analysis.exploratory_plots()
```

## Data
Demo data not included due to file size. Download GSE161529 from GEO and run `notebooks/00_data_loading.ipynb` to generate the processed h5ad file.

## Requirements
```
pip install -r requirements.txt
```

## Project Structure
```
src/          - modular pipeline functions and analysis class
notebooks/    - data loading and interactive analysis
results/      - figures, tables, enrichment outputs (gitignored)
data/demo/    - processed h5ad file (gitignored, generate locally)
```
