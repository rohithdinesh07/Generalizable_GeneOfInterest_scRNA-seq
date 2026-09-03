# Generalizable Gene-of-Interest scRNA-seq Pipeline

A modular single-cell RNA-seq analysis pipeline built to be dropped into any gene-of-interest study. The pipeline handles everything from raw 10x data through QC, batch correction, clustering, and GO enrichment. Swap your gene and go.

Demonstrated here on **CELSR2** in ER+ primary breast cancer tumors, using the GSE161529 dataset (Qian et al. 2020). 18 ER+ primary tumor samples, ~98,000 cells, 23,276 genes.

---

## Biological Background

CELSR2 is a planar cell polarity (PCP) receptor in the Flamingo/Starry Night GPCR family. It's been implicated in coordinating collective cell migration, and there's growing interest in whether PCP pathway dysregulation contributes to cancer cell dissemination in ER+ breast cancer, a subtype where the metastatic mechanism isn't as well understood as in more aggressive subtypes.

This analysis asks: where is CELSR2 expressed in the ER+ tumor microenvironment, what cell populations drive that expression, and what biological processes are enriched in those clusters?

### What We Found

After subsetting to epithelial/tumor clusters and re-clustering with Harmony batch correction, CELSR2 expression was highest in three populations:

- **Cluster 8 (21.6% expressing)** -- markers: H2AFZ, STMN1, TUBA1B, HMGN2. Proliferating luminal tumor cells (cycling). The highest-expressing cluster, consistent with PCP signaling playing a role in mitotic spindle orientation and cell division.
- **Cluster 12 (18.1%)** -- markers: MGP, PEG10, ZFP36L1, XBP1. Luminal progenitor-like cells. MGP is associated with epithelial differentiation and calcium signaling.
- **Cluster 0 (13.9%)** -- markers: RPL3, KRT19, KRT18, GATA3. Luminal epithelial tumor cells. GATA3 is a canonical luminal A marker.
- **Cluster 2 (11.2%)** -- markers: XBP1, TFF3, SCGB2A2, KRT19. Secretory luminal cells.

Non-epithelial populations (CAFs, macrophages, T cells, mast cells) had near-zero CELSR2 expression, confirming it's tumor-cell specific in this context.

GO enrichment on the high-CELSR2 clusters flagged mitotic spindle organization, chromosome segregation, and DNA replication, all consistent with the proliferating cluster driving expression. There was also enrichment for UPR/ER stress pathways in the luminal progenitor cluster, which lines up with XBP1 expression and luminal progenitor biology.

---

## Pipeline

```
Raw 10x MTX -> QC & filtering -> Normalization -> HVG selection
-> PCA -> Harmony batch correction -> UMAP -> Leiden clustering
-> Gene scoring -> Cluster subsetting -> DEG analysis -> GO enrichment
```

1. **QC** (`src/qc.py`) -- MT/ribosomal filtering, upper quantile cutoff, optional scVI doublet removal
2. **Preprocessing** (`src/preprocessing.py`) -- normalize to total counts, log1p, highly variable gene selection
3. **Clustering** (`src/clustering.py`) -- PCA, Harmony by sample, neighbors, UMAP, Leiden
4. **Analysis** (`src/analysis.py`) -- GeneOfInterestAnalysis class wrapping all plotting and subsetting
5. **Enrichment** (`src/enrichment.py`) -- preranked GSEA via gseapy on per-cluster DEG CSVs

---

## Reproducing This Analysis

### From raw data (full run)

1. Download GSE161529 from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161529
2. Place `GSE161529_RAW/` and `GSE161529_features.tsv` in `data/raw/`
3. Run all cells in `notebooks/01_analysis.ipynb` from the top

### From demo H5AD (skip to analysis)

The demo file (`data/demo/demo_adata.h5ad`) is the full concatenated dataset that's already been QC'd, filtered, and clustered. Sections 1-4 of the notebook are written out but annotated to skip. Load the H5AD in Section 5 and run from there.

The demo file is not committed to this repo (1.9GB). To get it, either run the full pipeline yourself or reach out.

---

## Usage

```python
GENE_OF_INTEREST = 'YOUR_GENE'   # swap this in the config cell

# pipeline
adata = plot_qc_metrics(adata)
adata = filter_cells(adata, mt_threshold=20, remove_doublets=True)
adata = preprocess(adata, HVG=2000)
adata = cluster(adata, batch_key='sample', res=0.5)

# analysis
analysis = GeneOfInterestAnalysis(adata, GENE_OF_INTEREST)
analysis.exploratory_plots()
analysis.score_genes(gene_list=[...], name_list=[...])
sub = analysis.subset_clusters(['0', '2', '8'])
sub.save_deg_results('deg_results_subset')
```

All functions are parameterized. Adjust thresholds, resolutions, and gene lists directly in the notebook.

---

## Requirements

```
pip install -r requirements.txt
```

Core dependencies: scanpy, anndata, harmonypy, gseapy, scvi-tools (optional, for doublet removal)

---

## Project Structure

```
src/              modular pipeline modules and analysis class
notebooks/        full analysis notebook (01_analysis.ipynb)
results/          figures, tables, enrichment outputs (gitignored)
data/demo/        processed h5ad demo file (gitignored)
data/raw/         raw 10x MTX files (gitignored, see data/raw/README.md)
```
