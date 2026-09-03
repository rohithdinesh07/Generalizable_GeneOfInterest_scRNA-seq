import scanpy as sc
import matplotlib.pyplot as plt

class GeneOfInterestAnalysis:
    """
    Interactive analysis and plotting for a gene of interest in scRNA-seq data.

    Parameters
    ----------
    adata : AnnData
    gene : str - gene of interest (e.g. 'CELSR2')
    """

    def __init__(self, adata, gene):
        self.adata = adata
        self.gene = gene

    # TODO: add methods for:
    # - exploratory_plots() - UMAP by leiden, sample, gene of interest
    # - plot_marker_genes() - dotplot by leiden
    # - score_genes() - proliferative, dormant, custom gene lists
    # - subset_clusters() - create new object with specific leiden inputs + PCA
