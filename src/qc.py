import scanpy as sc
import numpy as np

def run_qc(adata, mt_threshold=20, ribo_threshold=40, min_genes=100, min_cells=10, upper_quantile=0.98, remove_doublets=True):
    """
    Filter cells based on QC metrics.

    Parameters
    ----------
    adata : AnnData
    mt_threshold : float - max % mitochondrial reads (default 20)
    ribo_threshold : float - max % ribosomal reads (default 40)
    min_genes : int - minimum genes expressed per cell (default 100)
    min_cells : int - minimum cells a gene must appear in (default 10)
    upper_quantile : float - upper gene count quantile cutoff to remove likely doublets (default 0.98)
    remove_doublets : bool - whether to run scVI doublet detection (default True)

    Returns
    -------
    adata : AnnData - filtered
    """
    # TODO: implement
    pass
