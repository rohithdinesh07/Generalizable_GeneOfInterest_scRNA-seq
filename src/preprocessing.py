import scanpy as sc
from scipy.sparse import csr_matrix

def preprocess(adata, HVG=2000):
    """
    Normalize, log-transform, and optionally subset to highly variable genes.
    Run this after filtering, before PCA and clustering.

    Parameters:
        adata : AnnData
        HVG : int or None - how many highly variable genes to keep (default 2000), pass None to skip

    Returns:
        adata : AnnData
    """
    # stash raw counts before touching anything — needed later for DEG tools like DESeq2
    adata.layers["raw_counts"] = adata.X.copy()

    # normalize so every cell has the same total counts, then log transform
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    # keep only the most variable genes — reduces noise and makes PCA/clustering faster
    if HVG is not None:
        sc.pp.highly_variable_genes(adata, n_top_genes=HVG)
        adata = adata[:, adata.var.highly_variable]

    # convert to sparse format to save memory
    adata.X = csr_matrix(adata.X)

    return adata