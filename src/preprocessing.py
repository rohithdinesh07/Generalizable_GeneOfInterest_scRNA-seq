import scanpy as sc
from scipy.sparse import csr_matrix

def preprocess(adata, HVG=2000):
    """
    Normalize, log-transform, select highly variable genes, and compress AnnData.

    Parameters
    ----------
    adata : AnnData
    HVG : int or None - number of highly variable genes to keep (default 2000). Set to None to skip.

    Returns
    -------
    adata : AnnData
    """
    # Save raw counts — required by tools like scVI and DESeq2 that expect raw integer counts
    adata.layers["raw_counts"] = adata.X.copy()

    # Normalize and log transform
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    # Select highly variable genes — reduces noise and speeds up PCA/clustering
    if HVG is not None:
        sc.pp.highly_variable_genes(adata, n_top_genes=HVG)
        adata = adata[:, adata.var.highly_variable]

    # Compress to sparse matrix for local storage efficiency
    adata.X = csr_matrix(adata.X)

    return adata
