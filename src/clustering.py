import scanpy as sc
import scanpy.external as sce

def cluster(adata, batch_key='sample', resolution=0.5, n_pcs=50):
    """
    Run Harmony batch correction, neighbors, UMAP, and Leiden clustering.

    Parameters
    ----------
    adata : AnnData
    batch_key : str - obs column to correct for (default 'sample')
    resolution : float - Leiden resolution (default 0.5)
    n_pcs : int - number of PCs to use (default 50)

    Returns
    -------
    adata : AnnData
    """
    # TODO: implement
    pass
