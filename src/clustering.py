import scanpy as sc
import scanpy.external as sce
import matplotlib.pyplot as plt

def cluster(adata, batch_key='sample', res=0.5, n_pcs=50, state_name="PCA_Variance"):
    """
    Run PCA, Harmony batch correction, neighbors, UMAP, and Leiden clustering.

    Parameters
    ----------
    adata : AnnData
    batch_key : str - obs column to correct for (default 'sample')
    res : float - Leiden resolution (default 0.5)
    n_pcs : int - number of PCs to use (default 50)
    state_name : str - name for saving PCA variance plot

    Returns
    -------
    adata : AnnData
    """

    # Run PCA first — Harmony corrects PCA embeddings so this must come before
    sc.tl.pca(adata, svd_solver='arpack', n_comps=n_pcs)

    # Plot PCA variance ratio so user can evaluate how many PCs capture meaningful variance
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=n_pcs)
    plt.savefig(f"results/figures/exploration_figs/{state_name}.png", dpi=300, bbox_inches='tight')
    plt.clf()

    # Harmony batch correction — removes patient-driven variation so clusters reflect biology not donor
    sce.pp.harmony_integrate(adata, batch_key)

    # Build neighbor graph using Harmony-corrected embeddings
    sc.pp.neighbors(adata, use_rep='X_pca_harmony')
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=res)

    return adata