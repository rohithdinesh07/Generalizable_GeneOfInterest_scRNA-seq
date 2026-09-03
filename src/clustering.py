import scanpy as sc
import scanpy.external as sce
import matplotlib.pyplot as plt

def cluster(adata, batch_key='sample', res=0.5, n_pcs=50, state_name="PCA_Variance"):
    """
    Takes you from preprocessed data all the way to labeled clusters.
    Runs PCA, corrects for batch effects with Harmony, builds the neighbor
    graph, computes UMAP, and runs Leiden clustering.

    Parameters:
        adata : AnnData
        batch_key : str - obs column that identifies batches/samples (default 'sample')
        res : float - Leiden resolution, higher = more clusters (default 0.5)
        n_pcs : int - number of principal components to compute (default 50)
        state_name : str - filename for the PCA variance plot, no extension needed

    Returns:
        adata : AnnData
    """
    # PCA has to come before Harmony since Harmony works by correcting the PCA embeddings
    sc.tl.pca(adata, svd_solver='arpack', n_comps=n_pcs)

    # look at this plot to decide how many PCs are capturing signal vs noise
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=n_pcs)
    plt.savefig(f"results/figures/exploration_figs/{state_name}.png", dpi=300, bbox_inches='tight')
    plt.clf()

    # Harmony pulls apart patient/batch-driven variation so clusters reflect biology, not donor
    sce.pp.harmony_integrate(adata, batch_key)

    # build the neighbor graph from Harmony-corrected PCs, not raw PCA
    sc.pp.neighbors(adata, use_rep='X_pca_harmony')
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=res)

    return adata