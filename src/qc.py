import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt


def plot_qc_metrics(adata):
    """
    Run this first to get a look at your data before filtering anything.
    Flags MT and ribosomal genes, computes QC metrics, and plots violins
    so you can decide what thresholds make sense for your dataset.

    Parameters:
        adata : AnnData

    Returns:
        adata : AnnData - with QC metrics added to adata.obs
    """
    # tag mitochondrial and ribosomal genes by their name prefixes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))

    # populates adata.obs with pct_counts_mt, n_genes_by_counts etc
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True, log1p=True)

    # look at these plots before deciding your thresholds in filter_cells()
    sc.pl.violin(
        adata,
        ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo'],
        jitter=0.4,
        multi_panel=True
    )
    plt.savefig("results/figures/exploration_figs/MT_Rb_Plots.png", dpi=300, bbox_inches='tight')
    plt.clf()

    return adata


def filter_cells(adata, mt_threshold=20, ribo_threshold=40, min_genes=100, min_cells=10, upper_quantile=0.98, remove_doublets=True):
    """
    Filter out low quality cells based on the QC plots from plot_qc_metrics().
    Adjust thresholds based on what you saw in the violins — defaults are
    reasonable starting points but every dataset is different.

    Parameters:
        adata : AnnData - needs QC metrics already computed, run plot_qc_metrics first
        mt_threshold : float - max % mitochondrial reads allowed (default 20)
        ribo_threshold : float - max % ribosomal reads allowed (default 40)
        min_genes : int - minimum genes a cell needs to express (default 100)
        min_cells : int - minimum cells a gene needs to appear in (default 10)
        upper_quantile : float - drops cells in the top quantile of gene counts (default 0.98)
        remove_doublets : bool - run scVI doublet detection on top of quantile filtering (default True)

    Returns:
        adata : AnnData - filtered
    """
    print(f"Cells before filtering: {adata.n_obs}")

    # cells with way too many genes are probably two cells captured together
    upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, upper_quantile)
    adata = adata[adata.obs.n_genes_by_counts < upper_lim]

    # high MT% usually means the cell is dying or already dead
    adata = adata[adata.obs.pct_counts_mt < mt_threshold]
    adata = adata[adata.obs.pct_counts_ribo < ribo_threshold]

    # drop cells that barely express anything and genes that show up in almost no cells
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    # scVI doublet detection — more rigorous than the quantile cutoff alone
    if remove_doublets:
        try:
            import scvi
            scvi.model.SCVI.setup_anndata(adata)
            vae = scvi.model.SCVI(adata)
            vae.train()
            solo = scvi.external.SOLO.from_scvi_model(vae)
            solo.train()
            df = solo.predict()
            df['prediction'] = solo.predict(soft=False)
            adata.obs['doublet_score'] = df['doublet']
            adata.obs['predicted_doublet'] = df['prediction']
            adata = adata[adata.obs['predicted_doublet'] == 'singlet']
            print(f"Doublets removed via scVI")
        except ImportError:
            print("scvi-tools not installed, skipping doublet removal. Install with: pip install scvi-tools")

    print(f"Cells after filtering: {adata.n_obs}")
    return adata