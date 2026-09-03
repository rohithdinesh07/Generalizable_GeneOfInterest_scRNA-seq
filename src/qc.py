import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt


def plot_qc_metrics(adata):
    """
    Flag mitochondrial and ribosomal genes, compute QC metrics,
    and plot violin plots so the user can choose appropriate thresholds.

    Call this first before filter_cells().

    Parameters
    ----------
    adata : AnnData

    Returns
    -------
    adata : AnnData - with QC metrics added to adata.obs
    """
    # Flag mitochondrial genes (MT-) and ribosomal genes (RPS/RPL)
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))

    # Compute QC metrics — adds pct_counts_mt, pct_counts_ribo etc to adata.obs
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True, log1p=True)

    # Plot distributions so user can decide thresholds
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
    Filter cells based on QC metrics. Call this after plot_qc_metrics()
    once you have chosen appropriate thresholds from the plots.

    Parameters
    ----------
    adata : AnnData - must have QC metrics computed (run plot_qc_metrics first)
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
    print(f"Cells before filtering: {adata.n_obs}")

    # Remove top quantile by gene count — likely doublets or low quality
    upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, upper_quantile)
    adata = adata[adata.obs.n_genes_by_counts < upper_lim]

    # Filter on MT and ribo thresholds
    adata = adata[adata.obs.pct_counts_mt < mt_threshold]
    adata = adata[adata.obs.pct_counts_ribo < ribo_threshold]

    # Filter cells with too few genes and genes expressed in too few cells
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    # scVI doublet detection
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