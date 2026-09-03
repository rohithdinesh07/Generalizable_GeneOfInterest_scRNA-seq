import scanpy as sc
import matplotlib.pyplot as plt

class GeneOfInterestAnalysis:
    """
    Interactive analysis and plotting for a gene of interest in scRNA-seq data.

    Parameters
    ----------
    adata : AnnData
    gene : str - gene of interest
    """

    def __init__(self, adata, gene):
        self.adata = adata
        self.gene = gene

    def exploratory_plots(self, color=None, state_name="exploratory_UMAP.png"):
        if color is None:
            color = ["leiden", "sample"]
        """UMAP by leiden, sample, and optional gene of interest."""
        sc.pl.umap(self.adata, color=color, ncols=len(color))
        plt.savefig(f"results/figures/leiden_plots/{state_name}", dpi=300, bbox_inches='tight')
        plt.clf()

    def plot_marker_genes(self, n_genes=5, state_name="top5_genes_leiden.png"):
        """Plot top marker genes per leiden cluster using Wilcoxon rank sum test."""
        sc.tl.rank_genes_groups(self.adata, groupby='leiden', method='wilcoxon')
        sc.pl.rank_genes_groups(self.adata, n_genes=n_genes, sharey=False)
        plt.savefig(f"results/figures/leiden_plots/{state_name}", dpi=300, bbox_inches='tight')
        plt.clf()

    def plot_dotplot(self, var_names, state_name="dotplot.png"):
        """
        Dotplot of selected genes grouped by leiden cluster.

        Parameters
        ----------
        var_names : list - genes to plot
        state_name : str - filename for saved figure
        """
        sc.pl.dotplot(self.adata, var_names=var_names, groupby='leiden', standard_scale='var')
        plt.savefig(f"results/figures/dotplots/{state_name}", dpi=300, bbox_inches='tight')
        plt.clf()

    def score_genes(self, gene_list, name_list, vmax=2, state_name="UMAP_by_score.png"):
        """
        Score cells for gene signatures and plot UMAPs colored by score.

        Parameters
        ----------
        gene_list : list of lists - each inner list is a gene signature
        name_list : list of str - names for each signature score
        vmax : float - max value for color scale (default 2)
        state_name : str - filename for saved figure
        """
        if len(gene_list) != len(name_list):
            print("gene_list and name_list must be the same length.")
            return

        # Score each gene signature and store in adata.obs
        for i in range(len(gene_list)):
            sc.tl.score_genes(self.adata, gene_list=gene_list[i], score_name=name_list[i])

        sc.pl.umap(self.adata, color=name_list, ncols=2, color_map='magma', vmax=vmax, size=20)
        plt.savefig(f"results/figures/score_plots/{state_name}", dpi=300, bbox_inches='tight')
        plt.clf()

    def subset_clusters(self, clusters, res=0.5):
        """
        Subset to specific leiden clusters and rerun PCA + UMAP for focused analysis.

        Parameters
        ----------
        clusters : list of str - leiden cluster ids to keep (e.g. ['0', '2', '4'])
        res : float - Leiden resolution for reclustering (default 0.5)

        Returns
        -------
        GeneOfInterestAnalysis - new analysis object on the subset
        """
        # Subset to chosen clusters
        adata_sub = self.adata[self.adata.obs['leiden'].isin(clusters)].copy()

        # Rerun PCA and UMAP on subset — original embeddings reflect the full dataset
        sc.tl.pca(adata_sub, svd_solver='arpack')
        sc.pp.neighbors(adata_sub)
        sc.tl.umap(adata_sub)
        sc.tl.leiden(adata_sub, resolution=res)

        # Return a new analysis object so all methods are available on the subset
        return GeneOfInterestAnalysis(adata_sub, self.gene)
    
    def save_adata(self, name):
        """
        Save current AnnData object to h5ad file.

        Parameters
        ----------
        name : str
            Filename without extension (e.g. 'adata_epithelial')
        """
        self.adata.write_h5ad(f"data/{name}.h5ad")
        print(f"Saved to data/{name}.h5ad")