import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd

class GeneOfInterestAnalysis:
    """
    Analysis and plotting class built around a gene of interest.
    Initialize with your adata and the gene you care about, then
    call methods interactively from the notebook.

    Parameters:
        adata : AnnData
        gene : str - gene of interest
    """

    def __init__(self, adata, gene):
        self.adata = adata
        self.gene = gene

    def exploratory_plots(self, color=None, state_name="exploratory_UMAP"):
        """
        First look at your clusters — plots UMAP colored by leiden and sample by default.

        Parameters:
            color : list - obs columns to color by (default ['leiden', 'sample'])
            state_name : str - filename without extension
        """
        if color is None:
            color = ["leiden", "sample"]
        sc.pl.umap(self.adata, color=color, ncols=len(color))
        plt.savefig(f"results/figures/leiden_plots/{state_name}.png", dpi=300, bbox_inches='tight')
        plt.clf()

    def plot_marker_genes(self, n_genes=5, state_name="top5_genes_leiden"):
        """
        Runs Wilcoxon and plots the top marker genes per cluster.
        Good for getting a quick read on what each cluster might be.

        Parameters:
            n_genes : int - how many top genes to show per cluster (default 5)
            state_name : str - filename without extension
        """
        sc.tl.rank_genes_groups(self.adata, groupby='leiden', method='wilcoxon')
        sc.pl.rank_genes_groups(self.adata, n_genes=n_genes, sharey=False)
        plt.savefig(f"results/figures/leiden_plots/{state_name}.png", dpi=300, bbox_inches='tight')
        plt.clf()

    def plot_dotplot(self, var_names, state_name="dotplot"):
        """
        Dotplot of specific genes across clusters — useful for validating
        known marker genes or checking your gene of interest expression.

        Parameters:
            var_names : list - genes to plot
            state_name : str - filename without extension
        """
        sc.pl.dotplot(self.adata, var_names=var_names, groupby='leiden', standard_scale='var', layer='counts')
        plt.savefig(f"results/figures/dotplots/{state_name}.png", dpi=300, bbox_inches='tight')
        plt.clf()

    def score_genes(self, gene_list, name_list, vmax=2, state_name="UMAP_by_score"):
        """
        Score cells for gene signatures and plot UMAPs colored by those scores.
        Pass multiple signatures at once and they'll all get plotted together.

        Parameters:
            gene_list : list of lists - each inner list is one gene signature
            name_list : list of str - name for each signature score
            vmax : float - max value for the color scale (default 2)
            state_name : str - filename without extension
        """
        if len(gene_list) != len(name_list):
            print("gene_list and name_list must be the same length.")
            return

        for i in range(len(gene_list)):
            sc.tl.score_genes(self.adata, gene_list=gene_list[i], score_name=name_list[i])

        sc.pl.umap(self.adata, color=name_list, ncols=2, color_map='magma', vmax=vmax, size=20)
        plt.savefig(f"results/figures/score_plots/{state_name}.png", dpi=300, bbox_inches='tight')
        plt.clf()

    def subset_clusters(self, clusters, res=0.5):
        """
        Pull out specific clusters and rerun PCA + UMAP just on that subset.
        Returns a new GeneOfInterestAnalysis object so you can keep using all the same methods.

        Parameters:
            clusters : list of str - leiden cluster ids to keep (e.g. ['0', '2', '4'])
            res : float - Leiden resolution for reclustering the subset (default 0.5)

        Returns:
            GeneOfInterestAnalysis - new object on the subset
        """
        adata_sub = self.adata[self.adata.obs['leiden'].isin(clusters)].copy()

        # rerun PCA and UMAP — the original embeddings were fit on the full dataset
        sc.tl.pca(adata_sub, svd_solver='arpack')
        sc.pp.neighbors(adata_sub)
        sc.tl.umap(adata_sub)
        sc.tl.leiden(adata_sub, resolution=res)

        return GeneOfInterestAnalysis(adata_sub, self.gene)

    def save_adata(self, name):
        """
        Save the current adata to an h5ad file — useful for checkpointing
        after a long step so you don't have to rerun everything.

        Parameters:
            name : str - filename without extension (e.g. 'adata_epithelial')
        """
        self.adata.write_h5ad(f"data/{name}.h5ad")
        print(f"Saved to data/{name}.h5ad")

    def save_deg_results(self, name):
        """
        Runs Wilcoxon DEG across all clusters and saves everything to one CSV.
        Each row is a gene-cluster pair so you can filter by cluster downstream.

        Parameters:
            name : str - filename without extension (e.g. 'deg_results')
        """
        sc.tl.rank_genes_groups(self.adata, groupby='leiden', method='wilcoxon')

        all_clusters = self.adata.obs['leiden'].unique().tolist()

        dfs = []
        for cluster in all_clusters:
            df = sc.get.rank_genes_groups_df(self.adata, group=cluster)
            df['cluster'] = cluster
            dfs.append(df)

        results = pd.concat(dfs, ignore_index=True)
        results.to_csv(f"results/tables/{name}.csv", index=False)
        print(f"Saved DEG results to results/tables/{name}.csv")