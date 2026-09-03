import gseapy as gp
import pandas as pd
import matplotlib.pyplot as plt

def run_enrichment(deg_df, cluster_id, gene_sets=["GO_Biological_Process_2023"], organism="Human", outdir="results/enrichment/", state_name="GSEA_results"):
    """
    Runs preranked GSEA on a specific cluster using the DEG CSV from save_deg_results().
    Genes are ranked by log2FC — no arbitrary cutoff needed, every gene contributes.

    Parameters:
        deg_df : pd.DataFrame - DEG results loaded from CSV
        cluster_id : str - which leiden cluster to run enrichment on (e.g. '0')
        gene_sets : list - which gene set libraries to query (default GO BP 2023)
        organism : str - organism (default 'Human')
        outdir : str - where to save GSEA output files
        state_name : str - filename without extension for saved outputs

    Returns:
        prerank_results : gseapy prerank object - pass this to specific_enrichment()
        results_df : pd.DataFrame - enrichment results with NES, pval, fdr per pathway
    """
    # filter down to just the cluster we care about
    df = deg_df[deg_df['cluster'] == cluster_id].copy()

    # rank genes by log2FC — highest at the top
    ranked_genes = df.set_index('names')['logfoldchanges'].sort_values(ascending=False)

    # run GSEA
    prerank_results = gp.prerank(rnk=ranked_genes, gene_sets=gene_sets, organism=organism, outdir=outdir, seed=42)

    # pull out the results table and save it
    results_df = prerank_results.res2d
    results_df.to_csv(f"results/tables/{state_name}.csv", index=False)

    # quick bar plot of top 10 pathways by NES score
    top = results_df.sort_values('NES', ascending=False).head(10)
    top.plot(kind='barh', x='Term', y='NES')
    plt.savefig(f"results/figures/enrichment_plots/{state_name}_top_pathways.png", dpi=300, bbox_inches='tight')
    plt.clf()

    return prerank_results, results_df


def specific_enrichment(prerank_results, term, outdir="results/figures/enrichment_plots/"):
    """
    Plot the enrichment curve for one specific pathway from your GSEA results.
    Grab the term name from prerank_results.res2d['Term'].

    Parameters:
        prerank_results : gseapy prerank object - output from run_enrichment()
        term : str - exact pathway name you want to plot
        outdir : str - where to save the figure
    """
    from gseapy import gseaplot

    gseaplot(
        rank_metric=prerank_results.ranking,
        term=term,
        ofname=f"{outdir}{term.replace(' ', '_')}_curve.png",
        **prerank_results.results[term]
    )
    print(f"Saved enrichment curve for: {term}")