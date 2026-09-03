import gseapy as gp
import pandas as pd

def run_enrichment(deg_df, gene_sets=["GO_Biological_Process_2023"], organism="Human", outdir="results/enrichment/"):
    """
    Run GO enrichment on DEG results using gseapy.

    Parameters
    ----------
    deg_df : pd.DataFrame - DEG results with 'names' and 'logfoldchanges' columns
    gene_sets : list - gene set libraries to query (default GO BP 2023)
    organism : str - organism (default 'Human')
    outdir : str - output directory for results

    Returns
    -------
    results : dict - enrichment results per group
    """
    # TODO: implement
    pass
