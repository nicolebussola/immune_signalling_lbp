import gseapy as gp
import pandas as pd


def get_expressed_genes(adata, cell_type_col, cell_type, expr_prop=0.1):
    """
    Extract genes that are expressed in at least a specified proportion of cells for a given cell type.

    Parameters:
    - adata (AnnData): Annotated data matrix containing gene expression data.
    - cell_type (str): Specific cell type of interest.
    - expr_prop (float): Minimum proportion of cells (0-1) in which a gene must be expressed.
                         Default is 0.1 (10%).

    Returns:
    - list: List of genes that meet the expression proportion threshold in the specified cell type.
    """
    temp = adata[adata.obs[cell_type_col] == cell_type, :]
    a = temp.X.getnnz(axis=0) / temp.X.shape[0]
    stats = (
        pd.DataFrame({"genes": temp.var_names, "props": a})
        .assign(cell_type=cell_type)
        .sort_values("genes")
    )

    stats = stats[stats["props"] >= expr_prop]
    expressed_genes = stats["genes"].values

    return expressed_genes


def enr_df(gset, gene_set, background, out_path):
    """
    Perform enrichment analysis on a gene set and return significant results.

    Parameters:
    gset (list): List of genes to be analyzed.
    gene_set (str): The name of the gene set database to use for enrichment.
    background (list): Background genes
    out_path (str): Output directory path for enrichment results.

    Returns:
    pd.DataFrame: DataFrame containing significant enrichment results sorted by p-value.
    """
    enr = gp.enrichr(
        gene_list=gset,
        gene_sets=gene_set,
        background=background,
        outdir=out_path,
        cutoff=0.05,
        verbose=True,
    )
    return enr.results[enr.results["Adjusted P-value"] < 0.05].sort_values(
        by="P-value", ascending=True
    )


def process_common_paths(enr_df_ligands, enr_df_receptors, no_background=False):
    """
    Process common pathways between ligands and receptors based on enrichment results.

    Parameters:
    enr_df_ligands (pd.DataFrame): Enrichment results DataFrame for ligands.
    enr_df_receptors (pd.DataFrame): Enrichment results DataFrame for receptors.

    Returns:
    pd.DataFrame: DataFrame containing common pathways between ligands and receptors.
    """
    common_pathways_set = list(
        set(enr_df_ligands["Term"]).intersection(set(enr_df_receptors["Term"]))
    )

    common_pathways = pd.merge(
        enr_df_ligands[enr_df_ligands["Term"].isin(common_pathways_set)],
        enr_df_receptors[enr_df_receptors["Term"].isin(common_pathways_set)],
        on=["Term", "Gene_set"],
    )
    common_pathways["Adjusted P-value"] = common_pathways[
        ["Adjusted P-value_x", "Adjusted P-value_y"]
    ].mean(axis=1)
    common_pathways["Genes"] = (
        common_pathways["Genes_x"] + ";" + common_pathways["Genes_y"]
    )
    if no_background:
        common_pathways["Overlap_total"] = (
            common_pathways["Overlap_x"].apply(lambda x: int(x.split("/")[0]))
            + common_pathways["Overlap_y"].apply(lambda x: int(x.split("/")[0]))
        ).astype(str)
        common_pathways["total"] = (
            common_pathways["Overlap_x"]
            .apply(lambda x: int(x.split("/")[1]))
            .astype(str)
        )
        common_pathways["Overlap"] = (
            common_pathways["Overlap_total"] + "/" + common_pathways["total"]
        )

    return common_pathways
