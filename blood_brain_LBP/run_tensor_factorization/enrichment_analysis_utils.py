import gseapy as gp
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


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


def get_mean_expression(pdata, cell_type, celltype_col):
    """
    Get the mean expression of all genes in a specific cell type.

    Parameters:
    pdata (AnnData): The annotated data matrix.
    cell_type (str): The cell type for which to calculate the mean gene expression.
    celltype_col (str): The cell type columnn in adata.obs to group.

    Returns:
    pd.DataFrame: A DataFrame containing the mean expression of all genes in the specified cell type.
    """
    df_c = pd.DataFrame(pdata[pdata.obs[celltype_col] == cell_type].X.mean(axis=0)).T
    df_c.columns = list(pdata.var_names)
    df_c = df_c.rename(index={0: cell_type})
    return df_c


def compute_prod_dict(lr_pairs, pdata, df_micro, celltype_col, is_ligand=True):
    """
    Compute the product expression for ligand-receptor or receptor-ligand pairs.

    Parameters:
    lr_pairs (list): List of ligand-receptor pairs.
    pdata (AnnData): The annotated data matrix.
    df_micro (pd.DataFrame): DataFrame containing microarray data.
    celltype_col (str): The cell type columnn in adata.obs to group.
    is_ligand (bool): Flag to indicate if the computation is for ligand (True) or receptor (False).

    Returns:
    dict: Dictionary with the product of ligand and receptor expressions for each cell type.
    """
    prod_dict = {}
    cell_types = [cl for cl in pdata.obs[celltype_col].unique() if cl != "Microglia"]

    for lr_pair in lr_pairs:
        ligand, receptor = lr_pair.split("^")
        prod_dict[lr_pair] = {}
        for cell_type in cell_types:
            df_c = get_mean_expression(pdata, cell_type, celltype_col)
            if is_ligand:
                prod_dict[lr_pair][cell_type] = df_micro[ligand][0] * df_c[receptor][0]
            else:
                prod_dict[lr_pair][cell_type] = df_micro[receptor][0] * df_c[ligand][0]
    return prod_dict


def plot_prod_dict(prod_dict, palette, title, out_file, figsize=(20, 6)):
    """
    Plot the product expression.

    Parameters:
    prod_dict (dict): Dictionary with the product of ligand and receptor expressions for each cell type.
    palette (dict): Custom palette per cell types
    title (str): The title of the plot.
    out_file (str): The output file path for the plot.

    Returns:
    None
    """
    df = pd.DataFrame(prod_dict)
    df = df.reset_index().melt(
        id_vars="index", var_name="Gene Pair", value_name="Expression"
    )
    df["Expression"] = df["Expression"] / 1e4
    df.rename(columns={"index": "Cell Type"}, inplace=True)

    plt.figure(figsize=figsize)
    sns.scatterplot(
        data=df,
        x="Gene Pair",
        y="Cell Type",
        size="Expression",
        hue="Cell Type",
        palette=palette,
        sizes=(20, 1000),
        legend=False,
    )

    plt.title(title)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    plt.savefig(out_file)


def enr_df(gset, gene_set, background, out_path, organism="human"):
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
        organism=organism,
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
