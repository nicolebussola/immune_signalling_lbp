import cell2cell as c2c
import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns


def process_loadings(lr_loadings, factor, thresh):
    """
    Process LR loadings for a specific factor and filter based on a threshold.

    Parameters:
    lr_loadings (pd.DataFrame): DataFrame containing ligand-receptor loadings.
    factor (str): The specific factor to be processed.
    thresh (float): Threshold value for filtering loadings.

    Returns:
    pd.DataFrame: Filtered DataFrame with ligand and receptor columns.
    """
    lr_loadings_factor = pd.DataFrame(lr_loadings[factor].sort_values(ascending=False))
    lr_loadings_factor_filtered = lr_loadings_factor.loc[
        lr_loadings_factor[factor] > thresh
    ]
    lr_loadings_factor_filtered["ligand"] = [
        x.split("^")[0] for x in lr_loadings_factor_filtered.index
    ]
    lr_loadings_factor_filtered["receptor"] = [
        x.split("^")[1] for x in lr_loadings_factor_filtered.index
    ]

    lr_loadings_factor_filtered["receptor"] = lr_loadings_factor_filtered[
        "receptor"
    ].apply(lambda x: x.split("_")[0])
    lr_loadings_factor_filtered["ligand"] = lr_loadings_factor_filtered["ligand"].apply(
        lambda x: x.split("_")[0]
    )

    return lr_loadings_factor_filtered


def find_threshold_loadings(factors, plot=False, quantile=0.70):
    """
    Find the threshold loadings based on quantile and optionally plot the distribution.

    Parameters:
    factors (pd.DataFrame): DataFrame containing factors.
    plot (bool): Whether to plot the distribution of network values. Default is False.
    quantile (float): The quantile value to determine the threshold. Default is 0.70.

    Returns:
    float: The threshold value rounded to four decimal places.
    """
    networks = c2c.analysis.tensor_downstream.get_factor_specific_ccc_networks(
        factors,
        sender_label="Sender Cells",
        receiver_label="Receiver Cells",
    )

    network_by_factors = c2c.analysis.tensor_downstream.flatten_factor_ccc_networks(
        networks, orderby="receivers"
    )

    net_values = network_by_factors.values.flatten()
    thresh = np.quantile(net_values, quantile)
    if plot:
        plt.figure()
        _ = plt.hist(net_values, bins=50)
        plt.vlines(thresh, ymin=0, ymax=250, colors="red")
        plt.title(f"t={round(thresh, 4)}")
    return round(thresh, 4)


def extract_ligands_receptors(lr_loadings):
    """
    Extract unique ligands and receptors from loadings.

    Parameters:
    lr_loadings (pd.DataFrame): DataFrame containing ligand-receptor loadings.

    Returns:
    tuple: A tuple containing two lists - cleaned ligands and cleaned receptors.
    """
    ligands = list(lr_loadings["ligand"].unique())
    ligands_clean = [g.split("_")[0] for g in ligands]
    receptors = list(lr_loadings["receptor"].unique())
    receptors_clean = [g.split("_")[0] for g in receptors]
    return ligands_clean, receptors_clean


def ranked_genes_df(adata, geneset, key):
    """
    Create a DataFrame of ranked genes for a given gene set and key.

    Parameters:
    adata (AnnData): Annotated data matrix with DE genes computed.
    geneset (list): List of genes to be ranked.
    key (str): Key to access the ranked genes data.

    Returns:
    pd.DataFrame: DataFrame of ranked genes filtered by the provided gene set.
    """
    ddf = sc.get.rank_genes_groups_df(adata, group=None, key=key)
    ddf = ddf[ddf["names"].isin(geneset)]
    ddf["sig"] = (ddf["pvals_adj"] < 0.05).astype(int)
    ddf["count"] = 1
    return ddf


def process_de_genes_celltype(adata, geneset, n_celltypes, celltype, logfc=0.5, t=0.05):
    """
    Process differentially expressed genes for specific cell type.

    Parameters:
    adata (AnnData): Annotated data matrix.
    geneset (list): List of genes to consider.
    n_celltypes (int): Number of unique cell types.
    celltype (str): Key for specific cell type ('Microglia', 'CD14', or 'CD16').
    logfc (float): Log fold change threshold. Default is 0.5.
    t (float): Threshold for log fold change. Default is 0.05.

    Returns:
    list: List of genes that meet the criteria for differential expression.
    """
    ddf = ranked_genes_df(adata, geneset, key=celltype)

    if celltype == "CD14":
        ddf = ddf[ddf["group"] != "CD16.mono"]
        de_count = n_celltypes - 2
    elif celltype == "CD16":
        ddf = ddf[ddf["group"] != "CD14.mono"]
        de_count = n_celltypes - 2
    else:
        de_count = n_celltypes - 1

    condition1 = (ddf["logfoldchanges"] < -logfc) & (ddf["pvals_adj"] <= t)
    condition2 = (ddf["logfoldchanges"] >= -t) & (ddf["logfoldchanges"] <= t)
    de_df = (
        ddf[condition1 | condition2]
        .sort_values(by="names")
        .groupby("names")
        .agg({"count": "sum"})
        == de_count
    )
    return de_df[de_df["count"]].index.tolist()


def process_de(ddf, de_count, logfc_threshold=0.5, pval_threshold=0.05, t=0.05):
    """Filter differentially expressed genes based on log fold change and p-value."""
    condition1 = (ddf["logfoldchanges"] < -logfc_threshold) & (
        ddf["pvals_adj"] < pval_threshold
    )
    condition2 = (ddf["logfoldchanges"].between(-t, t))

    filtered_df = (
        ddf[condition1 | condition2]
        .sort_values(by="names")
        .groupby("names")
        .agg({"count": "sum"})
        == de_count
    )
    return filtered_df[filtered_df["count"]].index.tolist()


def filter_and_process_de(df, cell_types, logfc=0.5, pval=0.05):
    """Filter and process differentially expressed genes for given cell types."""
    filtered_df = df[df["group"].isin(cell_types)]
    return process_de(
        filtered_df,
        de_count=len(cell_types),
        logfc_threshold=logfc,
        pval_threshold=pval,
    )


def get_unique_receptors(lr_loadings, ligands):
    """
    Get unique receptors from loadings based on provided ligands.

    Parameters:
    lr_loadings (pd.DataFrame): DataFrame containing ligand-receptor loadings.
    ligands (list): List of ligands to filter the receptors.

    Returns:
    np.ndarray: Array of unique receptors.
    """
    return lr_loadings[lr_loadings["ligand"].isin(ligands)]["receptor"].unique()


def get_unique_ligands(lr_loadings, receptors):
    """
    Get unique receptors from loadings based on provided ligands.

    Parameters:
    lr_loadings (pd.DataFrame): DataFrame containing ligand-receptor loadings.
    receptors (list): List of receptors to filter the ligands.

    Returns:
    np.ndarray: Array of unique ligands.
    """
    return lr_loadings[lr_loadings["receptor"].isin(receptors)]["ligand"].unique()


def filter_lr(loadings, ligands, receptors):
    """
    Filter ligand-receptor loadings based on provided ligands and receptors.

    Parameters:
    loadings (pd.DataFrame): DataFrame containing ligand-receptor loadings.
    ligands (list): List of ligands to filter the loadings.
    receptors (list): List of receptors to filter the loadings.

    Returns:
    pd.DataFrame: Filtered DataFrame of ligand-receptor pairs.
    """
    df_lr = loadings[
        (loadings["ligand"].isin(ligands)) & (loadings["receptor"].isin(receptors))
    ]
    factor = loadings.columns[1]
    df_lr = df_lr.reset_index(drop=True).rename(
        {"Unnamed: 0": "L-R pairs", factor: f"'{factor}' Loadings"}, axis=1
    )
    return df_lr


def lr_from_micro(
    adata_brain,
    lr_loadings,
    genes_ligands,
    n_celltypes_brain,
):
    """
    Generate ligand-receptor interactions from microglia.

    Parameters:
    adata_brain (AnnData): Annotated adata for brain.
    lr_loadings (pd.DataFrame): DataFrame containing ligand-receptor loadings.
    genes_ligands (list): List of ligands.
    n_celltypes_brain (int): Number of cell types for brain.

    Returns:
    pd.DataFrame: Filtered DataFrame of ligand-receptor interactions.
    """

    micro_ligands = process_de_genes_celltype(
        adata_brain, genes_ligands, n_celltypes_brain, celltype="Microglia"
    )
    all_receptors = get_unique_receptors(lr_loadings, micro_ligands)
    return filter_lr(lr_loadings, micro_ligands, all_receptors)


def lr_to_micro(
    adata_brain,
    lr_loadings,
    genes_receptors,
    n_celltypes_brain,
):
    """
    Generate ligand-receptor interactions from microglia.

    Parameters:
    adata_brain (AnnData): Annotated adata for brain.
    lr_loadings (pd.DataFrame): DataFrame containing ligand-receptor loadings.
    genes_ligands (list): List of ligands.
    n_celltypes_brain (int): Number of cell types for brain.

    Returns:
    pd.DataFrame: Filtered DataFrame of ligand-receptor interactions.
    """

    micro_receptors = process_de_genes_celltype(
        adata_brain, genes_receptors, n_celltypes_brain, celltype="Microglia"
    )
    gene_ligands = get_unique_ligands(lr_loadings, micro_receptors)
    return filter_lr(lr_loadings, gene_ligands, micro_receptors)


def lr_from_micro_to_mono(
    adata_ligands,
    adata_receptors,
    lr_loadings,
    genes_ligands,
    n_celltypes_ligands,
    n_celltypes_receptors,
):
    """
    Generate ligand-receptor interactions from microglia to monocytes.

    Parameters:
    adata_ligands (AnnData): Annotated data matrix for ligands.
    adata_receptors (AnnData): Annotated data matrix for receptors.
    lr_loadings (pd.DataFrame): DataFrame containing ligand-receptor loadings.
    genes_ligands (list): List of ligands.
    n_celltypes_ligands (int): Number of cell types for ligands.
    n_celltypes_receptors (int): Number of cell types for receptors.

    Returns:
    pd.DataFrame: Filtered DataFrame of ligand-receptor interactions.
    """

    micro_ligands = process_de_genes_celltype(
        adata_ligands, genes_ligands, n_celltypes_ligands, celltype="Microglia"
    )
    blood_receptors = get_unique_receptors(lr_loadings, micro_ligands)
    mono_cd16_receptors = process_de_genes_celltype(
        adata_receptors, blood_receptors, n_celltypes_receptors, celltype="CD16"
    )
    mono_cd14_receptors = process_de_genes_celltype(
        adata_receptors, blood_receptors, n_celltypes_receptors, celltype="CD14"
    )
    mono_receptors = list(set(mono_cd16_receptors).union(set(mono_cd14_receptors)))
    return filter_lr(lr_loadings, micro_ligands, mono_receptors)


def lr_from_mono_to_micro(
    adata_ligands,
    adata_receptors,
    lr_loadings,
    genes_ligands,
    n_celltypes_ligands,
    n_celltypes_receptors,
):
    """
    Generate ligand-receptor interactions to microglia from monocytes.

    Parameters:
    adata_ligands (AnnData): Annotated data matrix for ligands.
    adata_receptors (AnnData): Annotated data matrix for receptors.
    lr_loadings (pd.DataFrame): DataFrame containing ligand-receptor loadings.
    genes_ligands (list): List of ligands.
    n_celltypes_ligands (int): Number of cell types for ligands.
    n_celltypes_receptors (int): Number of cell types for receptors.

    Returns:
    pd.DataFrame: Filtered DataFrame of ligand-receptor interactions.
    """

    mono_cd16_ligands = process_de_genes_celltype(
        adata_ligands, genes_ligands, n_celltypes_ligands, celltype="CD16"
    )
    mono_cd14_ligands = process_de_genes_celltype(
        adata_ligands, genes_ligands, n_celltypes_ligands, celltype="CD14"
    )
    mono_ligands = list(set(mono_cd16_ligands).union(set(mono_cd14_ligands)))

    brain_receptors = get_unique_receptors(lr_loadings, mono_ligands)
    micro_receptors = process_de_genes_celltype(
        adata_receptors, brain_receptors, n_celltypes_receptors, celltype="Microglia"
    )
    return filter_lr(lr_loadings, mono_ligands, micro_receptors)


def select_top_n(d, n=None):
    """
    Select top N items from a dictionary based on their absolute values.

    Parameters:
    d (dict): Dictionary of items to be sorted.
    n (int): Number of top items to select. Default is None, which selects all items.

    Returns:
    dict: Dictionary containing the top N items.
    """
    d = dict(sorted(d.items(), key=lambda item: abs(item[1]), reverse=True))
    return {k: v for i, (k, v) in enumerate(d.items()) if i < n}


def filter_de_genes_micro(adata, geneset, lr_df, ligand_or_receptor, logfc=0.5):
    """
    Filter differentially expressed genes for microglia from the ligand-receptor DataFrame.

    Parameters:
    adata (AnnData): Annotated data matrix.
    geneset (list): List of genes to be filtered.
    lr_df (pd.DataFrame): DataFrame containing ligand-receptor pairs.
    ligand_or_receptor (str): Column name to filter ('ligand' or 'receptor').
    logfc (float): Log fold change threshold. Default is 0.5.

    Returns:
    pd.DataFrame: Filtered DataFrame of ligand-receptor pairs.
    """
    ddf = ranked_genes_df(adata, geneset, key="Monocytes")
    ddf_micro = ddf[ddf["group"] == "Microglia"]
    genes_to_exclude = (
        ddf_micro[ddf_micro["logfoldchanges"] > logfc]
        .sort_values(by="names")["names"]
        .unique()
    )
    return lr_df[~lr_df[ligand_or_receptor].isin(genes_to_exclude)].reset_index(
        drop=True
    )


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
    return stats["genes"].values


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
