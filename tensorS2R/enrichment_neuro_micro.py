from pathlib import Path

import cell2cell as c2c
import decoupler as dc
import scanpy as sc

from lr_loadings_utils import (
    compute_prod_dict,
    enr_df,
    filter_and_process_de,
    filter_lr,
    find_threshold_loadings,
    get_expressed_genes,
    get_mean_expression,
    plot_prod_dict,
    process_common_paths,
    process_loadings,
    ranked_genes_df,
)
from plot_utils import plot_pathway_network

GSEA_DBS = ["WikiPathway_2023_Human", "Reactome_2022", "GO_Biological_Process_2023"]

custom_palette = {
    "Astrocyte": "#5D3534",
    "Epithelial": "#FADF81",
    "GABAergic Neuron": "#7BB6F9",
    "Glutamatergic Neuron": "#C2D3E7",
    "Neuron-low count": "#EC68B1",
    "OPC": "#EF8B72",
    "Oligo": "#E8B394",
}


project_path = Path(
    "/sc/arion/projects/psychgen/lbp/data/ScProcesses_brainBlood_nicole/"
)
adata_br = sc.read(
    project_path
    / "batch_2"
    / "c2c_liana_outputs"
    / "data_brain_scran_de_micro_mono.h5ad"
)
adata_br.obs["cell_type_coarse"] = (
    adata_br.obs["cell_type.v2"]
    .astype(str)
    .replace(
        {
            "GABAergic Neuron": "Neurons",
            "Glutamatergic Neuron": "Neurons",
            "Neuron-low count": "Neurons",
        }
    )
)

sc.tl.rank_genes_groups(
    adata_br,
    "cell_type.v2",
    reference="GABAergic Neuron",
    method="wilcoxon",
    layer="scran_normalization",
    pts=True,
    key_added="GABA",
)

sc.tl.rank_genes_groups(
    adata_br,
    "cell_type.v2",
    reference="Glutamatergic Neuron",
    method="wilcoxon",
    layer="scran_normalization",
    pts=True,
    key_added="Gluta",
)

sc.tl.rank_genes_groups(
    adata_br,
    "cell_type.v2",
    reference="Neuron-low count",
    method="wilcoxon",
    layer="scran_normalization",
    pts=True,
    key_added="low",
)

background_micro = get_expressed_genes(
    adata_br, "cell_type_coarse", "Microglia", expr_prop=0.1
)
background_neuro = get_expressed_genes(
    adata_br, "cell_type_coarse", "Neurons", expr_prop=0.1
)

sample_key = "pt"
groupby = "cell_type.v2"
pdata_br = dc.get_pseudobulk(
    adata_br,
    sample_col=sample_key,
    groups_col=groupby,
    layer="counts",
    mode="sum",
    min_cells=10,
    min_counts=10000,
)

df_micro = get_mean_expression(pdata_br, "Microglia", "cell_type.v2")


how = "inner"
cells = "brain_coarse"
batch = "batch_2"
rank = 20

loadings_folder = (
    project_path / batch / "c2c_liana_outputs" / cells / f"rank_{rank}" / how
)
enrichr_folder = str(loadings_folder / "enrichr_results")
factors = c2c.io.load_tensor_factors(loadings_folder / "Loadings.xlsx")
lr_loadings = factors["Ligand-Receptor Pairs"]
thresh = find_threshold_loadings(factors, quantile=0.70)

## Neurons --> Microglia
factor_n2mi = "9"  # "18"

lr_loads_neuro_to_micro = process_loadings(lr_loadings, f"Factor {factor_n2mi}", thresh)

ligands_neuro = lr_loads_neuro_to_micro["ligand"].unique()

target_cells_neuro = ["OPC"]
target_cells_micro = [
    c for c in adata_br.obs["cell_type.v2"].unique() if "Microglia" not in c
]

neuro_ligands = []
for k in ["Gluta", "GABA", "low"]:
    ddf = ranked_genes_df(adata_br, ligands_neuro, key=k)
    ligands = filter_and_process_de(ddf, target_cells_neuro)
    neuro_ligands.append(ligands)

neuro_ligands = list(
    set(neuro_ligands[0]) | set(neuro_ligands[1]) | set(neuro_ligands[2])
)

receptors_micro = (
    lr_loads_neuro_to_micro[lr_loads_neuro_to_micro["ligand"].isin(neuro_ligands)][
        "receptor"
    ]
    .unique()
    .tolist()
)

ddf_micro = ranked_genes_df(adata_br, receptors_micro, key="Microglia")
micro_receptors = filter_and_process_de(
    ddf_micro, ["Astrocyte", "Epithelial", "OPC", "Oligo"]
)

lr_loads_neuro_to_micro_filtered = filter_lr(
    lr_loads_neuro_to_micro.reset_index(), neuro_ligands, micro_receptors
)

lr_loads_neuro_to_micro_filtered.to_csv(
    loadings_folder / f"lr_loads_neuro_to_micro_filtered_{factor_n2mi}.csv"
)

gset_ligands_neuro_to_micro = list(lr_loads_neuro_to_micro_filtered["ligand"].unique())
gset_receptors_neuro_to_micro = list(
    lr_loads_neuro_to_micro_filtered["receptor"].unique()
)

enr_df_ligands_neuro_to_micro = enr_df(
    gset_ligands_neuro_to_micro, GSEA_DBS, background_neuro, enrichr_folder
)
enr_df_receptors_neuro_to_micro = enr_df(
    gset_receptors_neuro_to_micro, GSEA_DBS, background_micro, enrichr_folder
)
common_paths_neuro_to_micro = process_common_paths(
    enr_df_ligands_neuro_to_micro, enr_df_receptors_neuro_to_micro
)

common_paths_neuro_to_micro.to_csv(
    loadings_folder / f"Enriched_paths_neuro_to_micro_{factor_n2mi}.csv"
)

plot_pathway_network(
    common_paths_neuro_to_micro,
    "Neurons -> Microglia",
    loadings_folder / f"Net_paths_neuro_to_micro_{factor_n2mi}.png",
)

df_product = lr_loads_neuro_to_micro_filtered.copy()
df_product["index"] = df_product["ligand"] + "^" + df_product["receptor"]

prod_dict_neuro_to_micro = compute_prod_dict(
    df_product["index"].to_list(),
    pdata_br,
    df_micro,
    "cell_type.v2",
    is_ligand=False,
)

plot_prod_dict(
    prod_dict_neuro_to_micro,
    custom_palette,
    "LR pairs expression products",
    str(loadings_folder / f"filtered_lr_product_dotplot_neuro_to_micro_{factor_n2mi}.png"),
)

## Microglia --> Neurons
target_cells_neuro = [
    c for c in adata_br.obs["cell_type.v2"].unique() if "Neuron" not in c
]

factor_mi2n = "10"

lr_loads_micro_to_neuro = process_loadings(lr_loadings, f"Factor {factor_mi2n}", thresh)

ligands_micro = lr_loads_micro_to_neuro["ligand"].unique()

ddf_micro = ranked_genes_df(adata_br, ligands_micro, key="Microglia")

micro_ligands = ligands_micro  # filter_and_process_de(ddf_micro, target_cells_micro)

receptors_neuro = lr_loads_micro_to_neuro[
    lr_loads_micro_to_neuro["ligand"].isin(micro_ligands)
]["receptor"].unique()

neuro_receptors_ = []
for k in ["Gluta", "GABA", "low"]:
    ddf = ranked_genes_df(adata_br, receptors_neuro, key=k)
    receptors = filter_and_process_de(ddf, target_cells_neuro)
    neuro_receptors_.append(receptors)

neuro_receptors = list(
    set(neuro_receptors_[0]) | set(neuro_receptors_[1]) | set(neuro_receptors_[2])
)

lr_loads_micro_to_neuro_filtered = filter_lr(
    lr_loads_micro_to_neuro.reset_index(), micro_ligands, neuro_receptors
)
lr_loads_micro_to_neuro_filtered.to_csv(
    loadings_folder / f"lr_loads_micro_to_neuro_filtered_{factor_mi2n}.csv"
)

gset_ligands_micro_to_neuro = list(lr_loads_micro_to_neuro_filtered["ligand"].unique())
gset_receptors_micro_to_neuro = list(
    lr_loads_micro_to_neuro_filtered["receptor"].unique()
)

enr_df_ligands_micro_to_neuro = enr_df(
    gset_ligands_micro_to_neuro, GSEA_DBS, background_micro, enrichr_folder
)
enr_df_receptors_micro_to_neuro = enr_df(
    gset_receptors_micro_to_neuro, GSEA_DBS, background_neuro, enrichr_folder
)
common_paths_micro_to_neuro = process_common_paths(
    enr_df_ligands_micro_to_neuro, enr_df_receptors_micro_to_neuro
)

common_paths_micro_to_neuro.to_csv(
    loadings_folder / f"Enriched_paths_micro_to_neuro_{factor_mi2n}.csv"
)

plot_pathway_network(
    common_paths_micro_to_neuro,
    "Microglia -> Neurons",
    loadings_folder / f"Net_paths_micro_to_neuro_{factor_mi2n}.png",
)

df_product = lr_loads_micro_to_neuro_filtered.copy()
df_product["index"] = df_product["ligand"] + "^" + df_product["receptor"]

prod_dict_micro_to_neuro = compute_prod_dict(
    df_product["index"].to_list(),
    pdata_br,
    df_micro,
    "cell_type.v2",
    is_ligand=True,
)

plot_prod_dict(
    prod_dict_micro_to_neuro,
    custom_palette,
    "LR pairs expression products",
    str(loadings_folder / f"filtered_lr_product_dotplot_micro_to_neuro_{factor_mi2n}.png"),
)
