from pathlib import Path

import cell2cell as c2c
import scanpy as sc

from lr_loadings_utils import (
    enr_df,
    filter_lr,
    find_threshold_loadings,
    get_expressed_genes,
    process_common_paths,
    process_loadings,
    ranked_genes_df,
)
from plot_utils import plot_pathway_network


def process_de(ddf, de_count, logfc_threshold=0.5, pval_threshold=0.05):
    """Filter DE genes: keep those significantly down in non-target types OR with near-zero logFC and non-significant."""
    condition1 = (ddf["logfoldchanges"] < -logfc_threshold) & (ddf["pvals_adj"] < pval_threshold)
    condition2 = (ddf["logfoldchanges"].between(-0.05, 0.05)) & (ddf["pvals_adj"] >= pval_threshold)
    filtered_df = (
        ddf[condition1 | condition2]
        .sort_values(by="names")
        .groupby("names")
        .agg({"count": "sum"})
        == de_count
    )
    return filtered_df[filtered_df["count"]].index.tolist()


def filter_and_process_de(df, cell_types, logfc=0.5, pval=0.05):
    filtered_df = df[df["group"].isin(cell_types)]
    return process_de(filtered_df, de_count=len(cell_types), logfc_threshold=logfc, pval_threshold=pval)


GSEA_DBS = ["WikiPathway_2023_Human", "Reactome_2022", "GO_Biological_Process_2023"]
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
adata_bl = sc.read(
    project_path
    / "batch_2"
    / "c2c_liana_outputs"
    / "data_blood_scran_de_micro_mono.h5ad"
)
adata_bl.obs["cell_type_coarse"] = adata_bl.obs["cell_type.v2"].astype(str).apply(
    lambda x: "Monocytes" if "mono" in x.lower() else x
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

background_mono = get_expressed_genes(
    adata_bl, "cell_type_coarse", "Monocytes", expr_prop=0.1
)
background_neuro = get_expressed_genes(
    adata_br, "cell_type_coarse", "Neurons", expr_prop=0.1
)

how = "inner"
cells = "brain_blood_coarse"
batch = "batch_2"
rank = 20

loadings_folder = (
    project_path / batch / "c2c_liana_outputs" / cells / f"rank_{rank}" / how / "cell_type.v4"
)
enrichr_folder = str(loadings_folder / "enrichr_results")
factors = c2c.io.load_tensor_factors(loadings_folder / "Loadings.xlsx")
lr_loadings = factors["Ligand-Receptor Pairs"]
thresh = find_threshold_loadings(factors, quantile=0.70)

target_cells_neuro = [
    c for c in adata_br.obs["cell_type.v2"].unique() if "Neuron" not in c
]
target_cells_mono = [
    c for c in adata_bl.obs["cell_type.v2"].unique() if "mono" not in c
]

## Neurons --> Monocytes
lr_loads_neuro_to_mono = process_loadings(lr_loadings, "Factor 6", thresh)
ligands_neuro = lr_loads_neuro_to_mono["ligand"].unique()

neuro_ligands = []
for k in ["Gluta", "GABA", "low"]:
    ddf = ranked_genes_df(adata_br, ligands_neuro, key=k)
    ligands = filter_and_process_de(ddf, target_cells_neuro)
    neuro_ligands.append(ligands)

neuro_ligands = list(
    set(neuro_ligands[0]) | set(neuro_ligands[1]) | set(neuro_ligands[2])
)

receptors_mono = (
    lr_loads_neuro_to_mono[lr_loads_neuro_to_mono["ligand"].isin(neuro_ligands)]["receptor"]
    .unique()
    .tolist()
)

ddf_cd14 = ranked_genes_df(adata_bl, receptors_mono, key="CD14")
ddf_cd16 = ranked_genes_df(adata_bl, receptors_mono, key="CD16")

cd14_receptors = filter_and_process_de(ddf_cd14, target_cells_mono)
cd16_receptors = filter_and_process_de(ddf_cd16, target_cells_mono)

mono_receptors = list(set(cd14_receptors) | set(cd16_receptors))

lr_loads_neuro_to_mono_filtered = filter_lr(
    lr_loads_neuro_to_mono.reset_index(), neuro_ligands, mono_receptors
)

lr_loads_neuro_to_mono_filtered.to_csv(
    loadings_folder / "lr_loads_neuro_to_mono_filtered.csv"
)

gset_ligands_neuro_to_mono = list(lr_loads_neuro_to_mono_filtered["ligand"].unique())
gset_receptors_neuro_to_mono = list(
    lr_loads_neuro_to_mono_filtered["receptor"].unique()
)

enr_df_ligands_neuro_to_mono = enr_df(
    gset_ligands_neuro_to_mono, GSEA_DBS, background_neuro, enrichr_folder
)
enr_df_receptors_neuro_to_mono = enr_df(
    gset_receptors_neuro_to_mono, GSEA_DBS, background_mono, enrichr_folder
)
common_paths_neuro_to_mono = process_common_paths(
    enr_df_ligands_neuro_to_mono, enr_df_receptors_neuro_to_mono
)
common_paths_neuro_to_mono.to_csv(
    loadings_folder / "Enriched_paths_neuro_to_mono.csv"
)

plot_pathway_network(
    common_paths_neuro_to_mono,
    "Neurons -> Monocytes",
    loadings_folder / "Net_paths_neuro_to_mono.png",
)

## Monocytes --> Neurons
lr_loads_mono_to_neuro = process_loadings(lr_loadings, "Factor 16", thresh)
ligands_mono = lr_loads_mono_to_neuro["ligand"].unique()

ddf_cd14 = ranked_genes_df(adata_bl, ligands_mono, key="CD14")
ddf_cd16 = ranked_genes_df(adata_bl, ligands_mono, key="CD16")

cd14_ligands = filter_and_process_de(ddf_cd14, target_cells_mono)
cd16_ligands = filter_and_process_de(ddf_cd16, target_cells_mono)

mono_ligands = list(set(cd14_ligands) | set(cd16_ligands))

receptors_neuro = (
    lr_loads_mono_to_neuro[lr_loads_mono_to_neuro["ligand"].isin(mono_ligands)]["receptor"]
    .unique()
    .tolist()
)

neuro_receptors_ = []
for k in ["Gluta", "GABA", "low"]:
    ddf = ranked_genes_df(adata_br, receptors_neuro, key=k)
    receptors = filter_and_process_de(ddf, target_cells_neuro)
    neuro_receptors_.append(receptors)

neuro_receptors = list(
    set(neuro_receptors_[0]) | set(neuro_receptors_[1]) | set(neuro_receptors_[2])
)

lr_loads_mono_to_neuro_filtered = filter_lr(
    lr_loads_mono_to_neuro.reset_index(), mono_ligands, neuro_receptors
)

lr_loads_mono_to_neuro_filtered.to_csv(
    loadings_folder / "lr_loads_mono_to_neuro_filtered.csv"
)

gset_ligands_mono_to_neuro = list(lr_loads_mono_to_neuro_filtered["ligand"].unique())
gset_receptors_mono_to_neuro = list(
    lr_loads_mono_to_neuro_filtered["receptor"].unique()
)

enr_df_ligands_mono_to_neuro = enr_df(
    gset_ligands_mono_to_neuro, GSEA_DBS, background_mono, enrichr_folder
)
enr_df_receptors_mono_to_neuro = enr_df(
    gset_receptors_mono_to_neuro, GSEA_DBS, background_neuro, enrichr_folder
)
common_paths_mono_to_neuro = process_common_paths(
    enr_df_ligands_mono_to_neuro, enr_df_receptors_mono_to_neuro
)
common_paths_mono_to_neuro.to_csv(
    loadings_folder / "Enriched_paths_mono_to_neuro.csv"
)

plot_pathway_network(
    common_paths_mono_to_neuro,
    "Monocytes -> Neurons",
    loadings_folder / "Net_paths_mono_to_neuro.png",
)
