import logging
import os
import warnings
from pathlib import Path

import cell2cell as c2c
from .lr_loadings_utils import (
    enr_df,
    filter_and_process_de,
    filter_lr,
    find_threshold_loadings,
    get_expressed_genes,
    get_unique_receptors,
    process_common_paths,
    process_de_genes_celltype,
    process_loadings,
    ranked_genes_df,
)


GSEA_DBS = ["WikiPathway_2023_Human", "Reactome_2022", "GO_Biological_Process_2023"]

logging.disable(logging.WARNING)


warnings.filterwarnings("ignore")


CELLS_FACTOR_DICT = {
    "micro_blood_coarse": {
        "cohort_1": {"Micro_to_Mono": "Factor 8", "Mono_to_Micro": "Factor 3"},
        "cohort_2": {"Micro_to_Mono": "Factor 5", "Mono_to_Micro": "Factor 3"},
    },
    "brain_blood_coarse": {
        "cohort_1": {"Micro_to_Mono": "Factor 12", "Mono_to_Micro": "Factor 15"},
        "cohort_2": {"Micro_to_Mono": "Factor 6", "Mono_to_Micro": "Factor 17"},
    },
    "brain_coarse": {
        "cohort_1": {"Micro_to_Mono": "Factor 4", "Mono_to_Micro": "Factor 3"}
    },
}


def create_enrichment_directory(loadings_path) -> Path:
    """Create directory for enrichment results."""
    enrich_path = Path(loadings_path) / "enrichr_results"
    enrich_path.mkdir(parents=True, exist_ok=True)
    return enrich_path


def enrich_mono_micro_brain(
    cohort, loadings_path, adata_br, genedbs, design, out_folder
):

    background_micro = get_expressed_genes(
        adata_br, "cell_type_coarse", "Microglia", expr_prop=0.1
    )
    background_mono = get_expressed_genes(
        adata_br, "cell_type_coarse", "Monocytes", expr_prop=0.1
    )

    factors = c2c.io.load_tensor_factors(loadings_path / "Loadings.xlsx")
    lr_loadings = factors["Ligand-Receptor Pairs"]
    thresh = find_threshold_loadings(factors, quantile=0.70)

    ## Monocytes --> Microglia
    mo2mi_factor = CELLS_FACTOR_DICT[design][cohort]["Mono_to_Micro"]
    lr_loads_mono_to_micro = process_loadings(lr_loadings, mo2mi_factor, thresh)
    receiver_cells_mo2mi = factors["Receiver Cells"][mo2mi_factor] > thresh

    receptors_micro = lr_loads_mono_to_micro["receptor"].unique()
    ddf_micro = ranked_genes_df(adata_br, receptors_micro, key="Microglia")
    micro_receptors = filter_and_process_de(
        ddf_micro, receiver_cells_mo2mi[receiver_cells_mo2mi].index.tolist()
    )

    ligands_mono = lr_loads_mono_to_micro[
        lr_loads_mono_to_micro["receptor"].isin(micro_receptors)
    ]["ligand"].unique()

    target_cells_mono = [
        c for c in adata_br.obs["cell_type.v2"].unique() if "mono" not in c
    ]
    mono_ligands_ = []
    for k in ["CD16", "CD14"]:
        ddf = ranked_genes_df(adata_br, ligands_mono, key=k)
        ligands = filter_and_process_de(ddf, target_cells_mono)
        mono_ligands_.append(ligands)

    mono_ligands = list(set(mono_ligands_[0]) | set(mono_ligands_[1]))

    lr_loads_mono_to_micro_filtered = filter_lr(
        lr_loads_mono_to_micro.reset_index(), mono_ligands, micro_receptors
    )
    lr_loads_mono_to_micro_filtered.to_csv(
        loadings_path / "lr_loads_mono_to_micro_filtered.csv"
    )

    gset_ligands_mono_to_micro = list(
        lr_loads_mono_to_micro_filtered["ligand"].unique()
    )
    gset_receptors_mono_to_micro = list(
        lr_loads_mono_to_micro_filtered["receptor"].unique()
    )

    enr_df_ligands_mono_to_micro = enr_df(
        gset_ligands_mono_to_micro,
        genedbs,
        background_mono,
        out_folder,
    )
    enr_df_receptors_mono_to_micro = enr_df(
        gset_receptors_mono_to_micro,
        genedbs,
        background_micro,
        out_folder,
    )
    common_paths_mono_to_micro = process_common_paths(
        enr_df_ligands_mono_to_micro, enr_df_receptors_mono_to_micro
    )
    common_paths_mono_to_micro.to_csv(out_folder / "Enriched_paths_mono_to_micro.csv")

    plot_pathway_network(
        common_paths_mono_to_micro,
        "Monocytes -> Microglia",
        out_folder / "Net_paths_mono_to_micro.png",
    )

    ## Microglia --> Monocytes
    mi2mo_factor = CELLS_FACTOR_DICT[design][cohort]["Micro_to_Mono"]
    lr_loads_micro_to_mono = process_loadings(lr_loadings, mi2mo_factor, thresh)
    sender_cells_mi2mo = factors["Sender Cells"][mi2mo_factor] > thresh

    ligands_micro = lr_loads_micro_to_mono["ligand"].unique()
    ddf_micro = ranked_genes_df(adata_br, ligands_micro, key="Microglia")
    micro_ligands = filter_and_process_de(
        ddf_micro, sender_cells_mi2mo[sender_cells_mi2mo].index.tolist()
    )

    receptors_mono = lr_loads_micro_to_mono[
        lr_loads_micro_to_mono["ligand"].isin(micro_ligands)
    ]["receptor"].unique()

    mono_receptors_ = []
    for k in ["CD16", "CD14"]:
        ddf = ranked_genes_df(adata_br, receptors_mono, key=k)
        receptors = filter_and_process_de(ddf, target_cells_mono)
        mono_receptors_.append(receptors)

    mono_receptors = list(set(mono_receptors_[0]) | set(mono_receptors_[1]))

    lr_loads_micro_to_mono_filtered = filter_lr(
        lr_loads_micro_to_mono.reset_index(), micro_ligands, mono_receptors
    )
    lr_loads_micro_to_mono_filtered.to_csv(
        loadings_path / "lr_loads_micro_to_mono_filtered.csv"
    )

    gset_ligands_micro_to_mono = list(
        lr_loads_micro_to_mono_filtered["ligand"].unique()
    )
    gset_receptors_micro_to_mono = list(
        lr_loads_micro_to_mono_filtered["receptor"].unique()
    )

    enr_df_ligands_micro_to_mono = enr_df(
        gset_ligands_micro_to_mono, genedbs, background_micro, out_folder
    )
    enr_df_receptors_micro_to_mono = enr_df(
        gset_receptors_micro_to_mono, genedbs, background_mono, out_folder
    )
    common_paths_micro_to_mono = process_common_paths(
        enr_df_ligands_micro_to_mono, enr_df_receptors_micro_to_mono
    )
    common_paths_micro_to_mono.to_csv(out_folder / "Enriched_paths_micro_to_mono.csv")




def enrich_mono_micro_blood_micro(
    cohort, loadings_path, adata_br, adata_bl, adata_tf, design, out_folder, genedbs=GSEA_DBS
):

    n_celltypes_brain = adata_br.obs["cell_type.v2"].nunique()
    background_micro = get_expressed_genes(
        adata_br, "cell_type_coarse", "Microglia", expr_prop=0.1
    )
    background_mono = get_expressed_genes(
        adata_bl, "cell_type_coarse", "Monocytes", expr_prop=0.1
    )

    factors = c2c.io.load_tensor_factors(loadings_path / "Loadings.xlsx")
    thresh = find_threshold_loadings(factors, quantile=0.70)
    lr_loadings = factors["Ligand-Receptor Pairs"]

    # Monocytes to Microglia
    mo2mi_factor = CELLS_FACTOR_DICT[design][cohort]["Mono_to_Micro"]
    lr_loads_mono_to_micro = process_loadings(
        lr_loadings, mo2mi_factor, thresh
    ).reset_index()

    target_cells_mono_bl = [
        c for c in adata_bl.obs["cell_type.v2"].unique() if "mono" not in c
    ]

    all_mono_ligands = lr_loads_mono_to_micro["ligand"].unique().tolist()

    mono_ligands_bl_ = []
    for k in ["CD16", "CD14"]:
        ddf = ranked_genes_df(adata_bl, all_mono_ligands, key=k)
        ligands = filter_and_process_de(ddf, target_cells_mono_bl)
        mono_ligands_bl_.append(ligands)

    mono_ligands_bl = set(mono_ligands_bl_[0]) | set(mono_ligands_bl_[1])

    target_cells_mono_tf = [
        c for c in adata_tf.obs["cell_type.v2"].unique() if "mono" not in c
    ]
    mono_ligands_tf_ = []
    for k in ["CD16", "CD14"]:
        ddf = ranked_genes_df(adata_tf, all_mono_ligands, key=k)
        ligands = filter_and_process_de(ddf, target_cells_mono_tf)
        mono_ligands_tf_.append(ligands)

    mono_ligands_tf = set(mono_ligands_tf_[0]) | set(mono_ligands_tf_[1])

    mono_ligands = list(mono_ligands_bl & mono_ligands_tf)

    brain_receptors = get_unique_receptors(lr_loads_mono_to_micro, mono_ligands)
    micro_receptors = process_de_genes_celltype(
        adata_br, brain_receptors, n_celltypes_brain, celltype="Microglia"
    )
    lr_loads_mono_to_micro_filtered = filter_lr(
        lr_loads_mono_to_micro, mono_ligands, micro_receptors
    )
    lr_loads_mono_to_micro_filtered.to_csv(
        loadings_path / "lr_loads_mono_to_micro_filtered.csv"
    )

    gset_ligands_mono_to_micro = list(
        lr_loads_mono_to_micro_filtered["ligand"].unique()
    )
    gset_receptors_mono_to_micro = list(
        lr_loads_mono_to_micro_filtered["receptor"].unique()
    )

    enr_df_ligands_mono_to_micro = enr_df(
        gset_ligands_mono_to_micro,
        genedbs,
        background_mono,
        out_folder,
    )
    enr_df_receptors_mono_to_micro = enr_df(
        gset_receptors_mono_to_micro,
        genedbs,
        background_micro,
        out_folder,
    )
    common_paths_mono_to_micro = process_common_paths(
        enr_df_ligands_mono_to_micro, enr_df_receptors_mono_to_micro
    )
    common_paths_mono_to_micro.to_csv(out_folder / "Enriched_paths_mono_to_micro.csv")



    # Microglia to Monocytes

    mi2mo_factor = CELLS_FACTOR_DICT[design][cohort]["Micro_to_Mono"]
    lr_loads_micro_to_mono = process_loadings(
        lr_loadings, mi2mo_factor, thresh
    ).reset_index()

    brain_ligands = lr_loads_micro_to_mono["ligand"].unique().tolist()
    micro_ligands = process_de_genes_celltype(
        adata_br, brain_ligands, n_celltypes_brain, celltype="Microglia"
    )
    all_receptors = get_unique_receptors(lr_loads_micro_to_mono, micro_ligands)

    mono_receptors_bl_ = []
    for k in ["CD16", "CD14"]:
        ddf = ranked_genes_df(adata_bl, all_receptors, key=k)
        receptors = filter_and_process_de(ddf, target_cells_mono_bl)
        mono_receptors_bl_.append(receptors)

    mono_receptors_bl = set(mono_receptors_bl_[0]) | set(mono_receptors_bl_[1])

    target_cells_mono_tf = [
        c for c in adata_tf.obs["cell_type.v2"].unique() if "mono" not in c
    ]
    mono_receptors_tf_ = []
    for k in ["CD16", "CD14"]:
        ddf = ranked_genes_df(adata_tf, all_receptors, key=k)
        receptors = filter_and_process_de(ddf, target_cells_mono_tf)
        mono_receptors_tf_.append(receptors)

    mono_receptors_tf = set(mono_receptors_tf_[0]) | set(mono_receptors_tf_[1])

    mono_receptors = list(mono_receptors_bl & mono_receptors_tf)

    lr_loads_micro_to_mono_filtered = filter_lr(
        lr_loads_micro_to_mono, micro_ligands, mono_receptors
    )
    lr_loads_micro_to_mono_filtered.to_csv(
        loadings_path / "lr_loads_micro_to_mono_filtered.csv"
    )
    gset_ligands_micro_to_mono = list(
        lr_loads_micro_to_mono_filtered["ligand"].unique()
    )
    gset_receptors_micro_to_mono = list(
        lr_loads_micro_to_mono_filtered["receptor"].unique()
    )

    enr_df_ligands_micro_to_mono = enr_df(
        gset_ligands_micro_to_mono,
        genedbs,
        background_micro,
        out_folder,
    )
    enr_df_receptors_micro_to_mono = enr_df(
        gset_receptors_micro_to_mono,
        genedbs,
        background_mono,
        out_folder,
    )
    common_paths_micro_to_mono = process_common_paths(
        enr_df_ligands_micro_to_mono, enr_df_receptors_micro_to_mono
    )
    common_paths_micro_to_mono.to_csv(out_folder / "Enriched_paths_micro_to_mono.csv")




def enrich_mono_micro_blood_brain(
    cohort, loadings_path, adata_br, adata_bl, adata_tf, design, out_folder, genedbs=GSEA_DBS
):

    n_celltypes_brain = adata_br.obs["cell_type.v2"].nunique()
    background_micro = get_expressed_genes(
        adata_br, "cell_type_coarse", "Microglia", expr_prop=0.1
    )
    background_mono = get_expressed_genes(
        adata_bl, "cell_type_coarse", "Monocytes", expr_prop=0.1
    )

    factors = c2c.io.load_tensor_factors(loadings_path / "Loadings.xlsx")
    thresh = find_threshold_loadings(factors, quantile=0.70)
    lr_loadings = factors["Ligand-Receptor Pairs"]

    # Monocytes to Microglia
    mo2mi_factor = CELLS_FACTOR_DICT[design][cohort]["Mono_to_Micro"]
    lr_loads_mono_to_micro = process_loadings(
        lr_loadings, mo2mi_factor, thresh
    ).reset_index()

    target_cells_mono_bl = [
        c for c in adata_bl.obs["cell_type.v2"].unique() if "mono" not in c
    ]

    all_mono_ligands = lr_loads_mono_to_micro["ligand"].unique().tolist()

    mono_ligands_bl_ = []
    for k in ["CD16", "CD14"]:
        ddf = ranked_genes_df(adata_bl, all_mono_ligands, key=k)
        ligands = filter_and_process_de(ddf, target_cells_mono_bl)
        mono_ligands_bl_.append(ligands)

    mono_ligands_bl = set(mono_ligands_bl_[0]) | set(mono_ligands_bl_[1])

    mo2mi_receiver_cells_ = factors["Receiver Cells"][mo2mi_factor] > 0.05
    mo2mi_receiver_cells = mo2mi_receiver_cells_[mo2mi_receiver_cells_].index.tolist()

    if any(
        "brain" in sub
        for sub in [g for g in mo2mi_receiver_cells if "Monocytes" not in g]
    ):
        target_cells_tf = [f.split("_")[0] for f in mo2mi_receiver_cells if "brain" in f]
        mono_ligands_ = []
        for k in ["CD16", "CD14"]:
            ddf = ranked_genes_df(adata_tf, mono_ligands_bl, key=k)
            ligands = filter_and_process_de(ddf, target_cells_tf)
            mono_ligands_.append(ligands)
        mono_ligands = list(set(mono_ligands_[0]) | set(mono_ligands_[1]))
    else:
        mono_ligands = list(mono_ligands_bl)

    all_receptors = get_unique_receptors(lr_loads_mono_to_micro, mono_ligands)

    micro_receptors = process_de_genes_celltype(
        adata_br, all_receptors, n_celltypes_brain, celltype="Microglia"
    )
    lr_loads_mono_to_micro_filtered = filter_lr(
        lr_loads_mono_to_micro, mono_ligands, micro_receptors
    )
    lr_loads_mono_to_micro_filtered.to_csv(
        loadings_path / "lr_loads_mono_to_micro_filtered.csv"
    )

    gset_ligands_mono_to_micro = list(
        lr_loads_mono_to_micro_filtered["ligand"].unique()
    )
    gset_receptors_mono_to_micro = list(
        lr_loads_mono_to_micro_filtered["receptor"].unique()
    )

    enr_df_ligands_mono_to_micro = enr_df(
        gset_ligands_mono_to_micro,
        genedbs,
        background_mono,
        out_folder,
    )
    enr_df_receptors_mono_to_micro = enr_df(
        gset_receptors_mono_to_micro,
        genedbs,
        background_micro,
        out_folder,
    )
    common_paths_mono_to_micro = process_common_paths(
        enr_df_ligands_mono_to_micro, enr_df_receptors_mono_to_micro
    )
    common_paths_mono_to_micro.to_csv(out_folder / "Enriched_paths_mono_to_micro.csv")


    # Microglia to Monocytes

    mi2mo_factor = CELLS_FACTOR_DICT[design][cohort]["Micro_to_Mono"]
    lr_loads_micro_to_mono = process_loadings(
        lr_loadings, mi2mo_factor, thresh
    ).reset_index()

    micro_ligands = process_de_genes_celltype(
        adata_br,
        lr_loads_micro_to_mono["ligand"].unique().tolist(),
        n_celltypes_brain,
        celltype="Microglia",
    )

    all_receptors = get_unique_receptors(lr_loads_micro_to_mono, micro_ligands)

    target_cells_mono_bl = [
        c for c in adata_bl.obs["cell_type.v2"].unique() if "mono" not in c
    ]

    mono_receptors_bl_ = []
    for k in ["CD16", "CD14"]:
        ddf = ranked_genes_df(adata_bl, all_receptors, key=k)
        receptors = filter_and_process_de(ddf, target_cells_mono_bl)
        mono_receptors_bl_.append(receptors)

    mono_receptors_bl = set(mono_receptors_bl_[0]) | set(mono_receptors_bl_[1])

    mi2mo_receiver_cells_ = factors["Receiver Cells"][mi2mo_factor] > 0.05
    mi2mo_receiver_cells = mi2mo_receiver_cells_[mi2mo_receiver_cells_].index.tolist()
    mi2mo_receiver_cells_no_mono = [
        g for g in mi2mo_receiver_cells if "Monocytes" not in g
    ]

    if any("brain" in sub for sub in mi2mo_receiver_cells_no_mono):
        target_cells_tf = [
            f.split("_")[0] for f in mi2mo_receiver_cells_no_mono if "brain" in f
        ]
        mono_receptors_ = []
        for k in ["CD16", "CD14"]:
            ddf = ranked_genes_df(adata_tf, mono_receptors_bl, key=k)
            receptors = filter_and_process_de(ddf, target_cells_tf)
            mono_receptors_.append(receptors)
        mono_receptors = list(set(mono_receptors_[0]) | set(mono_receptors_[1]))
    else:
        mono_receptors = list(mono_receptors_bl)

    lr_loads_micro_to_mono_filtered = filter_lr(
        lr_loads_micro_to_mono, micro_ligands, mono_receptors
    )
    lr_loads_micro_to_mono_filtered.to_csv(
        loadings_path / "lr_loads_micro_to_mono_filtered.csv"
    )
    gset_ligands_micro_to_mono = list(
        lr_loads_micro_to_mono_filtered["ligand"].unique()
    )
    gset_receptors_micro_to_mono = list(
        lr_loads_micro_to_mono_filtered["receptor"].unique()
    )

    enr_df_ligands_micro_to_mono = enr_df(
        gset_ligands_micro_to_mono,
        genedbs,
        background_micro,
        out_folder,
    )
    enr_df_receptors_micro_to_mono = enr_df(
        gset_receptors_micro_to_mono,
        genedbs,
        background_mono,
        out_folder,
    )
    common_paths_micro_to_mono = process_common_paths(
        enr_df_ligands_micro_to_mono, enr_df_receptors_micro_to_mono
    )
    common_paths_micro_to_mono.to_csv(out_folder / "Enriched_paths_micro_to_mono.csv")

