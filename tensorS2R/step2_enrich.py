"""
Step 2: LR filtering and enrichment.

Dispatches to the appropriate enrichment functions based on cohort and factorization type:
  - All cohorts / all modes:         micro-mono enrichment
  - cohort_2 + brain_coarse:       + neuron-microglia enrichment
  - cohort_2 + brain_blood_coarse: + neuron-monocyte enrichment

Run as a standalone step (requires step 1 outputs to exist):
    python -m tensorS2R.step2_enrich -p <path> -b <cohort> -f <factorization_type> [options]

Factor numbers for each cohort and design are configured in CELLS_FACTOR_DICT and
NEURO_FACTOR_DICT at the top of this file. Update them when re-running factorization
produces a different factor ordering.
"""

import logging
import warnings
from pathlib import Path

import cell2cell as c2c
import decoupler as dc
import scanpy as sc

from .lr_loadings_utils import (
    compute_prod_dict,
    enr_df,
    filter_and_process_de,
    filter_lr,
    find_threshold_loadings,
    get_expressed_genes,
    get_mean_expression,
    get_unique_receptors,
    plot_prod_dict,
    process_common_paths,
    process_de_genes_celltype,
    process_loadings,
    ranked_genes_df,
)
from .plot_utils import plot_pathway_network

logging.disable(logging.WARNING)
warnings.filterwarnings("ignore")

GSEA_DBS = ["WikiPathway_2023_Human", "Reactome_2022", "GO_Biological_Process_2023"]

# Factor numbers for microglia-monocyte communication per cohort and design.
CELLS_FACTOR_DICT = {
    "micro_blood_coarse": {
        "cohort_1": {"Micro_to_Mono": "Factor 8",  "Mono_to_Micro": "Factor 3"},
        "cohort_2": {"Micro_to_Mono": "Factor 5",  "Mono_to_Micro": "Factor 3"},
    },
    "brain_blood_coarse": {
        "cohort_1": {"Micro_to_Mono": "Factor 12", "Mono_to_Micro": "Factor 15"},
        "cohort_2": {"Micro_to_Mono": "Factor 6",  "Mono_to_Micro": "Factor 17"},
    },
    "brain_coarse": {
        "cohort_1": {"Micro_to_Mono": "Factor 4",  "Mono_to_Micro": "Factor 3"},
    },
}

# Factor numbers for neuron-focused communication (Cohort 2 only).
NEURO_FACTOR_DICT = {
    "brain_coarse": {
        "cohort_2": {"Neuro_to_Micro": "Factor 9", "Micro_to_Neuro": "Factor 10"},
    },
    "brain_blood_coarse": {
        "cohort_2": {"Neuro_to_Mono": "Factor 6", "Mono_to_Neuro": "Factor 16"},
    },
}

_BRAIN_CUSTOM_PALETTE = {
    "Astrocyte": "#5D3534",
    "Epithelial": "#FADF81",
    "GABAergic Neuron": "#7BB6F9",
    "Glutamatergic Neuron": "#C2D3E7",
    "Neuron-low count": "#EC68B1",
    "OPC": "#EF8B72",
    "Oligo": "#E8B394",
}




def _filter_mono_genes(adata, genes, target_cells, filter_fn=filter_and_process_de):
    """Filter genes by CD16/CD14 DE vs. target_cells; return the union across both keys."""
    result = []
    for k in ["CD16", "CD14"]:
        ddf = ranked_genes_df(adata, genes, key=k)
        result.append(filter_fn(ddf, target_cells))
    return list(set(result[0]) | set(result[1]))


def _filter_neuro_genes(adata_br, genes, target_cells, filter_fn=filter_and_process_de):
    """Filter genes by Glutamatergic/GABAergic/Neuron-low DE vs. target_cells; return the union."""
    result = []
    for k in ["Gluta", "GABA", "low"]:
        ddf = ranked_genes_df(adata_br, genes, key=k)
        result.append(filter_fn(ddf, target_cells))
    return list(set(result[0]) | set(result[1]) | set(result[2]))


def _enrich_and_plot(filtered_df, background_lig, background_rec, title,
                     enrichr_dir, csv_path, plot_path, genedbs):
    """Run Enrichr on filtered LR sets, find common pathways, and save outputs."""
    enr_lig = enr_df(list(filtered_df["ligand"].unique()), genedbs, background_lig, str(enrichr_dir))
    enr_rec = enr_df(list(filtered_df["receptor"].unique()), genedbs, background_rec, str(enrichr_dir))
    common = process_common_paths(enr_lig, enr_rec)
    common.to_csv(csv_path)
    plot_pathway_network(common, title, plot_path)
    return common


def _process_de_strict(ddf, de_count, logfc_threshold=0.5, pval_threshold=0.05):
    """Like filter_and_process_de but condition2 also requires p >= threshold (stricter)."""
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


def _filter_and_process_de_strict(df, cell_types, logfc=0.5, pval=0.05):
    """Strict version of filter_and_process_de; used for neuron-monocyte analyses."""
    filtered_df = df[df["group"].isin(cell_types)]
    return _process_de_strict(filtered_df, de_count=len(cell_types), logfc_threshold=logfc, pval_threshold=pval)


# Micro-Mono Enrichment (all cohorts, all factorization types)

def enrich_mono_micro_brain(cohort, loadings_path, adata_br, genedbs, design, out_folder):
    """LR filtering and enrichment for Mono↔Micro using brain-only data (brain_coarse)."""
    background_micro = get_expressed_genes(adata_br, "cell_type_coarse", "Microglia", expr_prop=0.1)
    background_mono = get_expressed_genes(adata_br, "cell_type_coarse", "Monocytes", expr_prop=0.1)

    factors = c2c.io.load_tensor_factors(loadings_path / "Loadings.xlsx")
    lr_loadings = factors["Ligand-Receptor Pairs"]
    thresh = find_threshold_loadings(factors, quantile=0.70)
    target_cells_mono = [c for c in adata_br.obs["cell_type.v2"].unique() if "mono" not in c]

    ## Monocytes --> Microglia
    mo2mi_factor = CELLS_FACTOR_DICT[design][cohort]["Mono_to_Micro"]
    lr_mono_to_micro = process_loadings(lr_loadings, mo2mi_factor, thresh)
    top_receivers = factors["Receiver Cells"][mo2mi_factor] > thresh

    # receptors: preferentially expressed in Microglia vs. top-loading receiver cell types
    micro_receptors = filter_and_process_de(
        ranked_genes_df(adata_br, lr_mono_to_micro["receptor"].unique(), key="Microglia"),
        top_receivers[top_receivers].index.tolist(),
    )
    # ligands: preferentially expressed in monocytes (CD16 or CD14) vs. all non-monocyte types
    candidate_ligands = lr_mono_to_micro[lr_mono_to_micro["receptor"].isin(micro_receptors)]["ligand"].unique()
    mono_ligands = _filter_mono_genes(adata_br, candidate_ligands, target_cells_mono)

    lr_mono_to_micro_filt = filter_lr(lr_mono_to_micro.reset_index(), mono_ligands, micro_receptors)
    lr_mono_to_micro_filt.to_csv(loadings_path / "lr_loads_mono_to_micro_filtered.csv")
    _enrich_and_plot(
        lr_mono_to_micro_filt, background_mono, background_micro,
        "Monocytes -> Microglia", out_folder,
        out_folder / "Enriched_paths_mono_to_micro.csv",
        out_folder / "Net_paths_mono_to_micro.png",
        genedbs,
    )

    ## Microglia --> Monocytes
    mi2mo_factor = CELLS_FACTOR_DICT[design][cohort]["Micro_to_Mono"]
    lr_micro_to_mono = process_loadings(lr_loadings, mi2mo_factor, thresh)
    top_senders = factors["Sender Cells"][mi2mo_factor] > thresh

    # ligands: preferentially expressed in Microglia vs. top-loading sender cell types
    micro_ligands = filter_and_process_de(
        ranked_genes_df(adata_br, lr_micro_to_mono["ligand"].unique(), key="Microglia"),
        top_senders[top_senders].index.tolist(),
    )
    # receptors: preferentially expressed in monocytes (CD16 or CD14)
    candidate_receptors = lr_micro_to_mono[lr_micro_to_mono["ligand"].isin(micro_ligands)]["receptor"].unique()
    mono_receptors = _filter_mono_genes(adata_br, candidate_receptors, target_cells_mono)

    lr_micro_to_mono_filt = filter_lr(lr_micro_to_mono.reset_index(), micro_ligands, mono_receptors)
    lr_micro_to_mono_filt.to_csv(loadings_path / "lr_loads_micro_to_mono_filtered.csv")
    _enrich_and_plot(
        lr_micro_to_mono_filt, background_micro, background_mono,
        "Microglia -> Monocytes", out_folder,
        out_folder / "Enriched_paths_micro_to_mono.csv",
        out_folder / "Net_paths_micro_to_mono.png",
        genedbs,
    )


def enrich_mono_micro_blood_micro(
    cohort, loadings_path, adata_br, adata_bl, adata_tf, design, out_folder, genedbs=GSEA_DBS
):
    """LR filtering and enrichment for Mono↔Micro using blood+microglia data (micro_blood_coarse)."""
    n_celltypes_brain = adata_br.obs["cell_type.v2"].nunique()
    background_micro = get_expressed_genes(adata_br, "cell_type_coarse", "Microglia", expr_prop=0.1)
    background_mono = get_expressed_genes(adata_bl, "cell_type_coarse", "Monocytes", expr_prop=0.1)

    factors = c2c.io.load_tensor_factors(loadings_path / "Loadings.xlsx")
    lr_loadings = factors["Ligand-Receptor Pairs"]
    thresh = find_threshold_loadings(factors, quantile=0.70)

    target_cells_mono_bl = [c for c in adata_bl.obs["cell_type.v2"].unique() if "mono" not in c]
    target_cells_mono_tf = [c for c in adata_tf.obs["cell_type.v2"].unique() if "mono" not in c]

    ## Monocytes --> Microglia
    mo2mi_factor = CELLS_FACTOR_DICT[design][cohort]["Mono_to_Micro"]
    lr_mono_to_micro = process_loadings(lr_loadings, mo2mi_factor, thresh).reset_index()

    # ligands: must be DE in monocytes in BOTH blood adata and combined adata (intersection)
    all_ligands = lr_mono_to_micro["ligand"].unique().tolist()
    mono_ligands = list(
        set(_filter_mono_genes(adata_bl, all_ligands, target_cells_mono_bl)) &
        set(_filter_mono_genes(adata_tf, all_ligands, target_cells_mono_tf))
    )
    # receptors: must be DE in Microglia vs. all other brain cell types
    micro_receptors = process_de_genes_celltype(
        adata_br, get_unique_receptors(lr_mono_to_micro, mono_ligands), n_celltypes_brain, celltype="Microglia"
    )
    lr_mono_to_micro_filt = filter_lr(lr_mono_to_micro, mono_ligands, micro_receptors)
    lr_mono_to_micro_filt.to_csv(loadings_path / "lr_loads_mono_to_micro_filtered.csv")
    _enrich_and_plot(
        lr_mono_to_micro_filt, background_mono, background_micro,
        "Monocytes -> Microglia", out_folder,
        out_folder / "Enriched_paths_mono_to_micro.csv",
        out_folder / "Net_paths_mono_to_micro.png",
        genedbs,
    )

    ## Microglia --> Monocytes
    mi2mo_factor = CELLS_FACTOR_DICT[design][cohort]["Micro_to_Mono"]
    lr_micro_to_mono = process_loadings(lr_loadings, mi2mo_factor, thresh).reset_index()

    # ligands: must be DE in Microglia vs. all other brain cell types
    micro_ligands = process_de_genes_celltype(
        adata_br, lr_micro_to_mono["ligand"].unique().tolist(), n_celltypes_brain, celltype="Microglia"
    )
    # receptors: must be DE in monocytes in BOTH blood adata and combined adata (intersection)
    all_receptors = get_unique_receptors(lr_micro_to_mono, micro_ligands)
    mono_receptors = list(
        set(_filter_mono_genes(adata_bl, all_receptors, target_cells_mono_bl)) &
        set(_filter_mono_genes(adata_tf, all_receptors, target_cells_mono_tf))
    )
    lr_micro_to_mono_filt = filter_lr(lr_micro_to_mono, micro_ligands, mono_receptors)
    lr_micro_to_mono_filt.to_csv(loadings_path / "lr_loads_micro_to_mono_filtered.csv")
    _enrich_and_plot(
        lr_micro_to_mono_filt, background_micro, background_mono,
        "Microglia -> Monocytes", out_folder,
        out_folder / "Enriched_paths_micro_to_mono.csv",
        out_folder / "Net_paths_micro_to_mono.png",
        genedbs,
    )


def enrich_mono_micro_blood_brain(
    cohort, loadings_path, adata_br, adata_bl, adata_tf, design, out_folder, genedbs=GSEA_DBS
):
    """LR filtering and enrichment for Mono↔Micro using full blood+brain data (brain_blood_coarse)."""
    n_celltypes_brain = adata_br.obs["cell_type.v2"].nunique()
    background_micro = get_expressed_genes(adata_br, "cell_type_coarse", "Microglia", expr_prop=0.1)
    background_mono = get_expressed_genes(adata_bl, "cell_type_coarse", "Monocytes", expr_prop=0.1)

    factors = c2c.io.load_tensor_factors(loadings_path / "Loadings.xlsx")
    lr_loadings = factors["Ligand-Receptor Pairs"]
    thresh = find_threshold_loadings(factors, quantile=0.70)

    target_cells_mono_bl = [c for c in adata_bl.obs["cell_type.v2"].unique() if "mono" not in c]

    ## Monocytes --> Microglia
    mo2mi_factor = CELLS_FACTOR_DICT[design][cohort]["Mono_to_Micro"]
    lr_mono_to_micro = process_loadings(lr_loadings, mo2mi_factor, thresh).reset_index()

    # First pass: ligands DE in monocytes in the blood adata
    all_ligands = lr_mono_to_micro["ligand"].unique().tolist()
    mono_ligands_bl = set(_filter_mono_genes(adata_bl, all_ligands, target_cells_mono_bl))

    # Second pass: if brain cell types are top receivers, also confirm expression in
    # the brain context using the combined adata (keeps only genes that survive both filters)
    mo2mi_top_receivers_ = factors["Receiver Cells"][mo2mi_factor] > 0.05
    mo2mi_top_receivers = mo2mi_top_receivers_[mo2mi_top_receivers_].index.tolist()
    non_mono_brain_receivers = [r for r in mo2mi_top_receivers if "Monocytes" not in r and "brain" in r]
    if non_mono_brain_receivers:
        target_cells_tf = [r.split("_")[0] for r in non_mono_brain_receivers]
        mono_ligands = set(_filter_mono_genes(adata_tf, list(mono_ligands_bl), target_cells_tf))
    else:
        mono_ligands = mono_ligands_bl
    mono_ligands = list(mono_ligands)

    # Receptors: DE in Microglia vs. all other brain cell types
    micro_receptors = process_de_genes_celltype(
        adata_br, get_unique_receptors(lr_mono_to_micro, mono_ligands), n_celltypes_brain, celltype="Microglia"
    )
    lr_mono_to_micro_filt = filter_lr(lr_mono_to_micro, mono_ligands, micro_receptors)
    lr_mono_to_micro_filt.to_csv(loadings_path / "lr_loads_mono_to_micro_filtered.csv")
    _enrich_and_plot(
        lr_mono_to_micro_filt, background_mono, background_micro,
        "Monocytes -> Microglia", out_folder,
        out_folder / "Enriched_paths_mono_to_micro.csv",
        out_folder / "Net_paths_mono_to_micro.png",
        genedbs,
    )

    ## Microglia --> Monocytes
    mi2mo_factor = CELLS_FACTOR_DICT[design][cohort]["Micro_to_Mono"]
    lr_micro_to_mono = process_loadings(lr_loadings, mi2mo_factor, thresh).reset_index()

    # Ligands: DE in Microglia vs. all other brain cell types
    micro_ligands = process_de_genes_celltype(
        adata_br, lr_micro_to_mono["ligand"].unique().tolist(), n_celltypes_brain, celltype="Microglia"
    )

    # First pass: receptors DE in monocytes in the blood adata
    all_receptors = get_unique_receptors(lr_micro_to_mono, micro_ligands)
    mono_receptors_bl = set(_filter_mono_genes(adata_bl, all_receptors, target_cells_mono_bl))

    # Second pass: if brain cell types are top receivers, also confirm expression
    # in the brain context using the combined adata
    mi2mo_top_receivers_ = factors["Receiver Cells"][mi2mo_factor] > 0.05
    mi2mo_top_receivers = mi2mo_top_receivers_[mi2mo_top_receivers_].index.tolist()
    non_mono_brain_receivers = [r for r in mi2mo_top_receivers if "Monocytes" not in r and "brain" in r]
    if non_mono_brain_receivers:
        target_cells_tf = [r.split("_")[0] for r in non_mono_brain_receivers]
        mono_receptors = set(_filter_mono_genes(adata_tf, list(mono_receptors_bl), target_cells_tf))
    else:
        mono_receptors = mono_receptors_bl
    mono_receptors = list(mono_receptors)

    lr_micro_to_mono_filt = filter_lr(lr_micro_to_mono, micro_ligands, mono_receptors)
    lr_micro_to_mono_filt.to_csv(loadings_path / "lr_loads_micro_to_mono_filtered.csv")
    _enrich_and_plot(
        lr_micro_to_mono_filt, background_micro, background_mono,
        "Microglia -> Monocytes", out_folder,
        out_folder / "Enriched_paths_micro_to_mono.csv",
        out_folder / "Net_paths_micro_to_mono.png",
        genedbs,
    )


# Neuro-Micro Enrichment (Cohort 2, brain_coarse)

def enrich_neuro_micro(
    loadings_folder,
    adata_br,
    out_folder,
    factor_n2mi,
    factor_mi2n,
    genedbs=GSEA_DBS,
):
    """
    LR filtering, enrichment, and pseudo-bulk product plots for Neuron↔Microglia.

    factor_n2mi: numeric factor label for Neuron → Microglia (e.g. "9")
    factor_mi2n: numeric factor label for Microglia → Neuron (e.g. "10")
    """
    adata_br.obs["cell_type_coarse"] = (
        adata_br.obs["cell_type.v2"].astype(str).replace({
            "GABAergic Neuron": "Neurons",
            "Glutamatergic Neuron": "Neurons",
            "Neuron-low count": "Neurons",
        })
    )
    for ref, key in [
        ("GABAergic Neuron", "GABA"),
        ("Glutamatergic Neuron", "Gluta"),
        ("Neuron-low count", "low"),
    ]:
        sc.tl.rank_genes_groups(
            adata_br, "cell_type.v2", reference=ref, method="wilcoxon",
            layer="scran_normalization", pts=True, key_added=key,
        )

    background_micro = get_expressed_genes(adata_br, "cell_type_coarse", "Microglia", expr_prop=0.1)
    background_neuro = get_expressed_genes(adata_br, "cell_type_coarse", "Neurons", expr_prop=0.1)

    pdata_br = dc.get_pseudobulk(
        adata_br, sample_col="pt", groups_col="cell_type.v2",
        layer="counts", mode="sum", min_cells=10, min_counts=10000,
    )
    df_micro = get_mean_expression(pdata_br, "Microglia", "cell_type.v2")

    factors = c2c.io.load_tensor_factors(loadings_folder / "Loadings.xlsx")
    lr_loadings = factors["Ligand-Receptor Pairs"]
    thresh = find_threshold_loadings(factors, quantile=0.70)

    ## Neurons --> Microglia
    lr_n2mi = process_loadings(lr_loadings, f"Factor {factor_n2mi}", thresh)

    # neuron ligands: DE vs. OPC (closest non-neuronal reference in brain_coarse)
    neuro_ligands = _filter_neuro_genes(adata_br, lr_n2mi["ligand"].unique(), target_cells=["OPC"])
    # microglia receptors: DE vs. other glial cell types
    candidate_receptors = lr_n2mi[lr_n2mi["ligand"].isin(neuro_ligands)]["receptor"].unique()
    micro_receptors = filter_and_process_de(
        ranked_genes_df(adata_br, candidate_receptors, key="Microglia"),
        ["Astrocyte", "Epithelial", "OPC", "Oligo"],
    )

    lr_n2mi_filt = filter_lr(lr_n2mi.reset_index(), neuro_ligands, micro_receptors)
    lr_n2mi_filt.to_csv(loadings_folder / f"lr_loads_neuro_to_micro_filtered_{factor_n2mi}.csv")
    _enrich_and_plot(
        lr_n2mi_filt, background_neuro, background_micro,
        "Neurons -> Microglia", out_folder,
        loadings_folder / f"Enriched_paths_neuro_to_micro_{factor_n2mi}.csv",
        loadings_folder / f"Net_paths_neuro_to_micro_{factor_n2mi}.png",
        genedbs,
    )
    prod_dict_n2mi = compute_prod_dict(
        (lr_n2mi_filt["ligand"] + "^" + lr_n2mi_filt["receptor"]).tolist(),
        pdata_br, df_micro, "cell_type.v2", reference_cell="Microglia", is_ligand=False,
    )
    plot_prod_dict(prod_dict_n2mi, str(loadings_folder / f"lr_product_neuro_to_micro_{factor_n2mi}.pdf"))

    ## Microglia --> Neurons
    lr_mi2n = process_loadings(lr_loadings, f"Factor {factor_mi2n}", thresh)
    target_cells_neuro_all = [c for c in adata_br.obs["cell_type.v2"].unique() if "Neuron" not in c]

    # microglia ligands: not DE-filtered on the sender side
    micro_ligands = lr_mi2n["ligand"].unique()
    # neuron receptors: DE vs. all non-neuronal brain types across all three references
    candidate_receptors = lr_mi2n[lr_mi2n["ligand"].isin(micro_ligands)]["receptor"].unique()
    neuro_receptors = _filter_neuro_genes(adata_br, candidate_receptors, target_cells_neuro_all)

    lr_mi2n_filt = filter_lr(lr_mi2n.reset_index(), micro_ligands, neuro_receptors)
    lr_mi2n_filt.to_csv(loadings_folder / f"lr_loads_micro_to_neuro_filtered_{factor_mi2n}.csv")
    _enrich_and_plot(
        lr_mi2n_filt, background_micro, background_neuro,
        "Microglia -> Neurons", out_folder,
        loadings_folder / f"Enriched_paths_micro_to_neuro_{factor_mi2n}.csv",
        loadings_folder / f"Net_paths_micro_to_neuro_{factor_mi2n}.png",
        genedbs,
    )
    prod_dict_mi2n = compute_prod_dict(
        (lr_mi2n_filt["ligand"] + "^" + lr_mi2n_filt["receptor"]).tolist(),
        pdata_br, df_micro, "cell_type.v2", reference_cell="Microglia", is_ligand=True,
    )
    plot_prod_dict(prod_dict_mi2n, str(loadings_folder / f"lr_product_micro_to_neuro_{factor_mi2n}.pdf"))


# Neuro-Mono Enrichment (Cohort 2, brain_blood_coarse)

def enrich_neuro_mono(
    loadings_folder,
    adata_br,
    adata_bl,
    out_folder,
    factor_n2mo,
    factor_mo2n,
    genedbs=GSEA_DBS,
):
    """
    LR filtering and enrichment for Neuron↔Monocyte communication.

    factor_n2mo: numeric factor label for Neuron → Monocyte (e.g. "6")
    factor_mo2n: numeric factor label for Monocyte → Neuron (e.g. "16")

    Uses stricter DE filtering (_filter_and_process_de_strict) than the micro-mono analyses:
    condition2 additionally requires p >= threshold to confirm non-significance.
    """
    adata_br.obs["cell_type_coarse"] = (
        adata_br.obs["cell_type.v2"].astype(str).replace({
            "GABAergic Neuron": "Neurons",
            "Glutamatergic Neuron": "Neurons",
            "Neuron-low count": "Neurons",
        })
    )
    adata_bl.obs["cell_type_coarse"] = adata_bl.obs["cell_type.v2"].astype(str).apply(
        lambda x: "Monocytes" if "mono" in x.lower() else x
    )
    for ref, key in [
        ("GABAergic Neuron", "GABA"),
        ("Glutamatergic Neuron", "Gluta"),
        ("Neuron-low count", "low"),
    ]:
        sc.tl.rank_genes_groups(
            adata_br, "cell_type.v2", reference=ref, method="wilcoxon",
            layer="scran_normalization", pts=True, key_added=key,
        )

    background_mono = get_expressed_genes(adata_bl, "cell_type_coarse", "Monocytes", expr_prop=0.1)
    background_neuro = get_expressed_genes(adata_br, "cell_type_coarse", "Neurons", expr_prop=0.1)

    factors = c2c.io.load_tensor_factors(loadings_folder / "Loadings.xlsx")
    lr_loadings = factors["Ligand-Receptor Pairs"]
    thresh = find_threshold_loadings(factors, quantile=0.70)

    target_cells_neuro = [c for c in adata_br.obs["cell_type.v2"].unique() if "Neuron" not in c]
    target_cells_mono = [c for c in adata_bl.obs["cell_type.v2"].unique() if "mono" not in c]

    ## Neurons --> Monocytes
    lr_n2mo = process_loadings(lr_loadings, f"Factor {factor_n2mo}", thresh)

    # neuron ligands: strict DE vs. all non-neuronal brain types across all three references
    neuro_ligands = _filter_neuro_genes(
        adata_br, lr_n2mo["ligand"].unique(), target_cells_neuro,
        filter_fn=_filter_and_process_de_strict,
    )
    # monocyte receptors: strict DE in CD16 or CD14 vs. all other blood types
    candidate_receptors = lr_n2mo[lr_n2mo["ligand"].isin(neuro_ligands)]["receptor"].unique().tolist()
    mono_receptors = _filter_mono_genes(
        adata_bl, candidate_receptors, target_cells_mono,
        filter_fn=_filter_and_process_de_strict,
    )

    lr_n2mo_filt = filter_lr(lr_n2mo.reset_index(), neuro_ligands, mono_receptors)
    lr_n2mo_filt.to_csv(loadings_folder / "lr_loads_neuro_to_mono_filtered.csv")
    _enrich_and_plot(
        lr_n2mo_filt, background_neuro, background_mono,
        "Neurons -> Monocytes", out_folder,
        loadings_folder / "Enriched_paths_neuro_to_mono.csv",
        loadings_folder / "Net_paths_neuro_to_mono.png",
        genedbs,
    )

    ## Monocytes --> Neurons
    lr_mo2n = process_loadings(lr_loadings, f"Factor {factor_mo2n}", thresh)

    # monocyte ligands: strict DE in CD16 or CD14 vs. all other blood types
    mono_ligands = _filter_mono_genes(
        adata_bl, lr_mo2n["ligand"].unique(), target_cells_mono,
        filter_fn=_filter_and_process_de_strict,
    )
    # neuron receptors: strict DE vs. all non-neuronal brain types across all three references
    candidate_receptors = lr_mo2n[lr_mo2n["ligand"].isin(mono_ligands)]["receptor"].unique().tolist()
    neuro_receptors = _filter_neuro_genes(
        adata_br, candidate_receptors, target_cells_neuro,
        filter_fn=_filter_and_process_de_strict,
    )

    lr_mo2n_filt = filter_lr(lr_mo2n.reset_index(), mono_ligands, neuro_receptors)
    lr_mo2n_filt.to_csv(loadings_folder / "lr_loads_mono_to_neuro_filtered.csv")
    _enrich_and_plot(
        lr_mo2n_filt, background_mono, background_neuro,
        "Monocytes -> Neurons", out_folder,
        loadings_folder / "Enriched_paths_mono_to_neuro.csv",
        loadings_folder / "Net_paths_mono_to_neuro.png",
        genedbs,
    )



def run_enrichment(
    cohort,
    factorization_type,
    loadings_folder,
    adata_br,
    adata_bl=None,
    adata_tf=None,
    out_folder=None,
    genedbs=GSEA_DBS,
):
    """
    Run all enrichment analyses appropriate for the given cohort and factorization type.

    Step 2a — Micro-Mono (all cohorts, all factorization types):
        Filters LR pairs and runs pathway enrichment for Monocyte↔Microglia factors.

    Step 2b — Neuron analyses (Cohort 2 only):
        brain_coarse:       + Neuron↔Microglia LR filtering, enrichment, and product heatmaps.
        brain_blood_coarse: + Neuron↔Monocyte LR filtering and enrichment.
    """
    if out_folder is None:
        out_folder = loadings_folder / "enrichr_results"
    out_folder = Path(out_folder)
    out_folder.mkdir(parents=True, exist_ok=True)

    # Step 2a: Micro-Mono (skip if cohort not in CELLS_FACTOR_DICT for this design)
    if cohort in CELLS_FACTOR_DICT.get(factorization_type, {}):
        if factorization_type == "brain_coarse":
            enrich_mono_micro_brain(cohort, loadings_folder, adata_br, genedbs, factorization_type, out_folder)
        elif factorization_type == "micro_blood_coarse":
            enrich_mono_micro_blood_micro(cohort, loadings_folder, adata_br, adata_bl, adata_tf, factorization_type, out_folder, genedbs)
        elif factorization_type == "brain_blood_coarse":
            enrich_mono_micro_blood_brain(cohort, loadings_folder, adata_br, adata_bl, adata_tf, factorization_type, out_folder, genedbs)
    else:
        print(f"No CELLS_FACTOR_DICT entry for {cohort} + {factorization_type} — skipping micro-mono enrichment.")

    # Step 2b: Neuro analyses (Cohort 2 only)
    if cohort == "cohort_2":
        if factorization_type == "brain_coarse" and "cohort_2" in NEURO_FACTOR_DICT.get("brain_coarse", {}):
            nf = NEURO_FACTOR_DICT["brain_coarse"]["cohort_2"]
            enrich_neuro_micro(
                loadings_folder, adata_br, out_folder,
                factor_n2mi=nf["Neuro_to_Micro"].replace("Factor ", ""),
                factor_mi2n=nf["Micro_to_Neuro"].replace("Factor ", ""),
                genedbs=genedbs,
            )
        elif factorization_type == "brain_blood_coarse" and "cohort_2" in NEURO_FACTOR_DICT.get("brain_blood_coarse", {}):
            nf = NEURO_FACTOR_DICT["brain_blood_coarse"]["cohort_2"]
            enrich_neuro_mono(
                loadings_folder, adata_br, adata_bl, out_folder,
                factor_n2mo=nf["Neuro_to_Mono"].replace("Factor ", ""),
                factor_mo2n=nf["Mono_to_Neuro"].replace("Factor ", ""),
                genedbs=genedbs,
            )



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Step 2: LR filtering and enrichment")
    parser.add_argument("-p", "--project_path", required=True, help="Root data directory")
    parser.add_argument("-b", "--cohort", required=True, choices=["cohort_1", "cohort_2"])
    parser.add_argument(
        "-f", "--factorization_type", required=True,
        choices=["brain_coarse", "micro_blood_coarse", "brain_blood_coarse"],
    )
    parser.add_argument("-r", "--factorization_rank", type=int, default=20)
    parser.add_argument("-m", "--merging_mode", default="inner")
    parser.add_argument("-g", "--groupby", default="cell_type_coarse")
    args = parser.parse_args()

    project_path = Path(args.project_path)
    cohort_liana_path = project_path / args.cohort / "c2c_liana_outputs"
    loadings_folder = (
        cohort_liana_path
        / args.factorization_type
        / f"rank_{args.factorization_rank}"
        / args.merging_mode
        / args.groupby
    )

    try:
        adata_br = sc.read(cohort_liana_path / "data_brain_scran_de_micro_mono.h5ad")
        adata_bl = sc.read(cohort_liana_path / "data_blood_scran_de_micro_mono.h5ad")
    except FileNotFoundError:
        from .adata_processing_utils import prepare_scran_de_data_micro_mono
        adata_bl, adata_br = prepare_scran_de_data_micro_mono(
            args.cohort, project_path, "cell_type.v2", "cell_type.v2"
        )

    adata_tf = None
    if args.factorization_type == "micro_blood_coarse":
        adata_tf = sc.read(cohort_liana_path / "blood_microglia_tf.h5ad")
    elif args.factorization_type == "brain_blood_coarse":
        adata_tf = sc.read(cohort_liana_path / "blood_brain_tf.h5ad")

    run_enrichment(args.cohort, args.factorization_type, loadings_folder, adata_br, adata_bl, adata_tf)
