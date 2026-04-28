"""
Step 3: ORA and pseudo-bulk LR product.

Runs over-representation analysis (ORA) via decoupler on MSigDB gene sets (WikiPathways,
Reactome, GO) for brain and blood pseudo-bulk profiles. Produces per-pathway boxplot/
swarmplot figures and interactive UMAP embeddings coloured by ORA score.

When -f/--factorization_type is provided, also computes pseudo-bulk LR product heatmaps
for the filtered LR pairs written by step 2 (lr_loads_mono_to_micro_filtered.csv and
lr_loads_micro_to_mono_filtered.csv). Uses scran-normalized adatas so all LR genes are present.

Run as a standalone step:
    python -m tensorS2R.step3_downstream -p <path> -b <cohort> [-t <threshold>]
    python -m tensorS2R.step3_downstream -p <path> -b <cohort> -f <factorization_type> [-r <rank>] [-m <mode>] [-g <groupby>]
"""

import argparse
from pathlib import Path

import decoupler as dc
import numpy as np
import pandas as pd
import scanpy as sc
from bokeh.palettes import Inferno256
from bokeh.plotting import show

from .adata_processing_utils import (
    cell_type_mapping,
    filter_adata,
    load_brain_blood_data,
)
from .lr_loadings_utils import compute_prod_dict, get_mean_expression
from .plot_utils import interactive_embedding, plot_prod_dict, plot_pseudobulk_ora

msigdb = dc.get_resource("MSigDB")

msigdb_wiki = msigdb[msigdb["collection"] == "wikipathways"]
msigdb_reac = msigdb[msigdb["collection"] == "reactome_pathways"]
msigdb_go = msigdb[msigdb["collection"] == "go_biological_process"]

msigdb_wiki = msigdb_wiki[~msigdb_wiki.duplicated(["geneset", "genesymbol"])]
msigdb_reac = msigdb_reac[~msigdb_reac.duplicated(["geneset", "genesymbol"])]
msigdb_go = msigdb_go[~msigdb_go.duplicated(["geneset", "genesymbol"])]

DB_MAP = [
    ("wikipathways", msigdb_wiki),
    ("reactome", msigdb_reac),
    ("go", msigdb_go),
]

DEFAULT_WIKI_PATHS = [
    "WP_CELLS_AND_MOLECULES_INVOLVED_IN_LOCAL_ACUTE_INFLAMMATORY_RESPONSE",
    "WP_COMPLEMENT_SYSTEM_IN_NEURONAL_DEVELOPMENT_AND_PLASTICITY",
    "WP_TYROBP_CAUSAL_NETWORK_IN_MICROGLIA",
]

DEFAULT_REAC_PATHS = [
    "REACTOME_CD163_MEDIATING_AN_ANTI_INFLAMMATORY_RESPONSE",
    "REACTOME_ANTI_INFLAMMATORY_RESPONSE_FAVOURING_LEISHMANIA_PARASITE_INFECTION",
    "REACTOME_LEISHMANIA_INFECTION",
]

DEFAULT_GO_PATHS = [
    "GO_REGULATION_OF_INFLAMMATORY_RESPONSE",
    "GOBP_INFLAMMATORY_RESPONSE",
    "GOBP_POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION",
]

PATH_LISTS = {
    "wikipathways": DEFAULT_WIKI_PATHS,
    "reactome": DEFAULT_REAC_PATHS,
    "go": DEFAULT_GO_PATHS,
}

SOURCE_PALETTE = {
    "Brain (Cohort 1)": "#346C70",
    "Blood (Cohort 1)": "#B72234",
    "Brain (Cohort 2)": "#64BCC4",
    "Blood (Cohort 2)": "#FC4404",
}


def _run_ora_all_dbs(adata, threshold_targets):
    """Run ORA for all three databases, storing each result in obsm."""
    for name, msig_db in DB_MAP:
        dc.run_ora(
            mat=adata,
            net=msig_db,
            source="geneset",
            target="genesymbol",
            verbose=True,
            n_up=int(np.ceil(threshold_targets * adata.shape[1])),
            use_raw=False,
        )
        adata.obsm[f"ora_estimate_{name}"] = adata.obsm["ora_estimate"].copy()


def _get_acts_clipped(adata, db_name):
    """Extract ORA activity matrix, replacing inf with the finite max."""
    acts = dc.get_acts(adata, obsm_key=f"ora_estimate_{db_name}").copy()
    vals = acts.X.ravel()
    finite_max = np.nanmax(vals[np.isfinite(vals)])
    acts.X[~np.isfinite(acts.X)] = finite_max
    return acts


def ora_pseudobulk(project_path, cohort, threshold_targets):
    """
    Compute pseudo-bulk ORA for brain and blood across all three databases.

    Parameters:
    project_path (Path): Root project directory.
    cohort (str): 'cohort_1' or 'cohort_2'.
    threshold_targets (float): Top fraction of genes used as expressed input (e.g. 0.1).

    Returns:
    tuple: (pdata_br, pdata_bl) with ORA results stored in obsm.
    """
    adata_bl, adata_br, blood_adata_hvg, brain_adata_hvg = load_brain_blood_data(
        cohort, project_path
    )

    # Align full adatas to HVG-filtered indices and transfer cell type labels
    adata_bl = adata_bl[adata_bl.obs.index.isin(blood_adata_hvg.obs.index)].copy()
    adata_bl.obs["cell_type"] = blood_adata_hvg.obs.loc[adata_bl.obs.index, "cell_type"]

    adata_br = adata_br[adata_br.obs.index.isin(brain_adata_hvg.obs.index)].copy()
    adata_br.obs["cell_type"] = brain_adata_hvg.obs.loc[adata_br.obs.index, "cell_type"]

    for adata in [adata_br, adata_bl]:
        adata.obs["cell_type.v2"] = (
            adata.obs["cell_type"]
            .astype(str)
            .replace(cell_type_mapping)
            .astype("category")
        )

    pdata_br = dc.get_pseudobulk(
        adata_br,
        sample_col="pt",
        groups_col="cell_type.v2",
        layer="counts",
        mode="sum",
        min_cells=5,
        min_counts=10000,
    )
    pdata_bl = dc.get_pseudobulk(
        adata_bl,
        sample_col="pt",
        groups_col="cell_type.v2",
        layer="counts",
        mode="sum",
        min_cells=5,
        min_counts=10000,
    )

    _run_ora_all_dbs(pdata_br, threshold_targets)
    _run_ora_all_dbs(pdata_bl, threshold_targets)

    return pdata_br, pdata_bl


def compute_lr_products(
    project_path, cohort, loadings_folder, factorization_type, out_folder
):
    """
    Compute and save pseudo-bulk LR product heatmaps for filtered LR pairs from step 2.

    Loads scran-normalized adatas (full gene set, not HVG-filtered) to ensure all
    ligand/receptor genes are present. Microglia is the reference cell throughout.
    For brain_coarse the non-reference cell types come from the brain pdata; for
    micro_blood_coarse and brain_blood_coarse they come from the blood pdata.

    Parameters:
    project_path (Path): Root data directory.
    cohort (str): 'cohort_1' or 'cohort_2'.
    loadings_folder (Path): Folder with step 2 filtered LR CSV outputs.
    factorization_type (str): One of brain_coarse, micro_blood_coarse, brain_blood_coarse.
    out_folder (Path): Directory to save heatmap PDFs.
    """
    cohort_liana_path = Path(project_path) / cohort / "c2c_liana_outputs"
    celltype_col = "cell_type.v2"

    adata_br = sc.read(cohort_liana_path / "data_brain_scran_de_micro_mono.h5ad")
    pdata_br = dc.get_pseudobulk(
        adata_br,
        sample_col="pt",
        groups_col=celltype_col,
        layer="counts",
        mode="sum",
        min_cells=10,
        min_counts=10000,
    )
    df_micro = get_mean_expression(pdata_br, "Microglia", celltype_col)

    if factorization_type == "brain_coarse":
        # both sender and receiver are brain cell types
        pdata_nonref = pdata_br
    else:
        # monocytes (non-reference senders/receivers) are in blood
        adata_bl = sc.read(cohort_liana_path / "data_blood_scran_de_micro_mono.h5ad")
        pdata_nonref = dc.get_pseudobulk(
            adata_bl,
            sample_col="pt",
            groups_col=celltype_col,
            layer="counts",
            mode="sum",
            min_cells=10,
            min_counts=10000,
        )

    # Mono↔Micro: is_ligand=True → Microglia sends; is_ligand=False → Microglia receives
    lr_files = [
        ("lr_loads_mono_to_micro_filtered.csv", False, "lr_product_mono_to_micro.pdf"),
        ("lr_loads_micro_to_mono_filtered.csv", True, "lr_product_micro_to_mono.pdf"),
    ]
    for csv_name, is_ligand, out_name in lr_files:
        csv_path = Path(loadings_folder) / csv_name
        if not csv_path.exists():
            print(f"  {csv_name} not found — skipping.")
            continue
        lr_df = pd.read_csv(csv_path)
        lr_pairs = (lr_df["ligand"] + "^" + lr_df["receptor"]).tolist()
        prod_dict = compute_prod_dict(
            lr_pairs,
            pdata_nonref,
            df_micro,
            celltype_col,
            reference_cell="Microglia",
            is_ligand=is_ligand,
        )
        plot_prod_dict(prod_dict, str(out_folder / out_name))

    # Neuro↔Micro (brain_coarse only; factor number is embedded in the CSV filename)
    neuro_micro_pairs = [
        (
            "lr_loads_neuro_to_micro_filtered_*.csv",
            False,
            "lr_product_neuro_to_micro_{}.pdf",
        ),
        (
            "lr_loads_micro_to_neuro_filtered_*.csv",
            True,
            "lr_product_micro_to_neuro_{}.pdf",
        ),
    ]
    for pattern, is_ligand, out_template in neuro_micro_pairs:
        for csv_path in sorted(Path(loadings_folder).glob(pattern)):
            factor = csv_path.stem.rsplit("_", 1)[-1]
            lr_df = pd.read_csv(csv_path)
            lr_pairs = (lr_df["ligand"] + "^" + lr_df["receptor"]).tolist()
            prod_dict = compute_prod_dict(
                lr_pairs,
                pdata_br,
                df_micro,
                celltype_col,
                reference_cell="Microglia",
                is_ligand=is_ligand,
            )
            plot_prod_dict(prod_dict, str(out_folder / out_template.format(factor)))


def main(
    project_path,
    cohort,
    threshold_targets,
    loadings_folder=None,
    factorization_type=None,
):

    project_path = Path(project_path)
    out_folder = project_path / cohort / "ora_outputs"
    out_folder.mkdir(parents=True, exist_ok=True)

    # --- Cell-level ORA ---
    blood_adata, brain_adata, blood_adata_hvg, brain_adata_hvg = load_brain_blood_data(
        cohort, project_path
    )

    brain_adata_filt = filter_adata(brain_adata, brain_adata_hvg)
    blood_adata_filt = filter_adata(blood_adata, blood_adata_hvg)
    brain_adata_hvg_filt = filter_adata(brain_adata_hvg, brain_adata_hvg)
    blood_adata_hvg_filt = filter_adata(blood_adata_hvg, blood_adata_hvg)

    _run_ora_all_dbs(brain_adata_filt, threshold_targets)
    _run_ora_all_dbs(blood_adata_filt, threshold_targets)

    for db_name, _ in DB_MAP:
        brain_adata_hvg_filt.obsm[f"ora_estimate_{db_name}"] = brain_adata_filt.obsm[
            f"ora_estimate_{db_name}"
        ].copy()
        blood_adata_hvg_filt.obsm[f"ora_estimate_{db_name}"] = blood_adata_filt.obsm[
            f"ora_estimate_{db_name}"
        ].copy()

    for db_name, _ in DB_MAP:
        acts_br = _get_acts_clipped(brain_adata_hvg_filt, db_name)
        acts_bl = _get_acts_clipped(blood_adata_hvg_filt, db_name)

        for path in PATH_LISTS[db_name]:
            if path not in acts_br.var_names:
                continue
            acts_br.obs[path] = acts_br.obsm[f"ora_estimate_{db_name}"][path]
            acts_bl.obs[path] = acts_bl.obsm[f"ora_estimate_{db_name}"][path]
            show(
                interactive_embedding(
                    acts_br, label=path, palette_cont=Inferno256, title_font_size="20px"
                )
            )
            show(
                interactive_embedding(
                    acts_bl, label=path, palette_cont=Inferno256, title_font_size="20px"
                )
            )

    # --- Pseudo-bulk ORA + boxplots ---
    cohort_label = "Cohort 1" if cohort == "cohort_1" else "Cohort 2"
    pdata_br, pdata_bl = ora_pseudobulk(project_path, cohort, threshold_targets)

    sources = [
        (f"Brain ({cohort_label})", pdata_br),
        (f"Blood ({cohort_label})", pdata_bl),
    ]
    for db_name, _ in DB_MAP:
        plot_pseudobulk_ora(sources, db_name, PATH_LISTS[db_name], out_folder)

    if loadings_folder is not None and factorization_type is not None:
        compute_lr_products(
            project_path, cohort, loadings_folder, factorization_type, out_folder
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Step 3: ORA and pseudo-bulk LR product"
    )
    parser.add_argument(
        "-p", "--project_path", required=True, help="Root data directory"
    )
    parser.add_argument(
        "-b", "--cohort", required=True, choices=["cohort_1", "cohort_2"]
    )
    parser.add_argument(
        "-t",
        "--threshold_targets",
        type=float,
        default=0.1,
        help="Top fraction of genes for ORA (default: 0.1 = top 10%%)",
    )
    parser.add_argument(
        "-f",
        "--factorization_type",
        choices=["brain_coarse", "micro_blood_coarse", "brain_blood_coarse"],
        default=None,
        help="Factorization type; required for LR product heatmaps",
    )
    parser.add_argument(
        "-r",
        "--factorization_rank",
        type=int,
        default=20,
        help="Tensor rank used in step 1 (default: 20)",
    )
    parser.add_argument("-m", "--merging_mode", default="inner")
    parser.add_argument("-g", "--groupby", default="cell_type_coarse")
    args = parser.parse_args()

    project_path = Path(args.project_path)
    loadings_folder = None
    if args.factorization_type is not None:
        loadings_folder = (
            project_path
            / args.cohort
            / "c2c_liana_outputs"
            / args.factorization_type
            / f"rank_{args.factorization_rank}"
            / args.merging_mode
            / args.groupby
        )

    main(
        args.project_path,
        args.cohort,
        args.threshold_targets,
        loadings_folder,
        args.factorization_type,
    )
