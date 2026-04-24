import argparse
from pathlib import Path

import decoupler as dc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from bokeh.palettes import Inferno256
from bokeh.plotting import show

from .adata_processing_utils import cell_type_mapping, filter_adata, load_brain_blood_data
from .plot_utils import interactive_embedding


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


def _extract_ora_df(pdata, db_name, source_label):
    """
    Extract ORA estimates from pdata into a tidy DataFrame.

    Parameters:
    pdata (AnnData): Pseudo-bulk AnnData with ora_estimate_{db_name} in obsm.
    db_name (str): Database key.
    source_label (str): Label added as 'Source' column (e.g. 'Brain (Cohort 1)').

    Returns:
    pd.DataFrame with columns: pathway scores + 'cell_t' + 'Source'.
    """
    df = pdata.obsm[f"ora_estimate_{db_name}"].copy()
    df.index = df.index.str.replace("Astrocytes", "Astrocyte", regex=False)
    df["cell_t"] = [x.split("_")[1] if "_" in x else x for x in df.index]
    df = df[~df.index.str.contains("nan", na=False)]
    df["Source"] = source_label
    return df


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
            adata.obs["cell_type"].astype(str).replace(cell_type_mapping).astype("category")
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



def plot_pseudobulk_ora(
    sources,
    db_name,
    paths,
    out_folder,
):
    """
    Boxplot + swarmplot of pseudo-bulk ORA estimates for selected pathways.

    Parameters:
    sources (list of tuple): Each entry is (label, pdata) where label is a
        key in SOURCE_PALETTE (e.g. 'Brain (Cohort 1)').
    db_name (str): Database key ('wikipathways', 'reactome', or 'go').
    paths (list): Pathway names to plot.
    out_folder (Path): Directory to save figures.
    """
    dfs = [_extract_ora_df(pdata, db_name, label) for label, pdata in sources]
    combined_df = pd.concat(dfs, ignore_index=True)

    available = [p for p in paths if p in combined_df.columns]
    if not available:
        print(f"No matching pathways found in {db_name} ORA results — skipping plot.")
        return

    palette = {k: v for k, v in SOURCE_PALETTE.items() if k in combined_df["Source"].unique()}

    for pathway in available:
        fig, ax = plt.subplots(figsize=(max(8, combined_df["cell_t"].nunique() * 1.5), 4))

        sns.boxplot(
            data=combined_df,
            x="cell_t",
            y=pathway,
            hue="Source",
            palette=palette,
            ax=ax,
        )
        sns.swarmplot(
            data=combined_df,
            x="cell_t",
            y=pathway,
            hue="Source",
            palette=palette,
            size=3,
            dodge=True,
            edgecolor="#444445",
            linewidth=1,
            alpha=0.8,
            legend=False,
            ax=ax,
        )

        cell_types = combined_df["cell_t"].unique()
        for i in range(1, len(cell_types)):
            ax.axvline(i - 0.5, color="black", linestyle="--", linewidth=0.3)

        ax.set_title(pathway, fontsize=10)
        ax.set_xlabel("")
        ax.set_ylabel("ORA score")
        ax.tick_params(axis="x", rotation=0, labelsize=12)
        ax.legend(loc="upper right", fontsize=12, bbox_to_anchor=(1.2, 0.5))

        plt.tight_layout()
        safe_name = pathway[:60].replace("/", "_")
        fig.savefig(out_folder / f"pseudobulk_ora_{db_name}_{safe_name}.pdf", bbox_inches="tight")
        plt.close(fig)




def main(project_path, cohort, threshold_targets):

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
            show(interactive_embedding(acts_br, label=path, palette_cont=Inferno256, title_font_size="20px"))
            show(interactive_embedding(acts_bl, label=path, palette_cont=Inferno256, title_font_size="20px"))

    # --- Pseudo-bulk ORA + boxplots ---
    cohort_label = "Cohort 1" if cohort == "cohort_1" else "Cohort 2"
    pdata_br, pdata_bl = ora_pseudobulk(project_path, cohort, threshold_targets)

    sources = [
        (f"Brain ({cohort_label})", pdata_br),
        (f"Blood ({cohort_label})", pdata_bl),
    ]
    for db_name, _ in DB_MAP:
        plot_pseudobulk_ora(sources, db_name, PATH_LISTS[db_name], out_folder)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p",
        "--project_path",
        help="Path where data are stored",
        default="/path/to/h5ad/",
    )
    parser.add_argument(
        "-b", "--cohort", help="cohort identifier", choices=["cohort_1", "cohort_2"]
    )
    parser.add_argument(
        "-t",
        "--threshold_targets",
        help="Top fraction of genes used as expressed input for ORA (default: 0.1 = top 10%%)",
        type=float,
        default=0.1,
    )

    args = parser.parse_args()
    main(args.project_path, args.cohort, args.threshold_targets)
