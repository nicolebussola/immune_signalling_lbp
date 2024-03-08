import logging
import os

import anndata2ri
import numpy as np
import pandas as pd
import rpy2.robjects as ro
import scanpy as sc
from bokeh.models import TabPanel, Tabs
from bokeh.plotting import output_file, show
from rich.console import Console
from rich.logging import RichHandler
from rich.table import Table
from rpy2.robjects.packages import importr

from ..labels import QC_LABELS_SAMPLE
from ..utils import QC_metrics_UMAP_plot, interactive_embedding
from .doublets_detection import (doubletdetection_method, scdblfinder, scds,
                                 scrublet, solo)

logging.getLogger().setLevel(logging.INFO)
FORMAT = "%(message)s"
logging.basicConfig(
    level="INFO", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
)

log = logging.getLogger("rich")

anndata2ri.activate()
console = Console()

sc.settings.verbosity = 0

importr("Seurat")
importr("scater")
importr("BiocParallel")
importr("scry")


np.random.seed(42)


def doublet_methods_decision(adata):
    if adata.shape[0] < 2000:
        return ["scrublet", "scds", "scdblfinder"]  ### Missing doubletfinder!!
    elif adata.shape[0] >= 2000 and adata.shape[0] <= 10000:
        return ["scds", "scdblfinder"]
    elif adata.shape[0] > 10000 and adata.shape[0] < 20000:
        return ["scrublet", "scds", "scdblfinder", "doubletdetection"]
    else:
        return ["scrublet", "scds", "solo", "doubletdetection"]


def run_preprocessing(
    input_path,
    output_path_plot,
    tissue,
    n_top_genes,
    save_plots,
):
    log.critical(f"Tissue: {tissue}")
    data_path = input_path / tissue
    patients = sorted(
        pd.unique(
            [
                "PT-" + p.split("-")[1]
                for p in os.listdir(data_path)
                if os.path.isdir(data_path / p)
            ]
        )
    )
    log.info(f"Patients found ({len(patients)}): {patients}")
    for PT in patients:
        for side in ["R", "L"]:
            patient_path = data_path / f"{PT}-{tissue}-{side}"

            try:
                h5_name = f"{PT}-{side}-CellBender_filtered.h5"

                log.info(f"File name: {h5_name}")

                adata = sc.read_10x_h5(patient_path / h5_name)
                adata.var_names_make_unique()

                log.info(f"Data shape:{adata.shape}")

                adata.X = adata.X.astype(np.float32)

                log.info("Read ddqc metrics")
                ddqc_obs = pd.read_csv(
                    patient_path / f"{PT}-{side}_CellBender_filtered_ddqc.csv"
                )
                ddqc_obs = ddqc_obs.set_index("barcodekey")
                ddqc_obs.index.name = None

                log.info("Basic filtering (min cells and min genes)")
                sc.pp.filter_genes(adata, min_cells=3)
                sc.pp.filter_cells(
                    adata, min_genes=200
                )  # probably already done in ddqc

                log.info("Remove MALAT1")
                malat1 = adata.var_names.str.startswith("MALAT1")
                keep = np.invert(malat1)
                adata = adata[:, keep]

                log.info("Compute pct mito, ribo, hb")

                adata.var["mt"] = adata.var_names.str.startswith("MT-")
                adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
                adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))

                sc.pp.calculate_qc_metrics(
                    adata,
                    qc_vars=["mt", "ribo", "hb"],
                    inplace=True,
                    percent_top=[20],
                    log1p=True,
                )

                log.info("Doublet detection")
                name = f"{PT}-{tissue}-{side}"
                doublets_methods = {
                    "scdblfinder": (scdblfinder, [adata]),
                    "scrublet": (scrublet, [adata]),
                    "scds": (scds, [adata]),
                    "doubletdetection": (
                        doubletdetection_method,
                        [adata, output_path_plot, name],
                    ),
                    "solo": (solo, [adata]),
                }

                doublet_methods_sequence = doublet_methods_decision(adata)

                log.critical(f"{adata.shape} shape: {doublet_methods_sequence}")
                for f in doublet_methods_sequence:
                    log.info(f"Method {f}")
                    adata = doublets_methods[f][0](*doublets_methods[f][1])
                    pred_col = adata.obs[f"predicted_doublets_{f}"]
                    if f == "doubletdetection":
                        adata.obs[f"predicted_doublets_{f}"] = pred_col.astype(float)
                    elif f in ["solo", "scdblfinder"]:
                        adata.obs[f"predicted_doublets_{f}"] = pred_col.map(
                            {"singlet": 0, "doublet": 1}
                        )

                    adata.obs[f"predicted_doublets_{f}"] = adata.obs[
                        f"predicted_doublets_{f}"
                    ].astype(int)

                adata.obs["predicted_doublets_consensus"] = (
                    adata.obs[
                        [c for c in adata.obs.columns if "predicted_doublets" in c]
                    ]
                    .all(axis=1)
                    .astype(str)
                )

                table = Table()
                table.add_column("Method")
                table.add_column("singlet")
                table.add_column("doublet")

                for method in doublet_methods_sequence + ["consensus"]:
                    try:
                        doublets_count = (
                            adata.obs[f"predicted_doublets_{method}"]
                            .value_counts()
                            .tolist()
                        )
                        table.add_row(
                            method,
                            str(doublets_count[0]),
                            str(doublets_count[1]),
                        )
                    except Exception as e:
                        log.exception(e)

                console.print(table)

                log.info(f"Highly deviant feature selection, top {n_top_genes} genes")
                dfs = ro.r(
                    """
                    g <- function(adata){
                    sce=devianceFeatureSelection(adata, assay="X")
                    return(rowData(sce)$binomial_deviance)}
                    """
                )
                binomial_deviance = ro.globalenv["g"](adata).T
                idx = binomial_deviance.argsort()[-n_top_genes:]
                mask = np.zeros(adata.var_names.shape, dtype=bool)
                mask[idx] = True

                adata.var["highly_deviant"] = mask
                adata.var["binomial_deviance"] = binomial_deviance
                adata.var["highly_variable"] = adata.var["highly_deviant"].copy()

                log.info(f"Log1p normalization")
                adata.layers["counts"] = adata.X
                scales_counts = sc.pp.normalize_total(
                    adata, target_sum=None, inplace=False
                )
                adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)

                adata.X = adata.layers["log1p_norm"]

                log.info("Plot UMAP embedding")

                for method in doublet_methods_sequence:
                    adata.obs[f"predicted_doublets_{method}"] = adata.obs[
                        f"predicted_doublets_{method}"
                    ].astype("str")

                adata.var["highly_variable"] = adata.var["highly_deviant"]

                sc.tl.pca(
                    adata, svd_solver="arpack", n_comps=20, use_highly_variable=True
                )
                sc.pp.neighbors(adata)
                sc.tl.umap(adata)

                adata.obs.index = [s.split("-")[0] for s in list(adata.obs.index)]

                log.info(
                    f"adata shape: {len(set(adata.obs.index))}, ddqc_obs shape: {len(set(ddqc_obs.index))}"
                )
                adata = adata[
                    ~adata.obs.index.isin(
                        list(set(adata.obs.index) - (set(ddqc_obs.index)))
                    )
                ].copy()  # cells do not coincide?
                adata.obs = pd.merge(
                    adata.obs, ddqc_obs, left_index=True, right_index=True
                )

                p = QC_metrics_UMAP_plot(adata)

                if save_plots:
                    output_file(
                        output_path_plot
                        / f"{PT}-{tissue}-{side}_cellbender_QC_sliders.html"
                    )

                show(p)

                adata.obs["cluster_labels"] = adata.obs["cluster_labels"].astype("str")
                adata.obs["passed_qc"] = adata.obs["passed_qc"].astype("str")
                adata.obs["side"] = side
                metrics = QC_LABELS_SAMPLE + [
                    c for c in adata.obs.columns if "doublet" in c
                ]

                sc.pl.umap(adata, color=metrics)

                tabs = [
                    TabPanel(child=interactive_embedding(adata, label), title=label)
                    for label in metrics
                ]

                p = Tabs(tabs=tabs)

                if save_plots:
                    output_file(
                        output_path_plot / f"{PT}-{tissue}-{side}_cellbender_QC.html"
                    )
                show(p)

                log.info(f"Save data")
                adata.write(data_path / f"{PT}-{tissue}-{side}_QC.h5ad")

            except Exception as e:
                log.exception(e)
