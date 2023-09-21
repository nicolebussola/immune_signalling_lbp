import datetime
import logging
import os
from pathlib import Path

import anndata2ri
import numpy as np
import pandas as pd
import rpy2.robjects as ro
import scanpy as sc
from bokeh.models import TabPanel, Tabs
from bokeh.plotting import output_file, show
from rpy2.robjects.packages import importr
from scipy.stats import median_abs_deviation

from ..utils import QC_metrics_UMAP_plot, interactive_embedding

logging.getLogger().setLevel(logging.INFO)
logging.basicConfig(format="%(asctime)s - %(message)s", datefmt="%d-%b-%y %H:%M:%S")
anndata2ri.activate()
timestamp = datetime.datetime.now().strftime("%m%d")
LABELS = [
    "log1p_total_counts",
    "log1p_n_genes_by_counts",
    "pct_counts_in_top_20_genes",
    "outlier",
    "pct_counts_mt",
    "pct_counts_ribo",
    "scDblFinder_score",
    "doublet_scores_scrublet",
    "scDblFinder_class",
    "predicted_doublets_scrublet",
    "cluster_labels",
    "passed_qc",
]

sc.settings.verbosity = 0

importr("Seurat")
importr("scater")
importr("scDblFinder")
importr("BiocParallel")
importr("scry")


def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


np.random.seed(42)

## Batch 1 ==> narsad_cellRanger_outs
## Batch 2


def run_preprocessing(input_path, output_path_plot, tissue, n_top_genes):
    data_path = input_path / tissue
    patients = pd.unique(
        [
            "PT-" + p.split("-")[1]
            for p in os.listdir(data_path)
            if os.path.isdir(data_path / p)
        ]
    )
    logging.info(f"Patients found ({len(patients)}): {patients}")
    for PT in patients:
        for side in ["R", "L"]:
            patient_path = data_path / f"{PT}-{tissue}-{side}"

            try:
                h5_name = f"{PT}-{side}_CellBender_filtered.h5"

                logging.info("\n")
                logging.info(f"File name:{h5_name}")

                adata = sc.read_10x_h5(patient_path / h5_name)

                logging.info(f"Data shape:{adata.shape}")

                adata.X = adata.X.astype(np.float32)

                adata.var_names_make_unique()

                logging.info("Read ddqc output")
                ddqc_obs = pd.read_csv(
                    patient_path / f"{PT}-{side}_CellBender_filtered_ddqc.csv"
                )
                ddqc_obs = ddqc_obs.set_index("barcodekey")
                ddqc_obs.index.name = None

                logging.info("Basic filtering")
                sc.pp.filter_genes(adata, min_cells=3)
                sc.pp.filter_cells(adata, min_genes=200)

                # Remove MALAT1
                malat1 = adata.var_names.str.startswith("MALAT1")
                keep = np.invert(malat1)
                adata = adata[:, keep]

                logging.info("Compute low-quality read")

                # mitochondrial genes
                adata.var["mt"] = adata.var_names.str.startswith("MT-")
                # ribosomal genes
                adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
                # hemoglobin genes.
                adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))

                sc.pp.calculate_qc_metrics(
                    adata,
                    qc_vars=["mt", "ribo", "hb"],
                    inplace=True,
                    percent_top=[20],
                    log1p=True,
                )

                logging.info("Compute MADs outliers")
                adata.obs["outlier"] = (
                    is_outlier(adata, "log1p_total_counts", 2)
                    | is_outlier(adata, "log1p_n_genes_by_counts", 2)
                    | is_outlier(adata, "pct_counts_in_top_20_genes", 2)
                )

                logging.info(f"Outliers (MADs):{adata.obs.outlier.value_counts()}")

                logging.info("Doublet detection")
                data_mat = adata.X.T.copy()

                scdblfinder = ro.r(
                    """
                    f <- function(data_mat){set.seed(123)
                    sce = scDblFinder(
                        SingleCellExperiment(
                            list(counts=data_mat),
                        )
                    )
                    doublet_score = sce$scDblFinder.score
                    doublet_class = sce$scDblFinder.class
                    return(list(doublet_score, doublet_class))}
                        """
                )

                doublet_results = scdblfinder(data_mat)
                adata.obs["scDblFinder_score"] = doublet_results[0]
                adata.obs["scDblFinder_class"] = doublet_results[1]
                adata.obs["scDblFinder_class"] = adata.obs["scDblFinder_class"].astype(
                    "object"
                )

                print(
                    f"Doublets (scDblFinder): {adata.obs.scDblFinder_class.value_counts()}"
                )

                sc.external.pp.scrublet(adata, random_state=123)
                adata.obs.rename(
                    columns={
                        "doublet_score": "doublet_scores_scrublet",
                        "predicted_doublet": "predicted_doublets_scrublet",
                    },
                    inplace=True,
                )

                print(
                    f"Doublets (scrublet): {adata.obs.predicted_doublets_scrublet.value_counts()}"
                )

                logging.info(
                    f"Highly deviant feature selection, top {n_top_genes} genes"
                )
                dfs = ro.r(
                    """
                    f <- function(adata){
                    sce=devianceFeatureSelection(adata, assay="X")
                    return(rowData(sce)$binomial_deviance)}
                    """
                )
                binomial_deviance = dfs(adata).T
                idx = binomial_deviance.argsort()[-n_top_genes:]
                mask = np.zeros(adata.var_names.shape, dtype=bool)
                mask[idx] = True

                adata.var["highly_deviant"] = mask
                adata.var["binomial_deviance"] = binomial_deviance
                adata.var["highly_variable"] = adata.var["highly_deviant"]

                logging.info(f"Log1p normalization")
                adata.layers["counts"] = adata.X
                scales_counts = sc.pp.normalize_total(
                    adata, target_sum=None, inplace=False
                )
                # log1p transform
                adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)

                adata.X = adata.layers["log1p_norm"]
                logging.info(f"Plot UMAP embedding")
                adata.obs["scDblFinder_class"] = adata.obs["scDblFinder_class"].astype(
                    "str"
                )
                adata.obs["predicted_doublets_scrublet"] = adata.obs[
                    "predicted_doublets_scrublet"
                ].astype("str")
                adata.var["highly_variable"] = adata.var["highly_deviant"]

                sc.tl.pca(
                    adata, svd_solver="arpack", n_comps=20, use_highly_variable=True
                )
                # sc.tl.paga(adata_deviance)
                # sc.pl.paga(adata_deviance, plot=False)
                sc.pp.neighbors(adata)
                # sc.tl.umap(adata_deviance, init_pos='paga')
                sc.tl.umap(adata)

                adata.obs.index = [s.split("-")[0] for s in list(adata.obs.index)]
                adata = adata[
                    ~adata.obs.index.isin(
                        list(set(adata.obs.index) - (set(ddqc_obs.index)))
                    )
                ].copy()  # cells do not coincide?
                adata.obs = pd.merge(
                    adata.obs, ddqc_obs, left_index=True, right_index=True
                )

                p = QC_metrics_UMAP_plot(adata)

                output_file(
                    output_path_plot
                    / f"{PT}-{tissue}-{side}_cellbender_QC_sliders.html"
                )

                show(p)

                adata.obs["outlier"] = adata.obs["outlier"].astype("object")
                adata.obs["cluster_labels"] = adata.obs["cluster_labels"].astype("str")
                adata.obs["passed_qc"] = adata.obs["passed_qc"].astype("str")

                sc.pl.umap(adata, color=LABELS)

                tabs = [
                    TabPanel(child=interactive_embedding(adata, label), title=label)
                    for label in LABELS
                ]

                p = Tabs(tabs=tabs)
                output_file(
                    output_path_plot / f"{PT}-{tissue}-{side}_cellbender_QC.html"
                )
                show(p)

                adata.obs["side"] = side
                adata.obs["outlier"] = adata.obs["outlier"].astype("str")

                adata.write(data_path / f"{PT}-{tissue}-{side}_QC.h5ad")
            except Exception as e:
                print(e)
