import logging
from pathlib import Path

import anndata2ri
import numpy as np
import pandas as pd
import rpy2.robjects as ro
import scanpy as sc
import scipy.sparse
import scrublet as scr
import seaborn as sns
from bokeh.models import ColorBar, TabPanel, Tabs
from bokeh.plotting import output_file, output_notebook, show

# import rpy2.rinterfaceib.callbacks as rcb
from rpy2.robjects.packages import importr
from scipy.stats import median_abs_deviation

from utils import QC_metrics_UMAP_plot, interactive_embedding

sc.settings.verbosity = 0

anndata2ri.activate()

np.random.seed(42)

BATCH = 1
tissue = "blood"
# PT = "PT-212"
# side = "R"

if BATCH == 1:
    batch_path = Path("../LBP_brain_blood_pairs/data/narsad_cellRanger_outs/")
else:
    batch_path = Path()


for PT in [
    "PT-182",
    "PT-185",
    "PT-201",
    "PT-203",
    "PT-205",
    "PT-206",
    "PT-208",
    "PT-212",
    "PT-214",
]:
    for side in ["R", "L"]:
        try:
            if tissue == "blood":
                h5_name = f"{PT}-{side}-B_CellBender_filtered.h5"
            else:
                h5_name = f"{PT}-{side}_CellBender_filtered.h5"

            print("\n======================")
            print("File name:\n")
            print(h5_name)
            print("======================\n")
            adata = sc.read_10x_h5(
                batch_path / tissue / f"{PT}-{tissue}-{side}" / h5_name
            )
            print("\n======================")
            print("Data shape:\n")
            print(adata.shape)
            print("======================\n")

            adata.X = adata.X.astype(np.float32)

            adata.var_names_make_unique()

            ############## Basic Filtering
            sc.pp.filter_genes(adata, min_cells=3)
            sc.pp.filter_cells(adata, min_genes=200)
            # Remove MALAT1
            malat1 = adata.var_names.str.startswith("MALAT1")
            keep = np.invert(malat1)
            adata = adata[:, keep]

            ############ Low-quality reads

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

            ############ MADs outliers
            def is_outlier(adata, metric: str, nmads: int):
                M = adata.obs[metric]
                outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
                    np.median(M) + nmads * median_abs_deviation(M) < M
                )
                return outlier

            adata.obs["outlier"] = (
                is_outlier(adata, "log1p_total_counts", 2)
                | is_outlier(adata, "log1p_n_genes_by_counts", 2)
                | is_outlier(adata, "pct_counts_in_top_20_genes", 2)
            )
            print("\n======================")
            print("Outliers (MADs):\n")
            print(adata.obs.outlier.value_counts())
            print("======================\n")

            ############ Doublet detection
            importr("Seurat")

            importr("scater")
            importr("scDblFinder")
            importr("BiocParallel")

            data_mat = adata.X.T.copy()
            ro.globalenv["data_mat"] = data_mat

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

            print("\n======================")
            print("Doublets (scDblFinder):\n")
            print(adata.obs.scDblFinder_class.value_counts())
            print("======================\n")

            sc.external.pp.scrublet(adata, random_state=123)
            adata.obs.rename(
                columns={
                    "doublet_score": "doublet_scores_scrublet",
                    "predicted_doublet": "predicted_doublets_scrublet",
                },
                inplace=True,
            )

            print("\n======================")
            print("Doublets (scrublet):\n")
            print(adata.obs["predicted_doublets_scrublet"].value_counts())
            print("======================\n")

            ############ Feature selection
            ro.globalenv["adata"] = adata
            importr("scry")
            ro.r(
                """
                sce = devianceFeatureSelection(adata, assay="X")
                    """
            )
            binomial_deviance = ro.r("rowData(sce)$binomial_deviance").T
            idx = binomial_deviance.argsort()[-4000:]
            mask = np.zeros(adata.var_names.shape, dtype=bool)
            mask[idx] = True

            adata.var["highly_deviant"] = mask
            adata.var["binomial_deviance"] = binomial_deviance
            adata.var["highly_variable"] = adata.var["highly_deviant"]

            ############ Log1p normalization
            adata.layers["counts"] = adata.X
            scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
            # log1p transform
            adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)

            adata.X = adata.layers["log1p_norm"]

            ############## Plots
            adata.obs["scDblFinder_class"] = adata.obs["scDblFinder_class"].astype(
                "str"
            )
            adata.obs["predicted_doublets_scrublet"] = adata.obs[
                "predicted_doublets_scrublet"
            ].astype("str")
            adata.var["highly_variable"] = adata.var["highly_deviant"]

            sc.tl.pca(adata, svd_solver="arpack", n_comps=20, use_highly_variable=True)
            # sc.tl.paga(adata_deviance)
            # sc.pl.paga(adata_deviance, plot=False)
            sc.pp.neighbors(adata)
            # sc.tl.umap(adata_deviance, init_pos='paga')
            sc.tl.umap(adata)

            p = QC_metrics_UMAP_plot(adata)

            output_file(
                f"../LBP_brain_blood_pairs/PLOTS/QC_plots/{PT}-{tissue}-{side}_cellbender_QC_sliders.html"
            )

            show(p)
            adata.obs["outlier"] = adata.obs["outlier"].astype("object")

            labels = [
                "log1p_total_counts",
                "outlier",
                "pct_counts_mt",
                "pct_counts_ribo",
                "pct_counts_hb",
                "scDblFinder_score",
                "doublet_scores_scrublet",
                "scDblFinder_class",
                "predicted_doublets_scrublet",
            ]
            sc.pl.umap(adata, color=labels)

            tabs = []

            for label in labels:
                tabs.append(
                    TabPanel(child=interactive_embedding(adata, label), title=label)
                )

            p = Tabs(tabs=tabs)
            output_file(
                f"../LBP_brain_blood_pairs/PLOTS/QC_plots/{PT}-{tissue}-{side}_cellbender_QC.html"
            )
            show(p)

            # adata_filtered = adata[(adata.obs['predicted_doublets_scrublet'] == 'False') & (v.obs['pct_counts_mt']<20)].copy()

            adata.obs["side"] = side
            adata.obs["outlier"] = adata.obs["outlier"].astype("str")

            adata.write(batch_path / tissue / f"{PT}-{tissue}-{side}_QC.h5ad")
        except Exception as e:
            print(e)
