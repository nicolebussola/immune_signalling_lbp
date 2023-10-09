import datetime
import logging
import os
import warnings
from pathlib import Path

import anndata as ad
import anndata2ri
import numpy as np
import pandas as pd
import rpy2.robjects as ro
import scanpy as sc
import scanpy.external as sce
import scib
import scvi
from rich.logging import RichHandler
from rpy2.robjects.packages import importr
from scipy.sparse import csr_matrix, issparse
from tabulate import tabulate
from tqdm import tqdm

from ..labels import (
    CHEMISTRY_V2_PATIENTS,
    FAILED_QC_SAMPLES,
    QC_LABELS_BATCH_1,
    QC_LABELS_BATCH_2,
)
from ..utils import gridlayout

importr("scry")
importr("scran")
importr("BiocParallel")
warnings.filterwarnings("ignore")
anndata2ri.activate()
np.random.seed(42)

logging.getLogger().setLevel(logging.INFO)
FORMAT = "%(message)s"
logging.basicConfig(
    level="INFO", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
)

log = logging.getLogger("rich")
timestamp = datetime.datetime.now().strftime("%m%d")


def run_merge(
    input_path,
    batch,
    output_path,
    tissue,
    output_path_plot,
    method_hvg,
    n_top_genes,
):
    sample_path = input_path / tissue

    samples_adata = []
    sum_cells = 0
    all_genes = []

    log.critical("Read and merge QCed data")
    sample_list = [
        s for s in os.listdir(sample_path) if os.path.splitext(s)[1] == ".h5ad"
    ]
    for sample in tqdm(sorted(sample_list)):
        adata_sample = sc.read(sample_path / sample)
        adata_sample.obs["pt"] = "-".join(sample.split("-")[:2])
        all_genes.extend(list(adata_sample.var.index))
        if sample in FAILED_QC_SAMPLES:
            adata_sample = adata_sample[
                adata_sample.obs[adata_sample.obs["log1p_total_counts"] > 7.73].index, :
            ]
        samples_adata.append(adata_sample)
        sum_cells += adata_sample.shape[0]

    log.info(f"Samples: {sorted(sample_list)}")
    log.info(f"Total genes: {len(all_genes)}")
    log.info(f"Sum cells: {sum_cells}")

    adata_merged = ad.concat(samples_adata, join="outer", index_unique="_")
    del samples_adata
    adata_merged.obs_names_make_unique()
    adata_merged.var_names_make_unique()
    adata_merged.obs["chemistry"] = adata_merged.obs["pt"].apply(
        lambda x: "v2" if x in CHEMISTRY_V2_PATIENTS else "v3"
    )
    adata_merged.X = adata_merged.layers["counts"].copy()

    log.critical("Data Filtering")

    adata_filtered = adata_merged[
        (adata_merged.obs["doublet_scores_scrublet"] < 0.2)
        & (adata_merged.obs["passed_qc"] == "True")
        & (adata_merged.obs["pct_counts_mt"] < 40)
    ].copy()

    log.info(f"Merged data prefiltering: {adata_merged.shape}")
    log.info(f"Merged data after filtering: {adata_filtered.shape}")

    del adata_merged

    adata_filtered.X = adata_filtered.layers["counts"].copy()

    log.critical("Normalize data")
    log.info("Log1p normalization")
    scales_counts = sc.pp.normalize_total(
        adata_filtered,
        layer="counts",
        exclude_highly_expressed=True,
        target_sum=None,
        inplace=False,
    )
    adata_filtered.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)

    log.info("Scran normalization")

    adata_pp = adata_filtered.copy()

    logging.info("Preclustering data....")
    sc.pp.normalize_total(adata_pp)
    sc.pp.log1p(adata_pp)
    sc.pp.pca(adata_pp, n_comps=25)
    sc.pp.neighbors(adata_pp)
    sc.tl.leiden(adata_pp, key_added="groups")

    log.info("Computing size factors........")
    data_mat = adata_pp.X.T
    # convert to CSC if possible. See https://github.com/MarioniLab/scran/issues/70
    if issparse(data_mat):
        if data_mat.nnz > 2**31 - 1:
            data_mat = data_mat.tocoo()

        else:
            data_mat = data_mat.tocsc()

    input_groups = adata_pp.obs["groups"]
    del adata_pp

    size_factors = ro.r(
        """
        f <- function(data_mat, input_groups){sizeFactors(
        computeSumFactors(
            SingleCellExperiment(
                list(counts=data_mat)),
                clusters = input_groups,
                min.mean = 0.1,
                BPPARAM = MulticoreParam()
        )
    )}
            """
    )
    adata_filtered.obs["size_factors"] = size_factors(data_mat, input_groups)
    scran = adata_filtered.X / adata_filtered.obs["size_factors"].values[:, None]
    adata_filtered.layers["scran_normalization"] = csr_matrix(sc.pp.log1p(scran))

    log.info("Pearson residuals")
    adata_filtered.X = adata_filtered.layers["counts"].copy()

    analytic_pearson = sc.experimental.pp.normalize_pearson_residuals(
        adata_filtered, inplace=False
    )
    adata_filtered.layers["analytic_pearson_residuals"] = csr_matrix(
        analytic_pearson["X"]
    )

    log.info("Write filtered and normalized adata")
    adata_filtered.write(output_path / f"{tissue}_filtered_normalized.h5ad")

    log.critical("Feature selection")

    log.info(f"Compute n={n_top_genes} top genes with {method_hvg} approach")
    if method_hvg == "cell_ranger":
        sc.pp.highly_variable_genes(
            adata_filtered,
            n_top_genes=n_top_genes,
            flavor="cell_ranger",
            batch_key="pt",
        )
    elif method_hvg == "HighlyDeviant":
        dfs = ro.r(
            """
            f <- function(adata){
            sce=devianceFeatureSelection(adata, assay="X")
            return(rowData(sce)$binomial_deviance)}
            """
        )
        binomial_deviance = dfs(adata_filtered).T
        idx = binomial_deviance.argsort()[-n_top_genes:]
        mask = np.zeros(adata_filtered.var_names.shape, dtype=bool)
        mask[idx] = True
        adata_filtered.var["highly_deviant"] = mask
        adata_filtered.var["binomial_deviance"] = binomial_deviance
        adata_filtered.var["highly_variable"] = adata_filtered.var["highly_deviant"]

    adata_hvg = adata_filtered[:, adata_filtered.var["highly_variable"]].copy()
    adata_hvg.X = adata_hvg.layers["counts"].copy()
    sc.pp.normalize_total(adata_hvg)
    sc.pp.log1p(adata_hvg)
    sc.pp.pca(adata_hvg, svd_solver="arpack")
    adata_hvg.layers["logcounts"] = adata_hvg.X.copy()

    LABELS = QC_LABELS_BATCH_1 if batch == "1" else QC_LABELS_BATCH_2

    if output_path_plot != Path():
        logging.info("Plot data")
        sc.pp.neighbors(adata_hvg)
        sc.tl.umap(adata_hvg)
        sc.pl.umap(adata_hvg, color=LABELS)
        gridlayout(
            LABELS,
            adata_hvg,
            fname=output_path_plot
            / f"merged_{tissue}_cellbender_QC_filtered_{n_top_genes}{method_hvg}_{timestamp}.html",
        )

    log.info(f"Compute Unintegration metrics per batch key")
    batch_keys = ["chemistry", "pt", "side"] if batch == "1" else ["pt", "side"]

    unintegrated_metrics = []
    for bk in batch_keys:
        log.info(f"Compute metrics per batch {bk}")

        unintegrated_metrics.extend(
            (
                scib.metrics.pcr_comparison(adata_filtered, adata_hvg, covariate=bk),
                scib.me.cell_cycle(
                    adata_filtered, adata_hvg, batch_key=bk, organism="human"
                ),
                scib.me.hvg_overlap(adata_filtered, adata_hvg, batch_key=bk),
            )
        )

    del adata_filtered
    log.critical("Data integration")
    log.info("Harmony integration")
    adata_harmony = adata_hvg.copy()

    if batch == "1":
        log.info("Integrate chemistry...\n")
        sce.pp.harmony_integrate(
            adata_harmony,
            "chemistry",
            adjusted_basis="X_pca_harmony_chem",
        )
        sc.pp.neighbors(adata_harmony, use_rep="X_pca_harmony_chem")
        sc.tl.umap(adata_harmony)

    basis = "X_pca_harmony_chem" if batch == "1" else "X_pca"
    log.info("Integrate pt...\n")
    sce.pp.harmony_integrate(
        adata_harmony,
        "pt",
        basis=basis,
        adjusted_basis="X_pca_harmony_pt",
    )
    sc.pp.neighbors(adata_harmony, use_rep="X_pca_harmony_pt")
    sc.tl.umap(adata_harmony)

    if output_path_plot != Path():
        log.info("Plot data")
        sc.pl.umap(adata_harmony, color=LABELS, wspace=1)
        gridlayout(
            LABELS,
            adata_harmony,
            fname=output_path_plot
            / f"merged_{tissue}_cellbender_QC_filtered_{n_top_genes}{method_hvg}_harmonyPt_{timestamp}.html",
        )

    log.info("ScVI integration")
    adata_scvi = adata_hvg.copy()

    if batch == "1":
        scvi.model.SCVI.setup_anndata(
            adata_scvi,
            layer="counts",
            batch_key="chemistry",
            categorical_covariate_keys=["pt"],
        )
    else:
        scvi.model.SCVI.setup_anndata(
            adata_scvi,
            layer="counts",
            batch_key="pt",
        )
    model_scvi = scvi.model.SCVI(adata_scvi)
    model_scvi.train(use_gpu=False, max_epochs=None)
    adata_scvi.obsm["X_scVI"] = model_scvi.get_latent_representation()
    adata_scvi.layers["X_normalized_scVI"] = model_scvi.get_normalized_expression()
    sc.pp.neighbors(adata_scvi, use_rep="X_scVI")
    sc.tl.umap(adata_scvi)

    if output_path_plot != Path():
        log.info("Plot data")
        sc.pp.neighbors(adata_scvi, use_rep="X_scVI")
        sc.tl.umap(adata_scvi)
        sc.pl.umap(adata_scvi, color=LABELS)
        gridlayout(
            LABELS,
            adata_scvi,
            fname=output_path_plot
            / f"merged_{tissue}_cellbender_QC_filtered_{n_top_genes}{method_hvg}_scvi_{timestamp}.html",
        )

    logging.critical("Leiden clustering and saving adata")
    for adata in [adata_harmony, adata_scvi]:
        for k, res in zip(
            [
                "leiden_res0_5",
                "leiden_res1",
                "leiden_res1_25",
                "leiden_res1_5",
                "leiden_res1_75",
            ],
            [0.5, 1, 1.25, 1.5, 1.75],
        ):
            sc.tl.leiden(adata, key_added=k, resolution=res)

    adata_harmony.write(
        output_path / f"{tissue}_{n_top_genes}{method_hvg}_harmony_clustered.h5ad"
    )
    adata_scvi.write(
        output_path / f"{tissue}_{n_top_genes}{method_hvg}_scvi_clustered.h5ad"
    )

    log.critical("Compare integration methods")

    harmony_metrics = []
    scvi_metrics = []

    for bk in batch_keys:
        log.info(f"Compute metrics per batch key = {bk}")

        harmony_metrics.extend(
            (
                scib.metrics.pcr_comparison(adata_hvg, adata_harmony, covariate=bk),
                scib.me.cell_cycle(
                    adata_hvg, adata_harmony, batch_key=bk, organism="human"
                ),
                scib.me.hvg_overlap(adata_hvg, adata_harmony, batch_key=bk),
            )
        )

        scvi_metrics.extend(
            (
                scib.metrics.pcr_comparison(adata_hvg, adata_scvi, covariate=bk),
                scib.me.cell_cycle(
                    adata_hvg, adata_scvi, batch_key=bk, organism="human"
                ),
                scib.me.hvg_overlap(adata_hvg, adata_scvi, batch_key=bk),
            )
        )

    batch_tables = []
    for i in range(len(batch_keys)):
        table = [
            [f"Unintegrated_{batch_keys[i]}"] + unintegrated_metrics[3 * i : 3 * i + 3],
            [f"harmony_{batch_keys[i]}"] + harmony_metrics[3 * i : 3 * i + 3],
            [f"scvi_{batch_keys[i]}"] + scvi_metrics[:3],
        ]
        batch_tables.append(table)
        print(
            tabulate(
                table,
                headers=[batch_keys[i], "PCR comparison", "Cell cycle", "HVG"],
                tablefmt="fancy_outline",
            )
        )

    df = pd.concat([pd.DataFrame(t) for t in batch_tables])
    df.columns = ["method", "PCR comparison", "Cell cycle", "HVG"]

    df.style.background_gradient(cmap="Blues")
    df.to_csv(
        input_path
        / f"{tissue}_{n_top_genes}{method_hvg}_integration_comparison_{timestamp}.csv"
    )
