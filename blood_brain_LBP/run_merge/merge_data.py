import datetime
import logging
import os
import warnings

import anndata as ad
import anndata2ri
import numpy as np
import pandas as pd
import rpy2.robjects as ro
import scanpy as sc
import scanpy.external as sce
import scib
import scvi
from rpy2.robjects.packages import importr
from scipy.sparse import csr_matrix, issparse
from tqdm import tqdm

from ..utils import gridlayout

warnings.filterwarnings("ignore")
anndata2ri.activate()
importr("scry")
importr("scran")
importr("BiocParallel")
np.random.seed(42)

logging.getLogger().setLevel(logging.INFO)
logging.basicConfig(format="%(asctime)s - %(message)s", datefmt="%d-%b-%y %H:%M:%S")

timestamp = datetime.datetime.now().strftime("%m%d")
LABELS = [
    "chemistry",
    "pt",
    "side",
    "total_counts",
    "n_genes_by_counts",
    "pct_counts_mt",
    "pct_counts_ribo",
    "scDblFinder_score",
    "scDblFinder_class",
    "doublet_scores_scrublet",
    "predicted_doublets_scrublet",
    "scrublet_class_refined",
    "passed_qc",
]


def run_merge(
    input_path,
    output_path,
    tissue,
    save_plots,
    output_path_plot,
    method_hvg,
    n_top_genes,
):
    sample_path = input_path / tissue

    samples_adata = []
    sum_cells = 0
    all_genes = []

    logging.info("########## Read and merge data ##########")
    sample_list = [
        s for s in os.listdir(sample_path) if os.path.splitext(s)[1] == ".h5ad"
    ]
    # Read anndata of single sample after QC
    for sample in tqdm(sorted(sample_list)):
        adata_sample = sc.read(sample_path / sample)
        adata_sample.obs["pt"] = "-".join(sample.split("-")[:2])
        all_genes.extend(list(adata_sample.var.index))
        samples_adata.append(adata_sample)
        sum_cells += adata_sample.shape[0]

    logging.info(f"Samples: {sorted(sample_list)}")
    logging.info(f"Sum cells: {sum_cells}")
    # Merge data
    adata_merged = ad.concat(samples_adata, join="outer", index_unique="_")
    adata_merged.obs_names_make_unique()
    adata_merged.X = adata_merged.layers["counts"].copy()

    # Add information of different sequencing techniques
    adata_merged.obs["chemistry"] = adata_merged.obs["pt"].apply(
        lambda x: "v2" if x in ["PT-182", "PT-185"] else "v3"
    )

    adata_merged.obs["scrublet_class_refined"] = adata_merged.obs[
        "doublet_scores_scrublet"
    ].apply(lambda x: "singlet" if x < 0.2 else "doublet")

    logging.info("########## log1p transform ##########")

    scales_counts = sc.pp.normalize_total(
        adata_merged,
        layer="counts",
        exclude_highly_expressed=True,
        target_sum=None,
        inplace=False,
    )
    adata_merged.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)

    adata_log1p = adata_merged.copy()
    adata_log1p.X = adata_log1p.layers["log1p_norm"]
    sc.pp.pca(adata_log1p, svd_solver="arpack", n_comps=25)
    sc.pp.neighbors(adata_log1p)
    sc.tl.umap(adata_log1p)

    if save_plots:
        logging.info("########## Plot data ##########")
        sc.pl.umap(adata_log1p, color=LABELS, show=False, return_fig=False)
        gridlayout(
            LABELS,
            adata_log1p,
            fname=output_path_plot / f"merged_{tissue}_cellbender_QC_{timestamp}.html",
        )

    logging.info("########## Data Filtering ##########")

    adata_filtered = adata_merged[
        (adata_merged.obs["scrublet_class_refined"] == "singlet")
        & (adata_merged.obs["passed_qc"] == "True")
        & (adata_merged.obs["pct_counts_mt"] < 40)
    ].copy()

    logging.info(f"Merged data prefiltering: {adata_merged.shape}")
    logging.info(f"Merged data after filtering: {adata_filtered.shape}")

    adata_filtered.X = adata_filtered.layers["counts"].copy()

    logging.info("########## Normalize and log1p transform ##########")
    logging.info("Log1p normalization")
    scales_counts = sc.pp.normalize_total(
        adata_filtered,
        layer="counts",
        exclude_highly_expressed=True,
        target_sum=None,
        inplace=False,
    )

    adata_filtered.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)
    adata_filtered.X = adata_filtered.layers["log1p_norm"].copy()

    filtered_labels = LABELS[:-6]
    sc.pp.pca(adata_filtered, svd_solver="arpack", n_comps=25)
    sc.pp.neighbors(adata_filtered)
    sc.tl.umap(adata_filtered)

    if save_plots:
        logging.info("########## Plot data ##########")
        sc.pl.umap(adata_filtered, color=LABELS)
        gridlayout(
            filtered_labels,
            adata_filtered,
            fname=output_path_plot
            / f"merged_{tissue}_cellbender_QC_filtered_{timestamp}.html",
        )

    logging.info("Scran normalization")

    adata_filtered.X = adata_filtered.layers["counts"].copy()
    adata_pp = adata_filtered.copy()
    sc.pp.normalize_total(adata_pp)  # already normalized
    sc.pp.log1p(adata_pp)
    sc.pp.pca(adata_pp, n_comps=25)
    sc.pp.neighbors(adata_pp)
    sc.tl.leiden(adata_pp, key_added="groups")

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

    logging.info("Pearson residuals")
    adata_filtered.X = adata_filtered.layers["counts"].copy()

    analytic_pearson = sc.experimental.pp.normalize_pearson_residuals(
        adata_filtered, inplace=False
    )
    adata_filtered.layers["analytic_pearson_residuals"] = csr_matrix(
        analytic_pearson["X"]
    )

    logging.info("########## Feature selection ##########")

    logging.info(f"Compute n={n_top_genes} top genes with {method_hvg} approach")
    # adata_filtered.write(
    #     output_path
    #     / tissue
    #     / f"{tissue}_filtered.h5ad"
    # )
    if method_hvg == "cell_ranger":
        sc.pp.highly_variable_genes(
            adata_filtered,
            n_top_genes=n_top_genes,
            flavor="cell_ranger",
            batch_key="chemistry",
        )
    elif method_hvg == "HighlyDeviant":
        deviance_feat_sel = ro.r(
            """
        f <- function(adata_filtered){
        sce = devianceFeatureSelection(adata_filtered, assay="X")
        return(rowData(sce)$binomial_deviance)
        }
        """
        )
        binomial_deviance = deviance_feat_sel(adata_filtered).T
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

    sc.pp.neighbors(adata_hvg)
    sc.tl.umap(adata_hvg)
    if save_plots:
        logging.info("########## Plot data ##########")
        sc.pl.umap(adata_hvg, color=filtered_labels)
        gridlayout(
            filtered_labels,
            adata_hvg,
            fname=output_path_plot
            / f"merged_{tissue}_cellbender_QC_filtered_{n_top_genes}{method_hvg}_{timestamp}.html",
        )

    logging.info("########## Data integration ##########")
    logging.info("Harmony integration")
    adata_harmony = adata_hvg.copy()

    logging.info("Integrate technology...\n")
    # Correct per technology
    sce.pp.harmony_integrate(
        adata_harmony, key="chemistry", adjusted_basis="X_pca_harmony_chem"
    )
    sc.pp.neighbors(adata_harmony, use_rep="X_pca_harmony_chem")
    sc.tl.umap(adata_harmony)

    logging.info("Integrate pt...\n")
    # Correct per technology and pt
    sce.pp.harmony_integrate(
        adata_harmony,
        "pt",
        basis="X_pca_harmony_chem",
        adjusted_basis="X_pca_harmony_chempt",
    )
    sc.pp.neighbors(adata_harmony, use_rep="X_pca_harmony_chempt")
    sc.tl.umap(adata_harmony)

    logging.info("Integrate side...\n")
    # Correct per technology, pt, and and side
    sce.pp.harmony_integrate(
        adata_harmony,
        "side",
        basis="X_pca_harmony_chempt",
        adjusted_basis="X_pca_harmony_chemptside",
    )
    sc.pp.neighbors(adata_harmony, use_rep="X_pca_harmony_chemptside")
    sc.tl.umap(adata_harmony)
    if save_plots:
        logging.info("########## Plot data ##########")
        sc.pl.umap(adata_harmony, color=filtered_labels, wspace=1)
        gridlayout(
            filtered_labels,
            adata_harmony,
            fname=output_path_plot
            / f"merged_{tissue}_cellbender_QC_filtered_{n_top_genes}{method_hvg}_harmonyChemPtSide_{timestamp}.html",
        )

    logging.info("ScVI integration")
    adata_scvi = adata_hvg.copy()

    scvi.model.SCVI.setup_anndata(
        adata_scvi,
        layer="counts",
        batch_key="pt",
        categorical_covariate_keys=["chemistry", "side"],
    )
    model_scvi = scvi.model.SCVI(adata_scvi)
    model_scvi.train(use_gpu=False)
    adata_scvi.obsm["X_scVI"] = model_scvi.get_latent_representation()
    adata_scvi.layers["X_normalized_scVI"] = model_scvi.get_normalized_expression()

    sc.pp.neighbors(adata_scvi, use_rep="X_scVI")
    sc.tl.umap(adata_scvi)
    if save_plots:
        logging.info("########## Plot data##########")
        sc.pl.umap(adata_scvi, color=filtered_labels)
        gridlayout(
            filtered_labels,
            adata_scvi,
            fname=output_path_plot
            / f"merged_{tissue}_cellbender_QC_filtered_{n_top_genes}{method_hvg}_scvi_{timestamp}.html",
        )

    logging.info("########## Compare integration methods##########")
    unintegrated_metrics_pcr = [
        scib.metrics.pcr_comparison(adata_filtered, adata_hvg, covariate="chemistry"),
        scib.metrics.pcr_comparison(adata_filtered, adata_hvg, covariate="pt"),
        scib.metrics.pcr_comparison(adata_filtered, adata_hvg, covariate="side"),
    ]
    harmony_metrics_pcr = [
        scib.metrics.pcr_comparison(
            adata_hvg,
            adata_harmony,
            covariate="chemistry",
            embed="X_pca_harmony_chemptside",
        ),
        scib.metrics.pcr_comparison(
            adata_hvg, adata_harmony, covariate="pt", embed="X_pca_harmony_chemptside"
        ),
        scib.metrics.pcr_comparison(
            adata_hvg, adata_harmony, covariate="side", embed="X_pca_harmony_chemptside"
        ),
    ]
    scvi_metrics_pcr = [
        scib.metrics.pcr_comparison(
            adata_hvg, adata_scvi, covariate="chemistry", embed="X_scVI"
        ),
        scib.metrics.pcr_comparison(
            adata_hvg, adata_scvi, covariate="pt", embed="X_scVI"
        ),
        scib.metrics.pcr_comparison(
            adata_hvg, adata_scvi, covariate="side", embed="X_scVI"
        ),
    ]

    unintegrated_metrics_cellc = [
        scib.me.cell_cycle(
            adata_filtered, adata_hvg, batch_key="chemistry", organism="human"
        ),
        scib.me.cell_cycle(adata_filtered, adata_hvg, batch_key="pt", organism="human"),
        scib.me.cell_cycle(
            adata_filtered, adata_hvg, batch_key="side", organism="human"
        ),
    ]
    harmony_metrics_cellc = [
        scib.me.cell_cycle(
            adata_hvg,
            adata_harmony,
            batch_key="chemistry",
            organism="human",
            embed="X_pca_harmony_chemptside",
        ),
        scib.me.cell_cycle(
            adata_hvg,
            adata_harmony,
            batch_key="pt",
            organism="human",
            embed="X_pca_harmony_chemptside",
        ),
        scib.me.cell_cycle(
            adata_hvg,
            adata_harmony,
            batch_key="side",
            organism="human",
            embed="X_pca_harmony_chemptside",
        ),
    ]
    scvi_metrics_cellc = [
        scib.me.cell_cycle(
            adata_hvg,
            adata_scvi,
            batch_key="chemistry",
            organism="human",
            embed="X_scVI",
        ),
        scib.me.cell_cycle(
            adata_hvg, adata_scvi, batch_key="pt", organism="human", embed="X_scVI"
        ),
        scib.me.cell_cycle(
            adata_hvg, adata_scvi, batch_key="side", organism="human", embed="X_scVI"
        ),
    ]

    unintegrated_metrics_hvg = [
        scib.me.hvg_overlap(adata_filtered, adata_hvg, batch_key="chemistry"),
        scib.me.hvg_overlap(adata_filtered, adata_hvg, batch_key="pt"),
        scib.me.hvg_overlap(adata_filtered, adata_hvg, batch_key="side"),
    ]
    harmony_metrics_hvg = [
        scib.me.hvg_overlap(adata_hvg, adata_harmony, batch_key="chemistry"),
        scib.me.hvg_overlap(adata_hvg, adata_harmony, batch_key="pt"),
        scib.me.hvg_overlap(adata_hvg, adata_harmony, batch_key="side"),
    ]
    scvi_metrics_hvg = [
        scib.me.hvg_overlap(adata_hvg, adata_scvi, batch_key="chemistry"),
        scib.me.hvg_overlap(adata_hvg, adata_scvi, batch_key="pt"),
        scib.me.hvg_overlap(adata_hvg, adata_scvi, batch_key="side"),
    ]

    metrics_dict = {
        "Unintegrated": [
            unintegrated_metrics_pcr,
            unintegrated_metrics_cellc,
            unintegrated_metrics_hvg,
        ],
        "harmony": [harmony_metrics_pcr, harmony_metrics_cellc, harmony_metrics_hvg],
        "scvi": [scvi_metrics_pcr, scvi_metrics_cellc, scvi_metrics_hvg],
    }
    metrics = pd.DataFrame(metrics_dict)

    df = pd.DataFrame()

    # Iterate through the dictionary and create columns for each element
    for method, metrics_list in metrics_dict.items():
        for i, metrics in enumerate(metrics_list):
            col_name = f"{method}_metrics_{['pcr', 'cellc', 'hvg'][i]}"
            df[col_name] = metrics

    df = df.T
    df.columns = ["chemistry", "pt", "side"]
    df.style.background_gradient(cmap="Blues")
    df.to_csv(
        input_path
        / f"{tissue}_{n_top_genes}{method_hvg}_integration_comparison_{timestamp}.csv"
    )

    logging.info("########## Leiden clustering and saving ##########")
    for adata in [adata_harmony, adata_scvi]:
        sc.tl.leiden(adata, key_added="leiden_res0_5", resolution=0.5)
        sc.tl.leiden(adata, key_added="leiden_res1", resolution=1.0)
        sc.tl.leiden(adata, key_added="leiden_res1_25", resolution=1.25)
        sc.tl.leiden(adata, key_added="leiden_res1_5", resolution=1.5)
        sc.tl.leiden(adata, key_added="leiden_res1_75", resolution=1.75)
        sc.pl.umap(
            adata,
            color=[
                "leiden_res0_5",
                "leiden_res1",
                "leiden_res1_25",
                "leiden_res1_5",
                "leiden_res1_75",
            ],
            legend_loc="on data",
        )

    adata_harmony.write(
        output_path
        / tissue
        / f"{tissue}_{n_top_genes}{method_hvg}_harmony_clustered.h5ad"
    )
    adata_scvi.write(
        output_path / tissue / f"{tissue}_{n_top_genes}{method_hvg}_scvi_clustered.h5ad"
    )
