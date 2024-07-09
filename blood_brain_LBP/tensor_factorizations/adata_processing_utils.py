import gc
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc
import scranPY as sp
from scipy.sparse import csr_matrix

cell_type_mapping = {
    "DC2/DC3": "DC",
    "DC3": "DC",
    "DC1": "DC",
    "pDC": "DC",
    "B cells": "B.cells",
    "Plasma cells": "Plasma",
    "CD14+ Monocytes": "CD14.mono",
    "CD16+ Monocytes": "CD16.mono",
    "CD14+ Mono": "CD14.mono",
    "CD16+ Mono": "CD16.mono",
    "CD4+ T cells": "CD4.T",
    "CD4+ T naive": "CD4.T",
    "CD8+ T cells": "CD8.T",
    "CD8+ NKT cells": "NKT",
    "CD16+ NKT": "NKT",
    # "CD16+ NK": "CD16.NK",
    "CD16+ NK": "NK",
    # "CD16+ NK cells": "CD16.NK",
     "CD16+ NK cells": "NK",
     "CD56+ NK cells": "NK",
    "Neuron-low count": "Neurons",
    "GABAergic Neuron": "Neurons",
    "Glutamatergic Neuron": "Neurons",
    "Neuron-low count": "Neurons",
}

relevant_obs_columns = [
    "n_genes_by_counts",
    "log1p_n_genes_by_counts",
    "total_counts",
    "log1p_total_counts",
    "pct_counts_in_top_20_genes",
    "total_counts_mt",
    "log1p_total_counts_mt",
    "pct_counts_mt",
    "total_counts_ribo",
    "log1p_total_counts_ribo",
    "pct_counts_ribo",
    "predicted_doublets_consensus",
    "side",
    "pt",
    "pt-side",
    "tissue",
    "chemistry",
    "cell_type",
    "cell_type.v2",
    "cell_type_coarse",
]


ADATA_BATCH_1_BLOOD_HVG = (
    "batch_1_blood_filtered_4000HighlyDeviant_20_harmony_annotated.h5ad"
)
ADATA_BATCH_1_BRAIN_HVG = (
    "batch_1_brain_filtered_4000HighlyDeviant_20_harmony_regressedScaled.h5ad"
)

ADATA_BATCH_2_BLOOD_HVG = (
    "batch_2_blood_filtered_4000HighlyDeviant_20_filt15_scanvi_annotated.h5ad"
)
ADATA_BATCH_2_BRAIN_HVG = (
    "batch_2_brain_filtered_4000HighlyDeviant_20_harmony_scaled.h5ad"
)
BATCH_2_BRAIN_ANNOTATIONS_PATH = (
    "../annotations/batch_2_brain_filtered_4000HighlyDeviant_20_harmony_scaled.csv"
)


def load_brain_blood_data(batch, project_path):
    """
    Load brain and blood data for a specific batch;
    both data for full genes and and highly variable genes are loaded.

    Parameters:
    batch (str): The batch identifier (e.g., 'batch_1', 'batch_2').
    project_path (Path): Path to the project directory.

    Returns:
    tuple: Processed AnnData objects for blood and brain data (blood_adata, brain_adata, blood_adata_hvg, brain_adata_hvg).
    """
    data_path = project_path / batch
    blood_adata = sc.read(data_path / f"{batch}_blood_filtered_20.h5ad")
    brain_adata = sc.read(data_path / f"{batch}_brain_filtered_20.h5ad")

    if batch == "batch_1":
        blood_adata_hvg = sc.read(data_path / ADATA_BATCH_1_BLOOD_HVG)
        brain_adata_hvg = sc.read(data_path / ADATA_BATCH_1_BRAIN_HVG)
    elif batch == "batch_2":
        blood_adata_hvg = sc.read(data_path / ADATA_BATCH_2_BLOOD_HVG)
        blood_adata_hvg = blood_adata_hvg[
            ~blood_adata_hvg.obs["cell_type"].isin(["low count", "unknown"])
        ].copy()
        blood_adata_hvg.obs = blood_adata_hvg.obs.set_index("Unnamed: 0")
        brain_adata_hvg = sc.read(data_path / ADATA_BATCH_2_BRAIN_HVG)
        brain_adata_hvg.obs["cell_type"] = list(
            pd.read_csv(BATCH_2_BRAIN_ANNOTATIONS_PATH)["cell_type"]
        )
        blood_adata.obs["batch"] = "batch_2"
        blood_adata.obs["chemistry"] = "v3"

    return blood_adata, brain_adata, blood_adata_hvg, brain_adata_hvg


def filter_adata(adata, hvg_adata, thresh_count=7, min_cells=100):
    """
    Filter AnnData object based on total counts and cell type occurrences.

    Parameters:
    adata (AnnData): The AnnData object to filter.
    hvg_adata (AnnData): The highly variable genes AnnData object.
    thresh_count (int, optional): Threshold for filtering based on log1p total counts. Default is 7.
    min_cells (int, optional): Minimum number of cells for a cell type to be retained. Default is 100.

    Returns:
    AnnData: Filtered AnnData object.
    """
    adata_filt = adata[adata.obs["log1p_total_counts"] > thresh_count].copy()
    hvg_filt = hvg_adata[hvg_adata.obs["log1p_total_counts"] > thresh_count].copy()
    cell_type_counts = hvg_filt.obs["cell_type"].value_counts()
    cell_types_min_cells = cell_type_counts[cell_type_counts > min_cells].index.tolist()
    hvg_filt = hvg_filt[hvg_filt.obs["cell_type"].isin(cell_types_min_cells), :].copy()

    return adata_filt[adata_filt.obs.index.isin(hvg_filt.obs.index)].copy()


def create_blood_brain_tf(batch, project_path, write=False, micro_only=True):
    """
    Create a combined AnnData object of blood and brain cells, and optionally write to file.

    Parameters:
    batch (str): The batch identifier (e.g., 'batch_1', 'batch_2').
    project_path (Path): Path to the project directory.
    write (bool, optional): Whether to write the processed AnnData to file. Default is False.
    micro_only (bool, optional): Keep only Microglia cells in brain data.

    Returns:
    AnnData or None: Combined AnnData object if write is False, otherwise None.
    """
    blood_adata, brain_adata, blood_adata_hvg, brain_adata_hvg = load_brain_blood_data(
        batch, project_path
    )
    blood_adata_filt = filter_adata(blood_adata, blood_adata_hvg)
    del blood_adata
    gc.collect()
    brain_adata_filt = filter_adata(brain_adata, brain_adata_hvg)
    del brain_adata
    gc.collect()
    blood_adata_filt.obs["cell_type"] = blood_adata_hvg.obs["cell_type"].astype(str)
    brain_adata_filt.obs["cell_type"] = brain_adata_hvg.obs["cell_type"].astype(str)

    if micro_only:
        brain_type = "microglia"
        brain_adata_filt = brain_adata_filt[
            brain_adata_filt.obs["cell_type"] == "Microglia"
        ]
    else:
        brain_type = "brain"

    common_genes = set(brain_adata_filt.var.index).intersection(
        brain_adata_filt.var.index
    )

    blood_adata_filt_common = blood_adata_filt[
        :, blood_adata_filt.var.index.isin(common_genes)
    ]
    brain_adata_filt_common = brain_adata_filt[
        :, brain_adata_filt.var.index.isin(common_genes)
    ]

    adata = ad.concat(
        [blood_adata_filt_common, brain_adata_filt_common],
        join="inner",
        index_unique="_",
        merge="same",
    ).copy()

    del blood_adata_filt, brain_adata_filt
    gc.collect()

    adata.obs["cell_type.v2"] = adata.obs["cell_type"].astype(str)

    adata.obs["cell_type.v2"] = (
        adata.obs["cell_type.v2"].replace(cell_type_mapping).astype("category")
    )

    adata.obs["cell_type_coarse"] = (
        adata.obs["cell_type.v2"]
        .astype(str)
        .replace(
            {
                "CD14.mono": "Monocytes",
                "CD16.mono": "Monocytes",
                "CD4.T": "T-cells",
                "CD8.T": "T-cells",
            }
        )
    )

    adata.obs = adata.obs[relevant_obs_columns]

    if micro_only:
        adata = adata[
            adata.obs["cell_type.v2"].isin(
                [
                    "CD14.mono",
                    "CD16.mono",
                    "DC",
                    "NKT",
                    "B.cells",
                    "CD4.T",
                    "Microglia",
                    "CD16.NK",
                    "CD8.T",
                ]
            )
        ]
    adata.X = adata.layers["counts"].copy()
    scales_counts = sc.pp.normalize_total(adata, target_sum=1e4, inplace=False)
    adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)

    if write:
        output_path = project_path / batch / "c2c_liana_outputs"
        sc.write(output_path / f"blood_{brain_type}_tf.h5ad", adata)
        del adata
        gc.collect()
        print(f"Processed {batch} and written to {output_path}")
    else:
        return adata


def filter_and_process_anno_adata(adata, adata_hvg):
    """
    Filter and process AnnData object.

    Parameters:
    adata (AnnData): The AnnData object to filter and process.
    adata_hvg (AnnData): The highly variable genes AnnData object.

    Returns:
    AnnData: Filtered and processed AnnData object.
    """
    adata_filt = filter_adata(adata, adata_hvg)
    del adata
    gc.collect()
    adata_filt.obs["cell_type"] = adata_hvg.obs["cell_type"].astype(str)

    adata_filt.X = adata_filt.layers["counts"].copy()

    scales_counts = sc.pp.normalize_total(adata_filt, target_sum=1e4, inplace=False)
    adata_filt.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)

    del scales_counts

    adata_filt.obs["cell_type.v2"] = (
        adata_filt.obs["cell_type"]
        .astype(str)
        .replace(cell_type_mapping)
        .astype("category")
    )
    adata_filt.obs["cell_type_coarse"] = (
        adata_filt.obs["cell_type.v2"]
        .astype(str)
        .replace(
            {
                "CD14.mono": "Monocytes",
                "CD16.mono": "Monocytes",
                "CD4.T": "T-cells",
                "CD8.T": "T-cells",
            }
        )
    )
    adata_filt.X = adata_filt.X.toarray()

    return adata_filt


def scran_norm(adata, anno_col):
    """
    Normalize AnnData object using scran normalization method.

    Parameters:
    adata (AnnData): The AnnData object to normalize.
    anno_col (str): The annotation column to use for clustering during normalization.

    Returns:
    AnnData: Normalized AnnData object with scran normalization layer added.
    """
    sum_factors = sp.compute_sum_factors(adata, clusters=anno_col, plotting=False)
    scran = adata.X / adata.obs["size_factors"].values[:, None]
    adata.layers["scran_normalization"] = csr_matrix(sc.pp.log1p(scran))
    return adata


def prepare_scran_de_data_micro_mono(batch, project_path, anno_col_scran, anno_col_de):
    """
    Prepare scran normalized data for differential expression analysis between Microglia and Monocytes.

    Parameters:
    batch (str): The batch identifier (e.g., 'batch_1', 'batch_2').
    project_path (Path): Path to the project directory.
    anno_col_scran (str): The annotation column to use for scran normalization.
    anno_col_de (str): The annotation column to use for differential expression analysis.

    Returns:
    tuple: Scran normalized AnnData objects for blood and brain data (adata_bl, adata_br).
    """
    blood_adata, brain_adata, blood_adata_hvg, brain_adata_hvg = load_brain_blood_data(
        batch, project_path
    )

    adata_bl = scran_norm(
        filter_and_process_anno_adata(blood_adata, blood_adata_hvg), anno_col_scran
    )
    del blood_adata
    gc.collect()
    adata_br = scran_norm(
        filter_and_process_anno_adata(brain_adata, brain_adata_hvg), anno_col_scran
    )

    del brain_adata
    gc.collect()

    adata_br.X = adata_br.layers["scran_normalization"]
    adata_bl.X = adata_bl.layers["scran_normalization"]

    sc.tl.rank_genes_groups(
        adata_br,
        anno_col_de,
        reference="Microglia",
        method="wilcoxon",
        layer="scran_normalization",
        pts=True,
        key_added="Microglia",
    )
    sc.tl.rank_genes_groups(
        adata_br,
        anno_col_de,
        reference="CD16.mono",
        method="wilcoxon",
        layer="scran_normalization",
        pts=True,
        key_added="CD16",
    )
    sc.tl.rank_genes_groups(
        adata_br,
        anno_col_de,
        reference="CD14.mono",
        method="wilcoxon",
        layer="scran_normalization",
        pts=True,
        key_added="CD14",
    )
    sc.tl.rank_genes_groups(
        adata_bl,
        anno_col_de,
        reference="CD16.mono",
        method="wilcoxon",
        layer="scran_normalization",
        pts=True,
        key_added="CD16",
    )
    sc.tl.rank_genes_groups(
        adata_bl,
        anno_col_de,
        reference="CD14.mono",
        method="wilcoxon",
        layer="scran_normalization",
        pts=True,
        key_added="CD14",
    )
    return adata_bl, adata_br
