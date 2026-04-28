"""
Step 1: Tensor factorization.

Runs LIANA rank-aggregate by sample and Tensor-cell2cell factorization.
Supports three analytical designs for human cohorts:
  - brain_coarse:       brain cell types only (Cohort 1 or 2)
  - micro_blood_coarse: microglia (brain) + all blood cell types
  - brain_blood_coarse: all brain + all blood cell types

Run as a standalone step:
    python -m tensorS2R.step1_factorize -p <path> -b <cohort> -f <factorization_type> [options]
"""

import os
import warnings
from pathlib import Path

import cell2cell as c2c
import liana as li
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tensorly as tl
from tqdm import tqdm

from .adata_processing_utils import (
    create_blood_brain_tf,
    filter_and_process_anno_adata,
    load_brain_blood_data,
)
from .lr_loadings_utils import find_threshold_loadings

warnings.filterwarnings("ignore")
np.random.seed(42)


def run_tensorcell2cell(
    project_path,
    cohort,
    sample_key,
    groupby,
    factorization_type,
    merging_mode,
    factorization_rank,
    device,
):
    """
    Run LIANA rank-aggregate + Tensor-cell2cell factorization for one analytical design.

    Returns:
        output_folder (Path): Directory where all factorization outputs are written.
        adata: The AnnData object used for factorization (useful for downstream steps).
    """
    project_path = Path(project_path)
    output_path = project_path / cohort / "c2c_liana_outputs" / factorization_type
    os.makedirs(output_path, exist_ok=True)

    if device == "cuda":
        tl.set_backend("pytorch")

    if factorization_type == "micro_blood_coarse":
        adata = create_blood_brain_tf(cohort, project_path, write=False, micro_only=True)
    elif factorization_type == "brain_blood_coarse":
        adata = create_blood_brain_tf(cohort, project_path, write=False, micro_only=False)
        adata.obs["cell_type.v3"] = (
            adata.obs["cell_type_coarse"].astype(str) + "_" + adata.obs["tissue"].astype(str)
        )
        adata.obs["cell_type.v4"] = (
            adata.obs["cell_type"].astype(str) + "_" + adata.obs["tissue"].astype(str)
        )
    elif factorization_type == "brain_coarse":
        _, brain_adata, _, brain_adata_hvg = load_brain_blood_data(cohort, project_path)
        adata = filter_and_process_anno_adata(brain_adata, brain_adata_hvg)
    else:
        raise ValueError(
            f"Unknown factorization_type: '{factorization_type}'. "
            "Choose from: micro_blood_coarse, brain_blood_coarse, brain_coarse."
        )

    adata.X = adata.layers["log1p_norm"].copy()
    adata.obs[sample_key] = adata.obs["pt"]

    li.mt.rank_aggregate.by_sample(
        adata=adata,
        groupby=groupby,
        sample_key=sample_key,
        expr_prop=0.1,
        use_raw=False,
        n_perms=100,
        return_all_lrs=False,
        verbose=True,
    )

    target_celltypes = list(adata.obs[groupby].unique())

    p = li.pl.dotplot(
        adata=adata,
        colour="magnitude_rank",
        size="specificity_rank",
        inverse_colour=True,
        inverse_size=True,
        source_labels=target_celltypes,
        target_labels=target_celltypes,
        filter_fun=lambda x: x["specificity_rank"] <= 0.05,
        orderby="magnitude_rank",
        orderby_ascending=True,
        top_n=20,
        figure_size=(20, 20),
        size_range=(1, 6),
    )
    p.save(output_path / "liana_agg.pdf")

    sorted_samples = sorted(adata.uns["liana_res"]["sample"].unique())

    tensor = li.multi.to_tensor_c2c(
        liana_res=adata.uns["liana_res"],
        sample_key="sample",
        source_key="source",
        target_key="target",
        ligand_key="ligand_complex",
        receptor_key="receptor_complex",
        score_key="magnitude_rank",
        non_negative=True,
        inverse_fun=lambda x: 1 - x,
        non_expressed_fill=None,
        how=merging_mode,
        lr_sep="^",
        context_order=sorted_samples,
        sort_elements=True,
    )

    context_dict = adata.obs.set_index("sample")["chemistry"].sort_values().to_dict()
    meta_tensor = c2c.tensor.generate_tensor_metadata(
        interaction_tensor=tensor,
        metadata_dicts=[context_dict, None, None, None],
        fill_with_order_elements=True,
    )

    if factorization_rank is None:
        if len(target_celltypes) > 7:
            factorization_rank = 20
        else:
            print("Computing optimal rank via elbow selection...")
            tensor.to_device(device)
            fig, error = tensor.elbow_rank_selection(
                upper_rank=25,
                runs=100,
                init="random",
                svd="numpy_svd",
                automatic_elbow=True,
                metric="error",
                random_state=0,
                fontsize=14,
                filename=str(output_path),
                tol=1e-8,
                n_iter_max=500,
            )
            factorization_rank = tensor.rank

    output_folder = output_path / f"rank_{factorization_rank}" / merging_mode / groupby
    os.makedirs(output_folder, exist_ok=True)

    tensor_factorized = c2c.analysis.run_tensor_cell2cell_pipeline(
        tensor,
        meta_tensor,
        copy_tensor=True,
        rank=int(factorization_rank),
        tf_optimization="robust",
        random_state=0,
        device=device,
        elbow_metric="error",
        smooth_elbow=False,
        upper_rank=25,
        tf_init="random",
        tf_svd="numpy_svd",
        cmaps=None,
        sample_col="Element",
        group_col="Category",
        fig_fontsize=14,
        output_folder=str(output_folder),
        output_fig=False,
    )

    c2c.plotting.tensor_factors_plot(
        interaction_tensor=tensor_factorized,
        metadata=meta_tensor,
        sample_col="Element",
        group_col="Category",
        meta_cmaps=["plasma", "Dark2_r", "tab20", "tab20"],
        fontsize=10,
        plot_legend=True,
        filename=output_folder / f"Factorization_{cohort}.pdf",
    )

    factors = c2c.io.load_tensor_factors(output_folder / "Loadings.xlsx")

    version_colors = c2c.plotting.aesthetics.get_colors_from_labels(["v3", "v2"], cmap="plasma")
    color_dict = {k: version_colors[v] for k, v in context_dict.items()}
    col_colors = pd.Series(color_dict).to_frame()
    col_colors.columns = ["Chemistry"]

    df_context = factors["Contexts"].copy()
    df_context["chemistry"] = list(context_dict.values())
    df_context.index = df_context.index.astype(str)

    sample_cm = c2c.plotting.loading_clustermap(
        df_context.iloc[:, :-1],
        use_zscore=False,
        col_colors=col_colors,
        figsize=(16, 6),
        dendrogram_ratio=0.3,
        cbar_fontsize=12,
        tick_fontsize=14,
        filename=output_folder / f"clustermap_{cohort}.pdf",
    )
    plt.sca(sample_cm.ax_heatmap)
    c2c.plotting.aesthetics.generate_legend(
        color_dict=version_colors,
        loc="upper right",
        bbox_to_anchor=(-0.15, 0.5),
        title="CellRanger version",
    )

    c2c.plotting.loading_clustermap(
        factors["Ligand-Receptor Pairs"],
        loading_threshold=0.1,
        use_zscore=False,
        figsize=(35, 8),
        filename=output_folder / f"Clustermap-LRs_{cohort}.pdf",
        row_cluster=False,
    )

    for selected_factor in tqdm(factors["Contexts"].columns.tolist(), desc="Factor loading plots"):
        loading_product = c2c.analysis.tensor_downstream.get_joint_loadings(
            factors, dim1="Sender Cells", dim2="Receiver Cells", factor=selected_factor,
        )
        c2c.plotting.loading_clustermap(
            loading_product.T,
            use_zscore=False,
            figsize=(8, 8),
            filename=output_folder / f"Clustermap_CC_{selected_factor.replace(' ', '-')}.pdf",
            cbar_label="Loading Product",
        )
        lr_cell_product = c2c.analysis.tensor_downstream.get_lr_by_cell_pairs(
            factors,
            lr_label="Ligand-Receptor Pairs",
            sender_label="Sender Cells",
            receiver_label="Receiver Cells",
            order_cells_by="receivers",
            factor=selected_factor,
            cci_threshold=None,
            lr_threshold=0.05,
        )
        c2c.plotting.loading_clustermap(
            lr_cell_product,
            use_zscore=False,
            figsize=(15, 24),
            filename=output_folder / f"Clustermap_LR-CC_{selected_factor.replace(' ', '-')}.pdf",
            cbar_label="Loading Product",
            yticklabels=1,
        )

    thresh_loadings = find_threshold_loadings(factors, quantile=0.70)
    c2c.plotting.ccc_networks_plot(
        factors,
        included_factors=factors["Contexts"].columns.tolist(),
        ccc_threshold=thresh_loadings,
        nrows=1,
        panel_size=(16, 16),
        filename=output_folder / "Networks.pdf",
    )
    c2c.analysis.tensor_downstream.compute_gini_coefficients(
        factors, sender_label="Sender Cells", receiver_label="Receiver Cells"
    ).sort_values("Gini").to_csv(output_folder / "Gini.csv")

    return output_folder, adata


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Step 1: Tensor factorization (LIANA + Tensor-cell2cell)")
    parser.add_argument("-p", "--project_path", required=True, help="Root data directory")
    parser.add_argument("-b", "--cohort", required=True, choices=["cohort_1", "cohort_2"])
    parser.add_argument("-s", "--sample_key", default="sample", help="Patient column in adata.obs")
    parser.add_argument("-g", "--groupby", default="cell_type_coarse", help="Cell-type column for LIANA")
    parser.add_argument(
        "-f", "--factorization_type", required=True,
        choices=["brain_coarse", "micro_blood_coarse", "brain_blood_coarse"],
    )
    parser.add_argument("-m", "--merging_mode", default="inner")
    parser.add_argument("-r", "--factorization_rank", type=int, default=None)
    parser.add_argument("-d", "--device", default="cuda", choices=["cuda", "cpu"])
    args = parser.parse_args()

    run_tensorcell2cell(
        args.project_path,
        args.cohort,
        args.sample_key,
        args.groupby,
        args.factorization_type,
        args.merging_mode,
        args.factorization_rank,
        args.device,
    )
