import argparse
import os
import warnings
from pathlib import Path

import cell2cell as c2c
import liana as li
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import statsmodels.stats.multitest as smm
import tensorly as tl
from tqdm import tqdm

from adata_processing_utils import (create_blood_brain_tf,
                                    filter_and_process_anno_adata,
                                    load_brain_blood_data)
from lr_loadings_utils import find_threshold_loadings

warnings.filterwarnings("ignore")

#### REQUIRES GPU AND LIANA ENV

tl.set_backend("pytorch")

np.random.seed(42)


def main(
    project_path,
    batch,
    sample_key,
    groupby,
    factorization_type,
    merging_mode,
    factorization_rank,
):
    project_path = Path(project_path)
    output_path = project_path / batch / "c2c_liana_outputs" / factorization_type
    os.makedirs(output_path, exist_ok=True)

    if factorization_type == "micro_blood_coarse":
        adata = create_blood_brain_tf(batch, project_path, write=False, micro_only=True)

    if factorization_type == "brain_blood_coarse":
        adata = create_blood_brain_tf(
            batch, project_path, write=False, micro_only=False
        )

        adata.obs["cell_type.v3"] = (
            adata.obs["cell_type_coarse"].astype(str)
            + "_"
            + adata.obs["tissue"].astype(str)
        )

    elif factorization_type == "brain_coarse":

        blood_adata, brain_adata, blood_adata_hvg, brain_adata_hvg = (
            load_brain_blood_data(batch, project_path)
        )
        adata = filter_and_process_anno_adata(brain_adata, brain_adata_hvg)

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
        lr_fill=np.nan,
        cell_fill=np.nan,
        outer_fraction=1 / 3.0,
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
        print("Compute rank")
        device = "cuda"
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

    output_folder = output_path / f"rank_{factorization_rank}" / merging_mode
    os.makedirs(output_folder, exist_ok=True)
    tensor_factorized = c2c.analysis.run_tensor_cell2cell_pipeline(
        tensor,
        meta_tensor,
        copy_tensor=True,
        rank=int(factorization_rank),
        tf_optimization="robust",
        random_state=0,
        device="cuda",
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

    factor, _ = c2c.plotting.tensor_factors_plot(
        interaction_tensor=tensor_factorized,
        reorder_elements={
                    "Sender Cells": [
                        "B.cells_blood",
                        "Monocytes_blood",
                        "NKT_blood",
                        "NK_blood",
                        "T-cells_blood",
                        "Astrocyte_brain",
                        "Microglia_brain",
                        "Epithelial_brain",
                        "Neurons_brain",
                        "OPC_brain",
                        "Oligo_brain",
                    ],
                    "Receiver Cells": [
                        "B.cells_blood",
                        "Monocytes_blood",
                        "NKT_blood",
                        "NK_blood",
                        "T-cells_blood",
                        "Astrocyte_brain",
                        "Microglia_brain",
                        "Epithelial_brain",
                        "Neurons_brain",
                        "OPC_brain",
                        "Oligo_brain",
                    ],
                },
        # reorder_elements={
        #     "Sender Cells": [
        #         "B.cells_blood",
        #         "CD16.NK_blood",
        #         "DC_blood",
        #         "Monocytes_blood",
        #         "NKT_blood",
        #         "Platelets_blood",
        #         "T-cells_blood",
        #         "CD16.NK_brain",
        #         "Microglia_brain",
        #         "Monocytes_brain",
        #         "Pericyte_brain",
        #         "T-cells_brain",
        #     ],
        #     "Receiver Cells": [
        #         "B.cells_blood",
        #         "CD16.NK_blood",
        #         "DC_blood",
        #         "Monocytes_blood",
        #         "NKT_blood",
        #         "Platelets_blood",
        #         "T-cells_blood",
        #         "CD16.NK_brain",
        #         "Microglia_brain",
        #         "Monocytes_brain",
        #         "Pericyte_brain",
        #         "T-cells_brain",
        #     ],
        # },
        metadata=meta_tensor,
        sample_col="Element",
        group_col="Category",
        meta_cmaps=["plasma", "Dark2_r", "tab20", "tab20"],
        fontsize=10,
        plot_legend=True,
        filename=output_folder / f"Factorization_{batch}.pdf",
    )

    factors = c2c.io.load_tensor_factors(output_folder / "Loadings.xlsx")

    context_dict = adata.obs.set_index("sample")["chemistry"].sort_values().to_dict()

    version_colors = c2c.plotting.aesthetics.get_colors_from_labels(
        ["v3", "v2"], cmap="plasma"
    )
    version_colors
    color_dict = {k: version_colors[v] for k, v in context_dict.items()}
    col_colors = pd.Series(color_dict)
    col_colors = col_colors.to_frame()
    col_colors.columns = ["Chemistry"]

    df_context = factors["Contexts"].copy()
    df_context["chemistry"] = context_dict.values()
    df_context.index = df_context.index.astype(str)

    sample_cm = c2c.plotting.loading_clustermap(
        df_context.iloc[:, :-1],
        use_zscore=False,
        col_colors=col_colors,
        figsize=(16, 6),
        dendrogram_ratio=0.3,
        cbar_fontsize=12,
        tick_fontsize=14,
        filename=output_folder / f"clustermap_{batch}.pdf",
    )

    plt.sca(sample_cm.ax_heatmap)
    legend = c2c.plotting.aesthetics.generate_legend(
        color_dict=version_colors,
        loc="upper right",
        bbox_to_anchor=(-0.15, 0.5),
        title="CellRanger version",
    )

    lr_cm = c2c.plotting.loading_clustermap(
        factors["Ligand-Receptor Pairs"],
        loading_threshold=0.1,
        use_zscore=False,
        figsize=(35, 8),
        filename=output_folder / f"Clustermap-LRs_{batch}.pdf",
        row_cluster=False,
    )

    for selected_factor in tqdm(
        factors["Contexts"].columns.tolist(), desc="Factor loadings plots"
    ):

        loading_product = c2c.analysis.tensor_downstream.get_joint_loadings(
            factors,
            dim1="Sender Cells",
            dim2="Receiver Cells",
            factor=selected_factor,
        )
        lprod_cm = c2c.plotting.loading_clustermap(
            loading_product.T,
            use_zscore=False,
            figsize=(8, 8),
            filename=output_folder
            / f"Clustermap_CC_{selected_factor.replace(' ','-')}.pdf",
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

        lsr_prod_cm = c2c.plotting.loading_clustermap(
            lr_cell_product,
            use_zscore=False,
            figsize=(15, 24),
            filename=output_folder
            / f"Clustermap_LR-CC_{selected_factor.replace(' ','-')}.pdf",
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p",
        "--project_path",
        help="Path were data are stored",
        default="/sc/arion/projects/psychgen/lbp/data/ScProcesses_brainBlood_nicole/",
    )
    parser.add_argument(
        "-b", "--batch", help="Batch identifier", choices=["batch_1", "batch_2"]
    )
    parser.add_argument(
        "-s", "--sample_key", help="Patient column in adata.obs", default="sample"
    )
    parser.add_argument(
        "-g",
        "--groupby",
        help="Celltype column in adata.obs for liana",
        default="cell_type_coarse",
    )
    parser.add_argument(
        "-f",
        "--factorization_type",
        help="Strategy for tensor factorization. This name is used to create out folder",
    )
    parser.add_argument(
        "-m", "--merging_mode", help="Merging mode for c2c", default="inner"
    )
    parser.add_argument(
        "-r", "--factorization_rank", help="c2c factorization rank", default=None
    )
    args = parser.parse_args()
    main(
        args.project_path,
        args.batch,
        args.sample_key,
        args.groupby,
        args.factorization_type,
        args.merging_mode,
        args.factorization_rank,
    )
