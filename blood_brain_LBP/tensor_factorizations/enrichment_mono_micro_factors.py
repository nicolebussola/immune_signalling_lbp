import argparse
import gc
import logging
import math
import os
import warnings
from pathlib import Path

import cell2cell as c2c
import decoupler as dc
import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pygal
import scanpy as sc
import seaborn as sns
from gseapy import barplot, dotplot

from adata_processing_utils import (create_blood_brain_tf,
                                    prepare_scran_de_data_micro_mono,
                                    scran_norm)
from enrichment_analysis_utils import (enr_df, get_expressed_genes,
                                       process_common_paths)
from lr_loadings_utils import (extract_ligands_receptors,
                               filter_de_genes_micro, find_threshold_loadings,
                               lr_from_micro_to_mono, lr_from_mono_to_micro,
                               process_loadings)

logging.disable(logging.WARNING)


warnings.filterwarnings("ignore")

CELLS_FACTOR_DICT = {
    "micro_blood_coarse": {
        "b1": {"Factor 3": "To Microglia", "Factor 8": "From Microglia"},
        "b2": {"Factor 3": "To Microglia", "Factor 5": "From Microglia"},
    },
    "brain_coarse": {"b1": {"Factor 3": "To Microglia", "Factor 4": "From Microglia"}},
    "brain_blood_coarse": {
        "b1": {"Factor 12": "From Microglia", "Factor 15": "To Microglia"}
    },
}
BATCH_DICT = {"batch_1": "b1", "batch_2": "b2"}

BLOOD_SCRAN_NORMALIZED_FNAME = "data_blood_scran_de_micro_mono.h5ad"
BRAIN_SCRAN_NORMALIZED_FNAME = "data_brain_scran_de_micro_mono.h5ad"


def main(
    project_path,
    batch,
    merging_mode,
    factorization_type,
    factorization_rank,
    tensorized_data_fname,
):
    """
    Main function to execute the data processing workflow.

    Parameters:
    project_path (str): The path to the project directory.
    batch (str): The batch identifier for the data to be processed.
    merging_mode (str): The mode of merging data for tensor factorization as defined by c2c.
    factorization_type (str): The type of factorization used, note that this corresponds to the folder where tensor is saved
    factorization_rank (int): The rank or number of components used in c2c factorization.

    Returns:
    None
    """

    project_path = Path(project_path)
    data_path = (
        project_path
        / batch
        / "c2c_liana_outputs"
        / factorization_type
        / f"rank_{factorization_rank}"
        / merging_mode
    )
    out_path = data_path / "enrichr_results"
    os.makedirs(out_path, exist_ok=True)

    print(f"\n Data: {factorization_type} - {batch}\n")

    print("\nLoad and clean tensor Loadings\n")

    tensor_factors = c2c.io.load_tensor_factors(data_path / "Loadings.xlsx")

    for k in tensor_factors.keys():
        tensor_factors[k] = tensor_factors[k].rename(
            CELLS_FACTOR_DICT[factorization_type][BATCH_DICT[batch]], axis=1
        )

    lr_loadings = tensor_factors["Ligand-Receptor Pairs"]
    thresh_loadings = find_threshold_loadings(tensor_factors, quantile=0.70)
    lr_loadings_from_micro = process_loadings(
        lr_loadings, "From Microglia", thresh_loadings
    ).reset_index()
    lr_loadings_to_micro = process_loadings(
        lr_loadings, "To Microglia", thresh_loadings
    ).reset_index()

    if os.path.isfile(
        project_path / batch / "c2c_liana_outputs" / BLOOD_SCRAN_NORMALIZED_FNAME
    ) and os.path.isfile(
        project_path / batch / "c2c_liana_outputs" / BRAIN_SCRAN_NORMALIZED_FNAME
    ):
        print("Data already created")
        adata_bl = sc.read(
            project_path
            / "batch_1"
            / "c2c_liana_outputs"
            / BLOOD_SCRAN_NORMALIZED_FNAME
        )
        adata_br = sc.read(
            project_path
            / "batch_1"
            / "c2c_liana_outputs"
            / BRAIN_SCRAN_NORMALIZED_FNAME
        )
    else:
        print("Scran and data preparation")
        adata_bl, adata_br = prepare_scran_de_data_micro_mono(
            batch, project_path, "cell_type.v2", "cell_type.v2"
        )

    n_celltypes_brain = adata_br.obs["cell_type.v2"].nunique()
    n_celltypes_blood = adata_bl.obs["cell_type.v2"].nunique()

    all_from_microglia_ligands, _ = extract_ligands_receptors(lr_loadings_from_micro)
    all_to_microglia_ligands, _ = extract_ligands_receptors(lr_loadings_to_micro)

    if factorization_type == "brain_coarse":
        adata_receptors = adata_br.copy()
        adata_ligands = adata_br.copy()
        n_celltypes_ligands = n_celltypes_receptors = n_celltypes_brain

    else:
        adata_ligands = adata_br.copy()
        adata_receptors = adata_bl.copy()
        n_celltypes_ligands = n_celltypes_brain
        n_celltypes_receptors = n_celltypes_blood

    lr_from_micro_to_mono_df = lr_from_micro_to_mono(
        adata_ligands=adata_ligands,
        adata_receptors=adata_receptors,
        lr_loadings=lr_loadings_from_micro,
        genes_ligands=all_from_microglia_ligands,
        n_celltypes_ligands=n_celltypes_ligands,
        n_celltypes_receptors=n_celltypes_receptors,
    )

    lr_from_mono_to_micro_df = lr_from_mono_to_micro(
        adata_ligands=adata_receptors,
        adata_receptors=adata_ligands,
        lr_loadings=lr_loadings_to_micro,
        genes_ligands=all_to_microglia_ligands,
        n_celltypes_ligands=n_celltypes_receptors,
        n_celltypes_receptors=n_celltypes_ligands,
    )

    if factorization_type == "micro_blood_coarse":
        print("\nCheck and filter genes on factorized data\n")
        # try:
        adata_tf = sc.read(
            project_path / batch / "c2c_liana_outputs" / tensorized_data_fname
        )
        # else:
        # adata_ccc = create_blood_brain_tf(batch, project_path, write=False)

        if "Microglia" not in adata_ccc.uns.keys():
            adata_tf.X = adata_tf.X.toarray()
            adata_tf = scran_norm(adata_tf, "cell_type_coarse")
            sc.tl.rank_genes_groups(
                adata_tf,
                "cell_type_coarse",
                reference="Microglia",
                method="wilcoxon",
                layer="scran_normalization",
                pts=True,
                key_added="Microglia",
            )
            sc.tl.rank_genes_groups(
                adata_tf,
                "cell_type_coarse",
                reference="Monocytes",
                method="wilcoxon",
                layer="scran_normalization",
                pts=True,
                key_added="Monocytes",
            )

        lr_from_micro_to_mono_df = filter_de_genes_micro(
            adata_ccc,
            lr_from_micro_to_mono_df["receptor"].unique().tolist(),
            lr_from_micro_to_mono_df,
            "receptor",
        )

        lr_from_mono_to_micro_df = filter_de_genes_micro(
            adata_ccc,
            lr_from_mono_to_micro_df["ligand"].unique().tolist(),
            lr_from_mono_to_micro_df,
            "ligand",
        )

    lr_from_micro_to_mono_df.to_csv(
        data_path / f"filtered_lr_loadings_{BATCH_DICT[batch]}_from_micro.csv",
        index=None,
    )
    lr_from_mono_to_micro_df.to_csv(
        data_path / f"filtered_lr_loadings_{BATCH_DICT[batch]}_to_micro.csv",
        index=None,
    )

    print("\nEnrichment analysis\n")

    gset_ligands_from = lr_from_micro_to_mono_df["ligand"].unique().tolist()
    gset_receptors_from = lr_from_micro_to_mono_df["receptor"].unique().tolist()

    gset_ligands_to = lr_from_mono_to_micro_df["ligand"].unique().tolist()
    gset_receptors_to = lr_from_mono_to_micro_df["receptor"].unique().tolist()

    gsea_datasets = ["WikiPathway_2023_Human", "Reactome_2022"]
    background_micro = get_expressed_genes(
        adata_br, "cell_type_coarse", "Microglia", expr_prop=0.1
    )
    background_mono = get_expressed_genes(
        adata_bl, "cell_type_coarse", "Monocytes", expr_prop=0.1
    )

    enr_df_ligands_from = enr_df(
        gset_ligands_from, gsea_datasets[0], background_micro, out_path
    )
    enr_df_receptors_from = enr_df(
        gset_receptors_from, gsea_datasets[0], background_mono, out_path
    )

    enr_df_ligands_to = enr_df(
        gset_ligands_to, gsea_datasets[1], background_mono, out_path
    )
    enr_df_receptors_to = enr_df(
        gset_receptors_to, gsea_datasets[1], background_micro, out_path
    )

    print("\nFind ligand-receptor common pathways\n")

    common_pathways_from = process_common_paths(
        enr_df_ligands_from, enr_df_receptors_from
    )
    common_pathways_to = process_common_paths(enr_df_ligands_to, enr_df_receptors_to)

    for direction, common_pathways_df, gsea_data in zip(
        ["from", "to"], [common_pathways_from, common_pathways_to], gsea_datasets
    ):
        barplot(
            common_pathways_df,
            column="Adjusted P-value",
            group="Gene_set",  # set group, so you could do a multi-sample/library comparsion
            size=10,
            top_term=common_pathways_df.shape[0],
            figsize=(3, 5),
            color={"WikiPathway_2023_Human": "#F6AE2D", "Reactome_2022": "#64D3F1"},
            ofname=out_path / f"{gsea_data}_enriched_{direction}_barplot.pdf",
        )

        dotplot(
            common_pathways_df,
            column="Adjusted P-value",
            x="Gene_set",
            size=5,
            top_term=common_pathways_df.shape[0],
            figsize=(4, 10),
            title=f"Batch {batch} - {factorization_type}",
            xticklabels_rot=0,
            show_ring=False,
            marker="o",
            ofname=out_path / f"{gsea_data}_enriched_{direction}_dotplot.pdf",
        )

    print("\nPseudobulk expression products for LR pairs\n")

    sample_key = "pt"
    groupby = "cell_type.v2"
    pdata_br = dc.get_pseudobulk(
        adata_br,
        sample_col=sample_key,
        groups_col=groupby,
        layer="counts",
        mode="sum",
        min_cells=10,
        min_counts=10000,
    )
    pdata_bl = dc.get_pseudobulk(
        adata_bl,
        sample_col=sample_key,
        groups_col=groupby,
        layer="counts",
        mode="sum",
        min_cells=10,
        min_counts=10000,
    )
    df_micro = pd.DataFrame(
        pdata_br[pdata_br.obs[groupby] == "Microglia"].X.mean(axis=0)
    ).T
    df_micro.columns = list(pdata_br.var_names)
    df_micro = df_micro.rename(index={0: "Microglia"})

    if factorization_type == "brain_coarse":
        padata_receptor = pdata_br.copy()
    else:
        padata_receptor = pdata_bl.copy()

    for direction, lr_pairs in zip(
        ["from", "to"],
        [
            lr_from_micro_to_mono_df["index"].to_list(),
            lr_from_mono_to_micro_df["index"].to_list(),
        ],
    ):
        prod_dict = {}
        for lr_pair in lr_pairs:
            ligand = lr_pair.split("^")[0]
            receptor = lr_pair.split("^")[1]
            prod_dict[lr_pair] = {}
            receptor_cell_types = [
                cl
                for cl in padata_receptor.obs["cell_type.v2"].unique()
                if cl != "Microglia"
            ]
            for cell_type in receptor_cell_types:
                df_c = pd.DataFrame(
                    padata_receptor[
                        padata_receptor.obs["cell_type.v2"] == cell_type
                    ].X.mean(axis=0)
                ).T
                df_c.columns = list(padata_receptor.var_names)
                df_c = df_c.rename(index={0: cell_type})
                prod_dict[lr_pair][cell_type] = df_micro[ligand][0] * df_c[receptor][0]

        dot_chart = pygal.Dot(x_label_rotation=30)

        dot_chart.title = "LR pairs expression products"

        prod_dict_df = pd.DataFrame(prod_dict)
        dot_chart.x_labels = prod_dict_df.columns
        for cell_type in prod_dict_df.index:
            dot_chart.add(cell_type, prod_dict_df.loc[cell_type])

        dot_chart.render_to_png(
            str(out_path / f"filtered_lr_product_dotplot_{direction}.png")
        )


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
        "-m", "--merging_mode", help="Merging mode for c2c", default="inner"
    )
    parser.add_argument(
        "-f", "--factorization_type", help="Strategy for tensor factorization"
    )
    parser.add_argument(
        "-r", "--factorization_rank", help="c2c factorization rank", default=8
    )
    parser.add_argument(
        "-td",
        "--tensorized_data_fname",
        help="data used for c2c factorization",
        default="blood_microglia_tf.h5ad",
    )
    args = parser.parse_args()
    main(
        args.project_path,
        args.batch,
        args.merging_mode,
        args.factorization_type,
        args.factorization_rank,
        args.tensorized_data_fname,
    )
