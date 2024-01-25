import logging
import urllib
import warnings

import anndata2ri
import celltypist
import numpy as np
import pandas as pd
import scanpy as sc
import scarches as sca
from celltypist import models
from rich.logging import RichHandler
from rpy2.robjects.packages import importr

from ..labels import (CHEMISTRY_V2_PATIENTS, FAILED_QC_SAMPLES,
                      QC_LABELS_BATCH_1, QC_LABELS_BATCH_2)
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


def run_annotation_blood(
    project_path,
    batch,
    n_top_genes,
    method_hvg,
    integration,
):
    log.info(f"Read adata and metadata for blood - batch {batch}")
    adata_name = f"blood_{n_top_genes}{method_hvg}_{integration}_clustered"
    adata = sc.read(project_path / f"{adata_name}.h5ad")
    meta_df = pd.read_csv(
        project_path / f"Batch{batch}_celllevel_metadata.tsv", sep="\t"
    )
    meta_df["donor"] = meta_df["donor"].apply(
        lambda x: ("-").join((x.split("-")[0], x.split("-")[1][1:]))
    )
    age_dict = meta_df.groupby("donor")["age"].agg(lambda x: x.unique()[0]).to_dict()
    sex_dict = meta_df.groupby("donor")["sex"].agg(lambda x: x.unique()[0]).to_dict()
    diag_dict = (
        meta_df.groupby("donor")["diagnosis"].agg(lambda x: x.unique()[0]).to_dict()
    )

    adata.obs["age"] = adata.obs["pt"].map(age_dict).astype(int)
    adata.obs["sex"] = adata.obs["pt"].map(sex_dict)
    adata.obs["diagnosis"] = adata.obs["pt"].map(diag_dict)

    log.critical(f"CellTypist annotation")
    adata_celltypist = adata.copy()
    adata_celltypist.X = adata.layers["counts"]
    sc.pp.normalize_per_cell(adata_celltypist, counts_per_cell_after=10**4)
    sc.pp.log1p(adata_celltypist)
    adata_celltypist.X = adata_celltypist.X.toarray()

    models.download_models(
        force_update=True, model=["Immune_All_Low.pkl", "Immune_All_High.pkl"]
    )

    model_low = models.Model.load(model="Immune_All_Low.pkl")
    model_high = models.Model.load(model="Immune_All_High.pkl")

    predictions_high = celltypist.annotate(
        adata_celltypist, model=model_high, majority_voting=True
    )
    predictions_high_adata = predictions_high.to_adata()

    adata.obs["celltypist_cell_label_coarse"] = predictions_high_adata.obs.loc[
        adata.obs.index, "majority_voting"
    ]
    adata.obs["celltypist_conf_score_coarse"] = predictions_high_adata.obs.loc[
        adata.obs.index, "conf_score"
    ]

    predictions_low = celltypist.annotate(
        adata_celltypist, model=model_low, majority_voting=True
    )
    predictions_low_adata = predictions_low.to_adata()

    adata.obs["celltypist_cell_label_fine"] = predictions_low_adata.obs.loc[
        adata.obs.index, "majority_voting"
    ]
    adata.obs["celltypist_conf_score_fine"] = predictions_low_adata.obs.loc[
        adata.obs.index, "conf_score"
    ]

    log.critical(f"ScArches annotation on full data")
    adata_all = sc.read(project_path / "blood_filtered_normalized.h5ad")

    adata_to_map = adata_all.copy()
    for layer in list(adata_to_map.layers.keys()):
        if layer != "counts":
            del adata_to_map.layers[layer]
    adata_to_map.X = adata_to_map.layers["counts"]
    reference_model_features = pd.read_csv(
        "https://figshare.com/ndownloader/files/41436645", index_col=0
    )
    adata_to_map.var["gene_names"] = adata_to_map.var.index
    log.info(
        f"Total number of genes needed for mapping: {reference_model_features.shape[0]}"
    )

    log.info(
        f"Number of genes found in query dataset:{adata_to_map.var.index.isin(reference_model_features.index).sum()}"
    )
    log.info("Add missing genes")

    missing_genes = [
        gene_id
        for gene_id in reference_model_features.index
        if gene_id not in adata_to_map.var.index
    ]
    missing_gene_adata = sc.AnnData(
        X=csr_matrix(
            np.zeros(shape=(adata.n_obs, len(missing_genes))), dtype="float32"
        ),
        obs=adata.obs.iloc[:, :1],
        var=reference_model_features.loc[missing_genes, :],
    )
    missing_gene_adata.layers["counts"] = missing_gene_adata.X

    if "PCs" in adata_to_map.varm.keys():
        del adata_to_map.varm["PCs"]

    adata_to_map_augmented = sc.concat(
        [adata_to_map, missing_gene_adata],
        axis=1,
        join="outer",
        index_unique=None,
        merge="unique",
    )
    adata_to_map_augmented = adata_to_map_augmented[
        :, reference_model_features.index
    ].copy()

    assert (adata_to_map_augmented.var.index == reference_model_features.index).all()
    adata_to_map_augmented.var["gene_ids"] = adata_to_map_augmented.var.index
    adata_to_map_augmented.var.set_index("gene_names", inplace=True)

    log.info("Load trained model")
    if not os.path.exists("./reference_model"):
        os.mkdir("./reference_model")
    elif not os.path.exists("./reference_model/model.pt"):
        urllib.request.urlretrieve(
            "https://figshare.com/ndownloader/files/41436648",
            filename="reference_model/model.pt",
        )
    adata_to_map_augmented.obs["batch"] = adata_to_map_augmented.obs["pt"]
    scarches_model = sca.models.SCVI.load_query_data(
        adata=adata_to_map_augmented,
        reference_model="./reference_model",
        freeze_dropout=True,
    )
    scarches_model.train(max_epochs=500, plan_kwargs=dict(weight_decay=0.0))
    adata.obsm["X_scVI"] = scarches_model.get_latent_representation()
    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.umap(adata)

    ref_emb = sc.read(
        filename="reference_embedding.h5ad",
        backup_url="https://figshare.com/ndownloader/files/41376264",
    )
    ref_emb.obs["reference_or_query"] = "reference"

    adata_emb = sc.AnnData(X=adata.obsm["X_scVI"], obs=adata.obs)
    adata_emb.obs["reference_or_query"] = "query"

    emb_ref_query = sc.concat(
        [ref_emb, adata_emb], axis=0, join="outer", index_unique=None, merge="unique"
    )
