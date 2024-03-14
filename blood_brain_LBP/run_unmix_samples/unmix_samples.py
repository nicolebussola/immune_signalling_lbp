import logging

import anndata as ad
import numpy as np
import scanpy as sc
from rich.logging import RichHandler
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, matthews_corrcoef
from sklearn.model_selection import StratifiedKFold
from utils import gridlayout

logging.getLogger().setLevel(logging.INFO)
FORMAT = "%(message)s"
logging.basicConfig(
    level="INFO", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
)

log = logging.getLogger("rich")


def run_unmix_samples(
    project_path, blood_adata_name_mixed, blood_adata_name, output_path_plot
):
    np.random.seed(42)

    blood_adata = sc.read(project_path / blood_adata_name_mixed)

    adata_R = blood_adata[
        (blood_adata.obs["pt"].isin(["280", "289"])) & (blood_adata.obs["side"] == "R")
    ]
    adata_L = blood_adata[
        (blood_adata.obs["pt"].isin(["280", "289"])) & (blood_adata.obs["side"] == "L")
    ]
    log.info("adata - R", adata_R.shape)
    log.info("adata - L", adata_L.shape)

    sc.pp.normalize_total(adata_R)
    sc.pp.log1p(adata_R)

    sc.pp.normalize_total(adata_L)
    sc.pp.log1p(adata_L)

    log.info("Train RF model - Cross validation")
    X = adata_R.X.copy()
    y = np.array(adata_R.obs["pt"])

    clf = RandomForestClassifier(n_estimators=500, random_state=42)

    num_folds = 5

    for fold, (train_idx, test_idx) in enumerate(
        StratifiedKFold(n_splits=num_folds, random_state=42, shuffle=True).split(X, y)
    ):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        clf.fit(X_train, y_train)
        y_pred = clf.predict(X_test)

        log.info(f"Fold {fold + 1} MCC:\n", matthews_corrcoef(y_test, y_pred))

    log.info("Train RF model on R data and test on L")
    X_train = adata_R.X.copy()
    y_train = np.array(adata_R.obs["pt"])

    X_test = adata_L.X.copy()
    y_test = np.array(adata_L.obs["pt"])

    clf = RandomForestClassifier(n_estimators=500, random_state=42)
    clf.fit(X_train, y_train)

    y_pred = clf.predict(X_test)

    log.info("We expect poor performace on test")
    accuracy = accuracy_score(y_test, y_pred)
    mcc = matthews_corrcoef(y_test, y_pred)
    log.info(f"Accuracy: {accuracy:.2f}")
    log.info(f"MCC: {mcc:.2f}")

    adata_L.obs["pt_rf"] = y_pred
    log.info("Visualized swapped labels")
    sc.pp.pca(adata_L, n_comps=50)
    sc.pp.neighbors(adata_L)
    sc.tl.umap(adata_L)

    adata_L.obs["changed"] = adata_L.obs.apply(
        lambda row: row["pt"] if row["pt"] == row["pt_rf"] else "swap", axis=1
    )

    adata_R.obs["changed"] = adata_R.obs["pt"]

    log.info("Visualized concatenated data")
    adata_blood_concat = ad.concat([adata_L, adata_R], join="outer")
    adata_blood_concat.obs
    sc.pp.normalize_total(adata_blood_concat)
    sc.pp.log1p(adata_blood_concat)
    sc.pp.pca(adata_blood_concat, n_comps=50)
    sc.pp.neighbors(adata_blood_concat)
    sc.tl.umap(adata_blood_concat)

    sc.pl.umap(
        adata_blood_concat,
        color=["pt", "pt_rf", "changed", "side"],
    )

    gridlayout(
        ["pt", "side", "pt_rf", "changed"],
        adata_blood_concat,
        width=900,
        ncols=1,
        fname=output_path_plot / f"{blood_adata_name_mixed}_RF.html",
    )

    adata_L.obs[(adata_L.obs["pt"] == "280") & (adata_L.obs["changed"] == "swap")]
    adata_L.obs[(adata_L.obs["pt"] == "289") & (adata_L.obs["changed"] == "swap")].index

    log.info("Write .txt files for variant calling")
    L280_to_L289 = [
        ("-").join((cb.split("_")[0], "1"))
        for cb in list(
            adata_L.obs[
                (adata_L.obs["pt"] == "280") & (adata_L.obs["changed"] == "swap")
            ].index
        )
    ]
    L280_to_L280 = [
        ("-").join((cb.split("_")[0], "1"))
        for cb in list(
            adata_L.obs[
                (adata_L.obs["pt"] == "280") & (adata_L.obs["changed"] == "280")
            ].index
        )
    ]

    L289_to_L280 = [
        ("-").join((cb.split("_")[0], "1"))
        for cb in list(
            adata_L.obs[
                (adata_L.obs["pt"] == "289") & (adata_L.obs["changed"] == "swap")
            ].index
        )
    ]
    L289_to_L289 = [
        ("-").join((cb.split("_")[0], "1"))
        for cb in list(
            adata_L.obs[
                (adata_L.obs["pt"] == "289") & (adata_L.obs["changed"] == "289")
            ].index
        )
    ]

    for cell_set in [L280_to_L289, L280_to_L280, L289_to_L280, L289_to_L289]:
        file_name = project_path / f"{cell_set}.txt"
        with open(file_name, "w") as file:
            for item in cell_set:
                file.write("CB:Z:%s\n" % item)

    log.info("Save blood data with new labels")
    blood_adata_unmixed = blood_adata.copy()
    blood_adata_unmixed.obs.loc[
        blood_adata_unmixed.obs.index.isin(adata_L.obs.index), "pt"
    ] = adata_L.obs["pt_rf"]
    sc.write(project_path / blood_adata_name, blood_adata)
