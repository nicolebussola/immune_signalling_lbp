from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report, matthews_corrcoef
from sklearn.model_selection import StratifiedKFold
from utils import gridlayout

np.random.seed(42)


project_path = Path(
    "/sc/arion/projects/psychgen/lbp/data/ScProcesses_brainBlood_nicole/batch_2/"
)

blood_adata = sc.read(
    project_path / "batch_2_blood_filtered_4000HighlyDeviant_20_mixed.h5ad"
)


L_280 = blood_adata[(blood_adata.obs["pt"] == "280") & (blood_adata.obs["side"] == "L")]
R_280 = blood_adata[(blood_adata.obs["pt"] == "280") & (blood_adata.obs["side"] == "R")]
L_289 = blood_adata[(blood_adata.obs["pt"] == "289") & (blood_adata.obs["side"] == "L")]
R_289 = blood_adata[(blood_adata.obs["pt"] == "289") & (blood_adata.obs["side"] == "R")]
print("280 - L", L_280.shape)
print("289 - L", L_289.shape)
print("280 - R", R_280.shape)
print("289 - R", R_289.shape)
print("280:", L_280.shape[0] + R_280.shape[0])
print("289:", L_289.shape[0] + R_289.shape[0])
print("Left:", L_280.shape[0] + L_289.shape[0])
print("Right:", R_280.shape[0] + R_289.shape[0])

adata_R = blood_adata[
    (blood_adata.obs["pt"].isin(["280", "289"])) & (blood_adata.obs["side"] == "R")
]
adata_L = blood_adata[
    (blood_adata.obs["pt"].isin(["280", "289"])) & (blood_adata.obs["side"] == "L")
]
adata_L, adata_R.shape


sc.pp.normalize_total(adata_R)
sc.pp.log1p(adata_R)

sc.pp.normalize_total(adata_L)
sc.pp.log1p(adata_L)
adata_R, adata_L


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

    report = classification_report(y_test, y_pred)
    print(f"Fold {fold + 1} MCC:\n", matthews_corrcoef(y_test, y_pred))


X_train = adata_R.X.copy()
y_train = np.array(adata_R.obs["pt"])

X_test = adata_L.X.copy()
y_test = np.array(adata_L.obs["pt"])

clf = RandomForestClassifier(n_estimators=500, random_state=42)
clf.fit(X_train, y_train)

y_pred = clf.predict(X_test)


accuracy = accuracy_score(y_test, y_pred)
mcc = matthews_corrcoef(y_test, y_pred)
print(f"Accuracy: {accuracy:.2f}")
print(f"MCC: {mcc:.2f}")


sc.pp.pca(adata_R, n_comps=50)
sc.pp.neighbors(adata_R)
sc.tl.umap(adata_R)

adata_R.obs["pt_rf"] = y_train


sc.pl.umap(
    adata_R,
    color=["pt", "pt_rf", "log1p_total_counts", "pct_counts_mt", "percent_ribo"],
)

# %%
adata_L.obs["pt_rf"] = y_pred

sc.pp.pca(adata_L, n_comps=50)
sc.pp.neighbors(adata_L)
sc.tl.umap(adata_L)

adata_L.obs["changed"] = adata_L.obs.apply(
    lambda row: row["pt"] if row["pt"] == row["pt_rf"] else "swap", axis=1
)


sc.pl.umap(
    adata_L,
    color=[
        "pt",
        "pt_rf",
        "changed",
        "log1p_total_counts",
        "pct_counts_mt",
        "percent_ribo",
    ],
)


gridlayout(["pt", "pt_rf", "changed"], adata_L, width=900, ncols=1)

adata_R.obs["changed"] = adata_R.obs["pt"]

# %%
adata_blood_concat = ad.concat([adata_L, adata_R], join="outer")
adata_blood_concat.obs
sc.pp.normalize_total(adata_blood_concat)
sc.pp.log1p(adata_blood_concat)
sc.pp.pca(adata_blood_concat, n_comps=50)
sc.pp.neighbors(adata_blood_concat)
sc.tl.umap(adata_blood_concat)

# %%
adata_blood_concat.obs[["side", "pt", "pt_rf", "changed"]]

# %%
sc.pl.umap(
    adata_blood_concat,
    color=[
        "pt",
        "pt_rf",
        "changed",
        "side",
        "log1p_total_counts",
        "pct_counts_mt",
        "percent_ribo",
    ],
)


adata_L.obs[(adata_L.obs["pt"] == "280") & (adata_L.obs["changed"] == "swap")]
adata_L.obs[(adata_L.obs["pt"] == "289") & (adata_L.obs["changed"] == "swap")].index

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


file_name = project_path / "L289_to_L280_new.txt"
with open(file_name, "w") as file:
    for item in L289_to_L280:
        file.write("CB:Z:%s\n" % item)


blood_adata_unmixed = blood_adata.copy()
blood_adata_unmixed.obs.loc[
    blood_adata_unmixed.obs.index.isin(adata_L.obs.index), "pt"
] = adata_L.obs["pt_rf"]
sc.write(project_path / "batch_2_blood_filtered_4000HighlyDeviant_20.h5ad", blood_adata)
