import os
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc
import scanpy.external as sce


def pipeline_bonilla(adata):
    """Reproduce original publication UMAP: log1p_norm → HVG 3k → scale → PCA 40 → Harmony → UMAP."""
    adata = adata.copy()
    adata.X = adata.layers["log1p_norm"].copy()
    sc.pp.highly_variable_genes(adata, n_top_genes=3000)
    adata = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata)
    sc.tl.pca(adata, n_comps=40)
    adata.obs["Replicate"] = adata.obs["Replicate"].astype(str)
    sce.pp.harmony_integrate(adata, "Replicate", adjusted_basis="X_pca_harmony")
    sc.pp.neighbors(adata, use_rep="X_pca_harmony")
    sc.tl.umap(adata)
    return adata


def subpipeline_bonilla(adata, brain_blood=False):
    """
    Re-derive UMAP from raw counts: remove mt/ribo/Gm genes, normalise to 10k,
    HVG 2k → scale → PCA 15 → Harmony [→ tissue correction] → UMAP.
    """
    adata = adata.copy()
    adata.layers["counts"] = adata.X.copy()

    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
    adata.var["gm"] = adata.var_names.str.contains("Gm")
    rem = adata.var["mt"] | adata.var["ribo"] | adata.var["gm"]
    adata = adata[:, ~rem].copy()

    scales_counts = sc.pp.normalize_total(adata, target_sum=1e4, inplace=False)
    adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)
    adata.X = adata.layers["log1p_norm"].copy()
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata)
    sc.tl.pca(adata, n_comps=15)
    adata.obs["Replicate"] = adata.obs["Replicate"].astype(str)

    if brain_blood:
        sce.pp.harmony_integrate(adata, "Replicate", adjusted_basis="X_pca_harmony")
        sce.pp.harmony_integrate(
            adata, "tissue", basis="X_pca_harmony", adjusted_basis="X_pca_harmony_tissue"
        )
        sc.pp.neighbors(adata, use_rep="X_pca_harmony_tissue")
    else:
        sce.pp.harmony_integrate(adata, "Replicate", adjusted_basis="X_pca_harmony")
        sc.pp.neighbors(adata, use_rep="X_pca_harmony")

    sc.tl.umap(adata, min_dist=0.1)
    return adata


def filter_low_cells(adata, cell_col="parent", min_cells=100):
    """Remove cell types with fewer than min_cells cells."""
    counts = adata.obs[cell_col].value_counts()
    keep = counts[counts > min_cells].index
    return adata[adata.obs[cell_col].isin(keep), :].copy()


def process_and_filter_adata(adata, cell_col="parent", min_cells=100):
    """Store raw counts, log-normalise to 10k, filter low-count cell types."""
    adata.layers["counts"] = adata.X.copy()
    scales_counts = sc.pp.normalize_total(adata, target_sum=1e4, inplace=False)
    adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)
    return filter_low_cells(adata, cell_col, min_cells)


def build_sample_key(adata_br, adata_bl):
    """
    Map (treatment, Replicate) → mouse_i using a key consistent across tissues,
    so paired brain and blood samples from the same mouse share the same label.
    """
    combined = pd.concat([adata_br.obs, adata_bl.obs])
    pair_keys = (
        combined[["treatment", "Replicate"]]
        .drop_duplicates()
        .sort_values(["treatment", "Replicate"])
        .reset_index(drop=True)
    )
    return {
        (row["treatment"], str(row["Replicate"])): f"mouse_{i + 1}"
        for i, row in pair_keys.iterrows()
    }


project_path = Path("/path/to/downloaded/data/")
data_path = project_path / "data" / "GSE225948_RAW"
out_path = project_path / "data" / "mouse_processed"
out_path.mkdir(parents=True, exist_ok=True)


young_list = [
    f for f in sorted(os.listdir(data_path)) if "aged" not in f and f.endswith(".gz")
]

anndata_list = []
for i, cnt_meta in enumerate([young_list[j : j + 2] for j in range(0, len(young_list), 2)]):
    counts = pd.read_csv(data_path / cnt_meta[0], index_col=0)
    metadata = pd.read_csv(data_path / cnt_meta[1], index_col=0)
    adata_i = sc.AnnData(X=counts.T)
    adata_i.obs = metadata
    anndata_list.append(adata_i)
    print(i, cnt_meta, adata_i.shape)

# 22 samples: first 11 are brain, last 11 are blood
adata_br_raw = ad.concat(anndata_list[:11])
adata_bl_raw = ad.concat(anndata_list[11:])  

pair_to_sample = build_sample_key(adata_br_raw, adata_bl_raw)

adata = ad.concat([adata_br_raw, adata_bl_raw])
adata.obs["sample"] = [
    pair_to_sample[(t, str(r))]
    for t, r in zip(adata.obs["treatment"], adata.obs["Replicate"])
]


adata = process_and_filter_adata(adata)

adata.obs["cell_type.v2"] = (
    adata.obs["parent"].astype(str) + "_" + adata.obs["tissue"].astype(str)
)

adata_umap = pipeline_bonilla(adata)
adata.obsm = adata_umap.obsm
adata.obsp = adata_umap.obsp

sc.write(out_path / "paired_blood_brain_mouse.h5ad", adata)


adata_br_filt = adata_br_raw[adata_br_raw.obs.index.isin(adata.obs.index)].copy()
adata_br_umap = subpipeline_bonilla(adata_br_filt)

sc.settings.figdir = str(out_path)
sc.pl.umap(
    adata_br_umap,
    color=["sub.celltype"],
    wspace=1,
    palette="tab10",
    save="_mouse_brain.pdf",
    show=False,
)
