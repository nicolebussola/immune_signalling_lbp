import os
from pathlib import Path

import scanpy as sc
import ucdeconvolve as ucd

# Replace with your token or set the UNICELL_TOKEN environment variable
ucd.api.authenticate(os.environ.get("UNICELL_TOKEN", "UNICELL_TOKEN"))

BRAIN_HVG_FILES = {
    "cohort_1": "cohort_1_brain_filtered_4000HighlyDeviant_20_harmony.h5ad",
    "cohort_2": "cohort_2_brain_filtered_4000HighlyDeviant_20_harmony.h5ad",
}

COHORT = "cohort_1"  # set to "cohort_1" or "cohort_2"
project_path = Path(
    f"/path/to/ptoject/{COHORT}/"
)

adata = sc.read(project_path / BRAIN_HVG_FILES[COHORT])

sc.tl.leiden(adata, resolution=2)

adata.X = adata.layers["counts"].copy()
ucd.tl.base(adata)

ucd.pl.base_clustermap(
    adata,
    groupby="leiden",
    n_top_celltypes=200,
    figsize=(20, 20),
    annot_kws={"fontsize": 20},
)

adata.obs.to_csv(project_path / f"{COHORT}_brain_unicell_annotations.csv")
