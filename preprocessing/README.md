# Preprocessing

This folder contains the preprocessing code used to generate the processed `.h5ad` inputs for human (Cohorts 1 & 2) and mouse data.

---

## Human preprocessing (Cohorts 1 & 2)

### Structure

- `run_ddqc/`: compute ddqc metrics for each sample from CellBender-filtered 10x `.h5` inputs.
- `run_preprocessing/`: sample-level preprocessing, QC metric calculation, doublet detection, feature selection, and QC plotting.
- `run_merge/`: merge per-sample processed `.h5ad` files, filter cells, normalize counts, and select HVGs.
- `run_annotation/run_annotate.py`: annotate blood merged data using CellTypist and scArches-based reference mapping.
- `run_annotation/unicell_annotation.py`: annotate brain merged data using UniCell Deconvolve (v0.1.0).
- `utils.py`, `labels.py`: shared utilities and annotation helpers.

### Intended flow

1. Run `run_ddqc` on raw sample-level CellBender-filtered inputs.
2. Run `run_preprocessing` on the same sample directories together with the ddqc outputs.
3. Run `run_merge` on the resulting per-sample `.h5ad` files.
4. Annotate merged data:
   - **Blood**: `run_annotation/run_annotate.py` (CellTypist + scArches)
   - **Brain**: `run_annotation/unicell_annotation.py` (UniCell Deconvolve)

### CLI entrypoints

From the repository root:

```bash
python -m preprocessing.run_ddqc -i /path/to/batch_root -t blood
python -m preprocessing.run_preprocessing -i /path/to/batch_root -o /path/to/plots -t blood
python -m preprocessing.run_merge -i /path/to/preprocessed_root -b 1 -o /path/to/output -t blood
```

### Input assumptions

#### `run_ddqc`

Expected directory layout:

```text
<batch_root>/
  blood/
    PT-XXX-blood-L/
      PT-XXX-L-CellBender_filtered.h5
    PT-XXX-blood-R/
      PT-XXX-R-CellBender_filtered.h5
  brain/
    ...
```

Outputs:

- `<sample_dir>/<PT>-<side>_CellBender_filtered_ddqc.csv`

#### `run_preprocessing`

Expected per-sample files:

- `PT-XXX-<side>-CellBender_filtered.h5`
- `PT-XXX-<side>_CellBender_filtered_ddqc.csv`

Expected metadata/columns generated or used later include:

- ddqc output indexed by `barcodekey`
- sample naming compatible with `PT-XXX-<tissue>-<side>`

Outputs:

- QC plots
- in-memory preprocessing state used to prepare per-sample `.h5ad`

#### `run_merge`

Expected inputs:

- one `.h5ad` per sample inside `<input_path>/<tissue>/`

Expected sample filename convention:

- filenames are parsed to recover `pt`, `side`, and `tissue`

Main processing steps:

- merge all samples
- remove consensus doublets
- keep `passed_qc == "True"`
- mitochondrial filtering
- log1p normalization
- scran size-factor normalization
- feature selection using one of: `HighlyDeviant`, `cell_ranger`, `HighlyDeviant_cr`

---

## Mouse preprocessing

### `process_mouse_data.py`

Preprocesses the publicly available paired brain+blood mouse scRNA-seq dataset (GEO GSE225948; Garcia-Bonilla et al. 2024, *Nature Immunology*) for cross-species comparison.

**Dataset:** 22 paired brain and blood samples from 11 male wild-type mice (8–12 weeks old); conditions: sham surgery, ischemic stroke at Day 2 (D02, acute) and Day 14 (D14, recovery).

**Processing steps:**

1. Load raw count + metadata `.gz` files from `GSE225948_RAW/` (young mice only).
2. Concatenate brain (first 11 samples) and blood (last 11 samples) separately.
3. Assign consistent `sample` labels using `(treatment, Replicate)` pairing so brain and blood from the same mouse share the same ID.
4. Log-normalize to 10k counts per cell; store raw counts in `layers["counts"]`.
5. Filter cell types with fewer than 100 cells (`parent` column).
6. Add `cell_type.v2 = parent_tissue` composite label.
7. Compute UMAP following the original publication pipeline (`pipeline_bonilla`: HVG 3k → scale → PCA 40 → Harmony on Replicate).
8. Write combined `paired_blood_brain_mouse.h5ad` to `mouse_processed/`.
9. Compute separate brain-only UMAP (`subpipeline_bonilla`: mt/ribo/Gm removal → HVG 2k → PCA 15 → Harmony).

**Key metadata columns in output `adata.obs`:**

| Column | Description |
|---|---|
| `sample` | Mouse ID (`mouse_1` … `mouse_11`), consistent across brain and blood |
| `tissue` | Tissue of origin |
| `treatment` | `Sham`, `D02`, or `D14` |
| `Replicate` | Mouse replicate number |
| `parent` | Coarse cell type (e.g. `Mo+`, `Mo-`, `MdC`, `Bc`, `Tc`) |
| `sub.celltype` | Fine-grained cell type (e.g. `Mo1`, `MdC2`) |
| `cell_type.v2` | `parent_tissue` composite label |

**Usage:**

```bash
python preprocessing/process_mouse_data.py
```

Edit the `project_path` variable at the top of the script to point to your local data directory before running. Raw GEO files should be placed in `<project_path>/data/GSE225948_RAW/`. Output is written to `<project_path>/data/mouse_processed/`.

**Note:** Mouse–human ortholog mapping for cross-species DE comparison is performed separately in R using `biomaRt::getLDS` and is a prerequisite for running `PBICs_analysis/DEA_dreamlet.R` and `PBICs_analysis/signature_comparison.r` on mouse data.

---

## Annotation scripts

### `run_annotation/run_annotate.py` — Blood annotation (CellTypist + scArches)

Annotates blood merged data using CellTypist (Immune_All_Low and Immune_All_High models) and scArches reference mapping against a pre-trained PBMC model from Figshare.

**Input:** merged blood `.h5ad` with `counts` layer (output of `run_merge`).

**Output:** `adata.obs` updated with `celltypist_cell_label_coarse`, `celltypist_cell_label_fine`, and scVI latent embedding in `adata.obsm["X_scVI"]`.

---

### `run_annotation/unicell_annotation.py` — Brain annotation (UniCell Deconvolve)

Annotates brain merged data using UniCell Deconvolve (v0.1.0) via its cloud API.

**Requires:** a UniCell API token — set the `UNICELL_TOKEN` environment variable or edit the script directly.

**Input:** HVG-processed brain `.h5ad` (output of `run_merge` + batch correction), which must contain:
- `layers["counts"]` — raw counts
- `obsp["connectivities"]` — neighbor graph (used for Leiden clustering)

| Cohort | Input file |
|---|---|
| Cohort 1 | `cohort_1_brain_filtered_4000HighlyDeviant_20_harmony_regressedScaled.h5ad` |
| Cohort 2 | `cohort_2_brain_filtered_4000HighlyDeviant_20_harmony_scaled.h5ad` |

**Output:** `{COHORT}_brain_unicell_annotations.csv` — full `adata.obs` containing UniCell cell-type predictions and Leiden cluster labels.

**Usage:** set `COHORT` at the top of the script, then run:

```bash
UNICELL_TOKEN=<your_token> python preprocessing/run_annotation/unicell_annotation.py
```
