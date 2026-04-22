# Preprocessing

This folder contains the preprocessing code used to generate the processed `.h5ad` inputs.

## Structure

- `run_ddqc/`: compute ddqc metrics for each sample from CellBender-filtered 10x `.h5` inputs.
- `run_preprocessing/`: sample-level preprocessing, QC metric calculation, doublet detection, feature selection, and QC plotting.
- `run_merge/`: merge per-sample processed `.h5ad` files, filter cells, normalize counts, and select HVGs.
- `run_annotation/`: annotate merged data using CellTypist and scArches-based reference mapping.
- `utils.py`, `labels.py`, `immune_markers_blood.py`, `sc-type.R`: shared utilities and annotation helpers.

## Intended flow

The preprocessing is not a single script. The intended order is:

1. Run `run_ddqc` on raw sample-level CellBender-filtered inputs.
2. Run `run_preprocessing` on the same sample directories together with the ddqc outputs.
3. Run `run_merge` on the resulting per-sample `.h5ad` files.
4. Run annotation on the merged data.

The outputs from this folder are the processed `.h5ad` objects.

## CLI entrypoints

From the repository root, the current entrypoints are:

```bash
python -m preprocessing.run_ddqc -i /path/to/batch_root -t blood
python -m preprocessing.run_preprocessing -i /path/to/batch_root -o /path/to/plots -t blood
python -m preprocessing.run_merge -i /path/to/preprocessed_root -b 1 -o /path/to/output -t blood
```


## Input assumptions

### `run_ddqc`

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

### `run_preprocessing`

Expected per-sample files:

- `PT-XXX-<side>-CellBender_filtered.h5`
- `PT-XXX-<side>_CellBender_filtered_ddqc.csv`

Expected metadata/columns generated or used later include:

- ddqc output indexed by `barcodekey`
- sample naming compatible with `PT-XXX-<tissue>-<side>`

Outputs:

- QC plots
- in-memory preprocessing state used to prepare per-sample `.h5ad`



### `run_merge`

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
- feature selection using one of:
  - `HighlyDeviant`
  - `cell_ranger`
  - `HighlyDeviant_cr`




