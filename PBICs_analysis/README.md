# PBICs Analysis

This folder contains four analysis scripts to reproduce results and figures for the PBICs (B cells, CD4+ T cells, CD8+ T cells, CD16+ natural killer (NK) cells, CD14+ monocytes, and CD16+ monocytes) analysis.

## Scripts

### 1. `compositional_analysis.R`

Standalone compositional analysis using *sccomp* [1].

Inputs:
- human paired blood-brain `.h5ad`

Expected human metadata used by the script:
- `cell_type.v2`
- `pt`
- `tissue`
- `side`
- `chemistry`

Outputs:
- `figures/supp_table_1_cell_counts.csv`
- `figures/sccomp_boxplot_tissue.pdf`



### 2. `DEA_dreamlet.R`

Differential expression and pseudobulk contrast generation for human and mouse using *dreamlet* [2].

Inputs:
- human paired blood-brain `.h5ad`
- mouse paired blood-brain `.h5ad`

Expected human metadata used by the script:
- `pt`
- `tissue`
- `chemistry`
- `cell_type.v2`

Expected mouse metadata used by the script:
- `cell_type.v2`
- `parent`
- `sub.celltype`
- `tissue`
- `treatment`
- `sample`

Primary outputs written to the working directory:
- human:
  - `df_cd14_dreamlet.csv`
  - `df_cd16_dreamlet.csv`
  - `df_cd4_dreamlet.csv`
  - `df_cd8_dreamlet.csv`
  - `df_nk_dreamlet.csv`
  - `df_bc_dreamlet.csv`
  - `mat_cd14_dreamlet.csv`
  - `mat_cd16_dreamlet.csv`
  - `mat_cd4_dreamlet.csv`
  - `mat_cd8_dreamlet.csv`
  - `mat_nk_dreamlet.csv`
  - `mat_bc_dreamlet.csv`
- mouse:
  - `df_cd14_dreamlet_mouse.csv`
  - `df_cd4_dreamlet_mouse.csv`
  - `df_cd8_dreamlet_mouse.csv`
  - `df_bc_dreamlet_mouse.csv`

Notes:
- `TF_and_pathway_activities.py` depends on the human `df_*_dreamlet.csv` and `mat_*_dreamlet.csv` files from this script.
- `signature_comparison.r` depends on both the human and mouse `df_*_dreamlet*.csv` files from this script.

### 3. `TF_and_pathway_activities.py`

Transcription factor (TF) [3] and pathway activity analysis [4] using DEA outputs from `DEA_dreamlet.R`.

Required inputs in the working directory:
- `df_cd14_dreamlet.csv`
- `df_cd16_dreamlet.csv`
- `df_cd4_dreamlet.csv`
- `df_cd8_dreamlet.csv`
- `df_nk_dreamlet.csv`
- `df_bc_dreamlet.csv`
- `mat_cd14_dreamlet.csv`
- `mat_cd16_dreamlet.csv`
- `mat_cd4_dreamlet.csv`
- `mat_cd8_dreamlet.csv`
- `mat_nk_dreamlet.csv`
- `mat_bc_dreamlet.csv`

Outputs:
- `figures/CD14_volcano.pdf`
- `figures/collectri_heatmap_bq.png`
- `figures/progeny_heatmap_human_subset.png`


### 4. `signature_comparison.r`

Cross-dataset signature comparison using dreamlet outputs from `DEA_dreamlet.R`, bulk data, and external datasets [5,6,7].

Required local dreamlet inputs in the working directory:
- human:
  - `df_cd14_dreamlet.csv`
  - `df_cd16_dreamlet.csv`
  - `df_cd4_dreamlet.csv`
  - `df_cd8_dreamlet.csv`
  - `df_nk_dreamlet.csv`
  - `df_bc_dreamlet.csv`
- mouse:
  - `df_cd14_dreamlet_mouse.csv`
  - `df_cd4_dreamlet_mouse.csv`
  - `df_cd8_dreamlet_mouse.csv`
  - `df_bc_dreamlet_mouse.csv`

Additional external inputs:
- bulk DE `.RDS` with at least:
  - `symbol`
  - `logFC`
- CSF Excel file
- replication Excel file
- cohort 2 signature directory containing `de_<braincell>_<bloodcell>.csv`


## Recommended execution order

If running the full dreamlet-based analysis:

1. Run `DEA_dreamlet.R`
2. Run `TF_and_pathway_activities.py`
3. Run `signature_comparison.r`

If running the compositional analysis:

1. Run `compositional_analysis.R`

## Reproducibility caveats

- The annotation mapping approach in the R scripts depends on Bioconductor annotation databases and may vary if package versions differ.


## References
[1] Mangiola, S. et al. sccomp: Robust differential composition and variability anal-
ysis for single-cell data. Proceedings of the National Academy of Sciences 120,
e2203828120 (2023)

[2] Hoffman, G. E. et al. Efficient differential expression analysis of large-scale single
cell transcriptomics data using dreamlet. bioRxiv (2023).

[3] M¨uller-Dott, S. et al. Expanding the coverage of regulons from high-confidence
prior knowledge for accurate estimation of transcription factor activities. Nucleic
acids research 51, 10934–10949 (2023)

[4] Holland, C. H. et al. Robustness and applicability of transcription factor and
pathway analysis tools on single-cell rna-seq data. Genome biology 21, 1–19
(2020)

[5] Garcia-Bonilla, L. et al. Analysis of brain and blood single-cell transcriptomics
in acute and subacute phases after experimental stroke. Nature Immunology 25,
357–370 (2024)

[6] Pasciuto, E. et al. Microglia require cd4 t cells to complete the fetal-to-adult
transition. Cell 182, 625–640 (2020)

[7] Schafflick, D. et al. Integrated single cell analysis of blood and cerebrospinal
fluid leukocytes in multiple sclerosis. Nature communications 11, 247 (2020)
