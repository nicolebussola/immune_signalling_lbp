# Immune signaling links blood to brain in living humans

Reproducibility code for:

> **Immune signaling links blood to brain in living humans**
> Nicole Bussola, Lora E. Liharska, John F. Fullard, Ryan C. Thompson, Eric Vornholt, Girish N. Nadkarni, Eric E. Schadt, Noam D. Beckmann, Vibhuti Patel, Ombretta Melaiu, Panos Roussos, Brian H. Kopell†, Alexander W. Charney†

---

## Overview

This repository contains all custom code used to reproduce the analyses described in the paper. Data from 259 blood-brain pairs from 169 Living Brain Project (LBP) participants were analyzed using single-cell and bulk transcriptomics to characterize how peripheral immune cells communicate with the living human brain.

The paper makes four main contributions:
1. A scalable paired blood-brain sampling strategy in a neurosurgical setting.
2. A publicly available transcriptomic resource from living humans (available on Synapse).
3. Characterization of peripheral blood immune cell (PBIC) transcriptional states in the brain.
4. A multi-step neuroimmune communication pipeline.

---

## Study cohorts

| Cohort | Pairs | Participants | Tissue | Technology |
|--------|-------|--------------|--------|------------|
| Cohort 1 | 13 | 9 | Immune-sorted PFC (scRNA-seq) + blood PBMCs | 10x 3′ v2/v3 |
| Cohort 2 | 13 | 7 | Unsorted PFC nuclei (snRNA-seq) + blood PBMCs | 10x 3′ v3.1 |
| Cohort 3 | 233 | 153 | Bulk RNA-seq (blood + brain) | [Liharska et al.](https://www.nature.com/articles/s41380-025-03163-1) |

---

## Repository structure

```
immune_signalling_lbp/
├── preprocessing/          # Raw data processing: QC, filtering, normalization, annotation
├── PBICs_analysis/         # DE analysis, TF/pathway activity, compositional analysis, cross-dataset comparison
├── tensorS2R/              # Tensor-S2R CCC pipeline (LIANA + Tensor-cell2cell + LR filtering + enrichment + downstream)
├── LBP_env.yml             # Conda environment
└── pyproject.toml
```

---

## Analyses

### 1. Preprocessing (`preprocessing/`)

Processes raw 10x CellBender-filtered outputs into merged, annotated `.h5ad` objects (human cohorts), and preprocesses the public mouse scRNA-seq dataset (GEO GSE225948) for cross-species comparison.

**Human (Cohorts 1 & 2):** ambient RNA removal (CellBender) → per-sample QC + doublet detection → merge → normalization (log1p + scran) → HVG selection → batch correction (Harmony / scANVI) → cell type annotation (blood: CellTypist + scArches; brain: UniCell Deconvolve).

**Mouse:** `process_mouse_data.py` — loads GEO GSE225948 (22 paired brain + blood samples, 11 mice), assigns consistent sample keys by (treatment, Replicate) pairing, log-normalizes, filters cell types with <100 cells, and computes UMAP following the original publication's pipeline.

See [`preprocessing/README.md`](preprocessing/README.md) for full details.

---

### 2. PBIC analysis (`PBICs_analysis/`)

Analyzes the six peripheral blood immune cell types (CD4+ T, CD8+ T, CD14+ Monocytes, CD16+ Monocytes, NK, B cells) shared between blood and brain (PBICs).

| Script | Description |
|--------|-------------|
| `DEA_dreamlet.R` | Pseudo-bulk DE analysis between blood and brain per cell type (dreamlet, human + mouse) |
| `compositional_analysis.R` | PBIC abundance differences between tissues (sccomp) |
| `TF_and_pathway_activities.py` | TF activity (CollecTRI) and pathway activity (PROGENy) from DE signatures |
| `signature_comparison.r` | Cross-dataset concordance: Cohort 2, mouse, bulk Cohort 3, external data |

See [`PBICs_analysis/README.md`](PBICs_analysis/README.md) for execution order and input/output details.

---

### 3. Tensor-S2R communication pipeline (`tensorS2R/`)

End-to-end ligand-receptor communication pipeline structured as three sequential steps:

| Step | Script | Description |
|------|--------|-------------|
| 1 | `step1_factorize.py` | LIANA+ rank-aggregate by sample → Tensor-cell2cell factorization |
| 2 | `step2_enrich.py` | LR pair filtering (top 70th-percentile loadings; scran DE) → Enrichr pathway enrichment (WikiPathways, Reactome, GO) → pseudo-bulk LR product heatmaps |
| 3 | `step3_downstream.py` | Pseudo-bulk ORA via decoupler on MSigDB gene sets at cell and pseudo-bulk level; optionally computes LR product heatmaps |

Step 2 dispatches automatically based on cohort and design. For Cohort 2, it additionally runs Neuron↔Microglia (brain_coarse) and Neuron↔Monocyte (brain_blood_coarse) enrichment analyses.

**Human pipeline** — analytical designs:

| Mode | Data | Key interactions |
|------|------|----------------|
| `brain_coarse` | Brain only | Monocytes ↔ Microglia, Neurons ↔ Microglia |
| `micro_blood_coarse` | Blood PBMCs + microglia (brain) | Monocytes ↔ Microglia |
| `brain_blood_coarse` | Full blood + all brain cell types | All PBICs ↔ all brain types |

```bash
python -m tensorS2R \
    -p <project_path> \
    -b <cohort_1|cohort_2> \
    -f <brain_coarse|micro_blood_coarse|brain_blood_coarse> \
    [--step <1|2|3|all>] \
    [-s <sample_key>] \
    [-g <groupby>] \
    [-m <inner|outer>] \
    [-r <rank>] \
    [-d <cuda|cpu>] \
    [-t <ora_threshold>]
```

**Mouse pipeline** (`factorization_run_mouse.py`) — LIANA mouse consensus + Tensor-cell2cell on the GEO GSE225948 paired brain/blood dataset (sham, D02, D14 conditions; rank fixed at 20):

```bash
python -m tensorS2R.factorization_run_mouse \
    -p <mouse_processed_path> \
    [-s <sample_key>] \
    [-g <groupby>] \
    [-f <factorization_type>] \
    [-m <inner|outer>] \
    [-r <rank>] \
    [-d <cuda|cpu>]
```

See [`tensorS2R/README.md`](tensorS2R/README.md) for full argument reference and output structure.

---

## Environment

A conda environment file is provided:

```bash
conda env create -f LBP_env.yml
conda activate LBP
```

Key Python dependencies: `scanpy`, `liana` (v1.1.0), `cell2cell` (v0.7.3), `decoupler`, `gseapy` (v1.1.3), `tensorly`, `bokeh`.

Key R dependencies: `dreamlet`, `sccomp`, `decoupleR`, `zellkonverter`.

---

## Data availability

Processed `.h5ad` objects and raw counts are deposited on Synapse. [TODO: add Synapse ID]

---

## Citation

[TODO: add citation once published]

---

## License

[TODO]
