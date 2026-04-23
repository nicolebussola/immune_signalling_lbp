# Immune signaling links blood to brain in living humans

Reproducibility code for:

> **Immune signaling links blood to brain in living humans**
> Nicole Bussola, Lora E. Liharska, John F. Fullard, Ryan C. Thompson, Eric Vornholt, Girish N. Nadkarni, Eric E. Schadt, Noam D. Beckmann, Ombretta Melaiu, Panos Roussos, Brian H. Kopell†, Alexander W. Charney†

---

## Overview

This repository contains all custom code used to reproduce the analyses described in the paper. Data from 276 blood-brain pairs from 169 Living Brain Project (LBP) participants were analyzed using single-cell and bulk transcriptomics to characterize how peripheral immune cells communicate with the living human brain.

The paper makes four main contributions:
1. A scalable paired blood-brain sampling strategy in a neurosurgical setting.
2. A publicly available transcriptomic resource from living humans (available on Synapse).
3. Characterization of peripheral blood immune cell (PBIC) transcriptional states in the brain.
4. A multi-step neuroimmune communication architecture linking monocytes, microglia, and neurons.

---

## Study cohorts

| Cohort | Pairs | Participants | Tissue | Technology |
|--------|-------|--------------|--------|------------|
| Cohort 1 | 13 | 9 | Immune-sorted PFC (scRNA-seq) + blood PBMCs | 10x 3′ v2/v3 |
| Cohort 2 | 13 | 7 | Unsorted PFC nuclei (snRNA-seq) + blood PBMCs | 10x 3′ v3.1 |
| Cohort 3 | 233 | 153 | Bulk RNA-seq (blood + brain) | — |

---

## Repository structure

```
blood-brain-LBP/
├── preprocessing/          # Raw data processing: QC, filtering, normalization, annotation
├── PBICs_analysis/         # DE analysis, TF/pathway activity, compositional analysis, cross-dataset comparison
├── tensorS2R/              # Tensor-S2R CCC pipeline (LIANA + Tensor-cell2cell + enrichment)
├── notebooks/              # Exploratory notebooks
├── LBP_env.yml             # Conda environment
└── pyproject.toml
```

---

## Analyses

### 1. Preprocessing (`preprocessing/`)

Processes raw 10x CellBender-filtered outputs into merged, annotated `.h5ad` objects.

Steps: ambient RNA removal (CellBender) → per-sample QC + doublet detection → merge → normalization (log1p + scran) → HVG selection → batch correction (Harmony / scANVI) → cell type annotation.

See [`preprocessing/README.md`](preprocessing/README.md) for full details.

---

### 2. PBIC analysis (`PBICs_analysis/`)

Analyzes the six peripheral blood immune cell types (CD4+ T, CD8+ T, CD14+ Mono, CD16+ Mono, NK, B cells) shared between blood and brain.

| Script | Description |
|--------|-------------|
| `DEA_dreamlet.R` | Pseudo-bulk DE analysis between blood and brain per cell type (dreamlet, human + mouse) |
| `compositional_analysis.R` | PBIC abundance differences between tissues (sccomp) |
| `TF_and_pathway_activities.py` | TF activity (CollecTRI via decoupleR) and pathway activity (PROGENy) from DE signatures |
| `signature_comparison.r` | Cross-dataset concordance: Cohort 2, mouse, bulk Cohort 3, external CSF data |

See [`PBICs_analysis/README.md`](PBICs_analysis/README.md) for execution order and input/output details.

---

### 3. Tensor-S2R communication pipeline (`tensorS2R/`)

End-to-end ligand-receptor communication pipeline combining LIANA+ and Tensor-cell2cell.

The pipeline implements four steps:
1. **LR pair filtering** — top 70th percentile loadings; ligands/receptors filtered by scran-based DE (preferentially expressed on sender/receiver vs. other cell types).
2. **Pathway enrichment** — Enrichr (GSEApy) on filtered LR sets against WikiPathways, Reactome, and GO Biological Process; common pathways extracted.
3. **Over-representation analysis (ORA)** — decoupleR ORA on common enriched pathways at cell and pseudo-bulk level.
4. **Pseudo-bulk LR product** — mean pseudo-bulk expression of ligand in sender × receptor in receiver cells.

Three analytical designs are available:

| Mode | Data | Key interaction |
|------|------|----------------|
| `micro_blood_coarse` | Blood PBMCs + microglia (brain) | Monocytes ↔ Microglia |
| `brain_coarse` | Brain only (Cohort 2) | Monocytes ↔ Microglia (brain) |
| `brain_blood_coarse` | Full blood + all brain cell types | All PBICs ↔ all brain types |

```bash
python -m tensorS2R \
    -p <project_path> \
    -b <cohort_1|cohort_2> \
    -f <micro_blood_coarse|brain_coarse|brain_blood_coarse> \
    [-s <sample_key>] \
    [-g <groupby>] \
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
