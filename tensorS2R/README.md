# tensorS2R Communication Pipeline

End-to-end pipeline for ligand-receptor communication analysis between Monocytes (blood) and Microglia (brain) using Tensor-cell2cell and LIANA.

## Running the pipeline

```bash
python -m tensorS2R \
    -p <project_path> \
    -b <cohort_1|cohort_2> \
    -f <factorization_type> \
    [-s <sample_key>] \
    [-g <groupby>] \
    [-m <merging_mode>] \
    [-r <rank>] \
    [-d <cpu|cuda>]
```

| Argument | Flag | Default | Description |
|---|---|---|---|
| `project_path` | `-p` | HPC path | Root directory where data are stored |
| `cohort` | `-b` | — | `cohort_1` or `cohort_2` |
| `factorization_type` | `-f` | — | See modes below |
| `sample_key` | `-s` | `sample` | Patient column in `adata.obs` |
| `groupby` | `-g` | `cell_type_coarse` | Cell-type column for LIANA predictions |
| `merging_mode` | `-m` | `inner` | Sample merging mode for tensor construction |
| `factorization_rank` | `-r` | auto (elbow) | Tensor rank; inferred automatically if omitted |
| `device` | `-d` | `cuda` | `cuda` or `cpu` |

---

## Factorization modes

### `brain_coarse` — Monocytes and Microglia in brain only

Uses brain AnnData only (both cell types present in brain data).

**Step 1** (`factorization_run.py`): LIANA rank-aggregate by sample → Tensor-cell2cell factorization on brain cell types.

**Step 2** (`enrichment_run.py` → `enrich_mono_micro_brain`):
- Extract LR pairs from Mono→Micro and Micro→Mono factors
- Filter ligands/receptors by differential expression on Microglia and Monocytes
- Enrichment background: genes expressed in ≥10% of Microglia or Monocytes in brain data
- Output: enriched pathways shared between ligand and receptor gene sets

---

### `micro_blood_coarse` — Microglia (brain) and Monocytes (blood)

Uses a combined blood+microglia AnnData (microglia only from brain, `micro_only=True`).

**Step 1** (`factorization_run.py`): LIANA rank-aggregate by sample → Tensor-cell2cell factorization on blood+microglia.

**Step 2** (`enrichment_run.py` → `enrich_mono_micro_blood_micro`):
- Requires pre-computed DEG AnnData (`blood_microglia_tf.h5ad`) in `c2c_liana_outputs/`
- Filter ligands expressed on CD16/CD14 Monocytes in both blood and combined adatas
- Filter receptors expressed on Microglia in brain data
- Enrichment background: genes expressed in ≥10% of Microglia (brain) or Monocytes (blood)

---

### `brain_blood_coarse` — Full brain and blood

Uses the full combined blood+brain AnnData (`micro_only=False`).

**Step 1** (`factorization_run.py`): LIANA rank-aggregate by sample → Tensor-cell2cell factorization on all blood and brain cell types.

**Step 2** (`enrichment_run.py` → `enrich_mono_micro_blood_brain`):
- Requires pre-computed DEG AnnData (`blood_brain_tf.h5ad`) in `c2c_liana_outputs/`
- Applies additional filtering based on receiver/sender cell context (brain vs. blood)
- Enrichment background: genes expressed in ≥10% of Microglia (brain) or Monocytes (blood)

---

## Module structure

| File | Purpose |
|---|---|
| `__main__.py` | Pipeline entry point — runs factorization then enrichment |
| `factorization_run.py` | `run_tensorcell2cell()` — LIANA + Tensor-cell2cell |
| `enrichment_run.py` | Enrichment functions for each factorization mode |
| `adata_processing_utils.py` | Data loading, filtering, and preprocessing |
| `lr_loadings_utils.py` | LR pair processing, DE filtering, enrichment utilities |
| `ora_analysis.py` | Over-representation analysis (ORA) with decoupler |
| `plot_utils.py` | Interactive Bokeh embedding plots |

## Output structure

```
<project_path>/<cohort>/c2c_liana_outputs/<factorization_type>/rank_<r>/<merging_mode>/<groupby>/
    Loadings.xlsx
    Factorization_<cohort>.pdf
    clustermap_<cohort>.pdf
    Clustermap-LRs_<cohort>.pdf
    Gini.csv
    Networks.pdf
    lr_loads_mono_to_micro_filtered.csv
    lr_loads_micro_to_mono_filtered.csv
    enrichr_results/
        Enriched_paths_mono_to_micro.csv
        Enriched_paths_micro_to_mono.csv
