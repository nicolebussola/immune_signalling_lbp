# tensorS2R Communication Pipeline

End-to-end pipeline for ligand-receptor communication analysis using Tensor-cell2cell and LIANA. Supports three human factorization modes (Cohorts 1 & 2) and a dedicated mouse pipeline (GEO GSE225948).

---

## Human pipeline

The full pipeline runs three sequential steps. You can run all steps at once or each independently.

### Run all steps

```bash
python -m tensorS2R \
    -p <project_path> \
    -b <cohort_1|cohort_2> \
    -f <factorization_type> \
    [--step <1|2|3|all>] \
    [-s <sample_key>] \
    [-g <groupby>] \
    [-m <merging_mode>] \
    [-r <rank>] \
    [-d <cpu|cuda>] \
    [-t <ora_threshold>]
```

| Argument | Flag | Default | Description |
|---|---|---|---|
| `project_path` | `-p` | — | Root directory where data are stored |
| `cohort` | `-b` | — | `cohort_1` or `cohort_2` |
| `factorization_type` | `-f` | — | See modes below (required for steps 1 and 2) |
| `step` | `--step` | `all` | `1`, `2`, `3`, or `all` |
| `sample_key` | `-s` | `sample` | Patient column in `adata.obs` |
| `groupby` | `-g` | `cell_type_coarse` | Cell-type column for LIANA predictions |
| `merging_mode` | `-m` | `inner` | Sample merging mode for tensor construction |
| `factorization_rank` | `-r` | auto (elbow) | Tensor rank; required for `--step 2` when step 1 was run separately |
| `device` | `-d` | `cuda` | `cuda` or `cpu` |
| `threshold_targets` | `-t` | `0.1` | Top fraction of genes for ORA (step 3) |

---

## Pipeline steps

### Step 1 — Tensor factorization

```bash
python -m tensorS2R.step1_factorize \
    -p <project_path> -b <cohort> -f <factorization_type> [-s] [-g] [-m] [-r] [-d]
```

Runs LIANA rank-aggregate by sample then Tensor-cell2cell factorization. Writes all outputs (loadings, factor plots, network, Gini) to:

```
<project_path>/<cohort>/c2c_liana_outputs/<factorization_type>/rank_<r>/<merging_mode>/<groupby>/
```

### Step 2 — LR filtering and enrichment

```bash
python -m tensorS2R.step2_enrich \
    -p <project_path> -b <cohort> -f <factorization_type> -r <rank> [-m] [-g]
```

Filters LR pairs from tensor loadings and runs Enrichr-based pathway enrichment. Dispatches based on cohort and design:

- **All cohorts / all modes** — Monocyte↔Microglia enrichment
- **Cohort 2 + `brain_coarse`** — also Neuron↔Microglia enrichment
- **Cohort 2 + `brain_blood_coarse`** — also Neuron↔Monocyte enrichment

Factor numbers are configured in `CELLS_FACTOR_DICT` and `NEURO_FACTOR_DICT` at the top of `step2_enrich.py`. Update them when re-running factorization produces a different ordering.

### Step 3 — ORA and LR product

```bash
python -m tensorS2R.step3_downstream \
    -p <project_path> -b <cohort> [-t <threshold>] \
    [-f <factorization_type> -r <rank> [-m <merging_mode>] [-g <groupby>]]
```

Runs pseudo-bulk ORA via decoupler on MSigDB gene sets (WikiPathways, Reactome, GO Biological Process). Produces boxplot/swarmplot figures per pathway and interactive UMAP embeddings coloured by ORA score. When `-f` is provided, also computes pseudo-bulk LR product heatmaps for the filtered Mono↔Micro LR pairs written by step 2 (`lr_loads_mono_to_micro_filtered.csv`, `lr_loads_micro_to_mono_filtered.csv`). Outputs go to `<project_path>/<cohort>/ora_outputs/`.

---

## Human factorization modes

### `brain_coarse` — Monocytes and Microglia in brain only

Uses brain AnnData only (both cell types present in brain data).

**Step 1** (`step1_factorize.py`): LIANA rank-aggregate by sample → Tensor-cell2cell factorization on brain cell types.

**Step 2** (`step2_enrich.py` → `enrich_mono_micro_brain`):
- Extract LR pairs from Mono→Micro and Micro→Mono factors
- Filter ligands/receptors by differential expression on Microglia and Monocytes
- Enrichment background: genes expressed in ≥10% of Microglia or Monocytes in brain data
- Output: enriched pathways shared between ligand and receptor gene sets

**Step 2 — Cohort 2 only** (`enrich_neuro_micro`):
- Extract LR pairs from Neuron→Microglia (Factor 9) and Microglia→Neuron (Factor 10)
- Filter neuron ligands by DE vs. OPC; microglia receptors by DE vs. non-microglial types
- Enrichment (LR product heatmaps computed in step 3)

---

### `micro_blood_coarse` — Microglia (brain) and Monocytes (blood)

Uses a combined blood+microglia AnnData (microglia only from brain, `micro_only=True`).

**Step 1** (`step1_factorize.py`): LIANA rank-aggregate by sample → Tensor-cell2cell factorization on blood+microglia.

**Step 2** (`step2_enrich.py` → `enrich_mono_micro_blood_micro`):
- Requires pre-computed DEG AnnData (`blood_microglia_tf.h5ad`) in `c2c_liana_outputs/`
- Filter ligands expressed on CD16/CD14 Monocytes in both blood and combined adatas
- Filter receptors expressed on Microglia in brain data
- Enrichment background: genes expressed in ≥10% of Microglia (brain) or Monocytes (blood)

---

### `brain_blood_coarse` — Full brain and blood

Uses the full combined blood+brain AnnData (`micro_only=False`).

**Step 1** (`step1_factorize.py`): LIANA rank-aggregate by sample → Tensor-cell2cell factorization on all blood and brain cell types.

**Step 2** (`step2_enrich.py` → `enrich_mono_micro_blood_brain`):
- Requires pre-computed DEG AnnData (`blood_brain_tf.h5ad`) in `c2c_liana_outputs/`
- Applies additional filtering based on receiver/sender cell context (brain vs. blood)
- Enrichment background: genes expressed in ≥10% of Microglia (brain) or Monocytes (blood)

**Step 2 — Cohort 2 only** (`enrich_neuro_mono`):
- Extract LR pairs from Neuron→Monocyte (Factor 6) and Monocyte→Neuron (Factor 16)
- Filter neuron ligands by DE across all three neuron references (GABAergic, Glutamatergic, Neuron-low count)
- Filter monocyte receptors/ligands by DE in CD14 and CD16 monocytes vs. all other blood cell types

---

## Mouse pipeline

Runs LIANA mouse consensus + Tensor-cell2cell on the paired brain/blood mouse dataset (GEO GSE225948; 11 mice; sham, D02, D14 conditions). Requires the output of `preprocessing/process_mouse_data.py` (`paired_blood_brain_mouse.h5ad`).

```bash
python -m tensorS2R.factorization_run_mouse \
    -p <mouse_processed_path> \
    [-s <sample_key>] \
    [-g <groupby>] \
    [-f <factorization_type>] \
    [-m <merging_mode>] \
    [-r <rank>] \
    [-d <cpu|cuda>]
```

| Argument | Flag | Default | Description |
|---|---|---|---|
| `project_path` | `-p` | — | Path to `mouse_processed/` directory |
| `sample_key` | `-s` | `sample` | Sample column in `adata.obs` |
| `groupby` | `-g` | `parent` | Cell-type column for LIANA predictions (`parent` = coarse, `sub.celltype` = fine) |
| `factorization_type` | `-f` | `mouse_coarse` | Label used as output subfolder name |
| `merging_mode` | `-m` | `inner` | Tensor construction merging mode |
| `factorization_rank` | `-r` | auto | Fixed at 20 when >7 cell types (paper default) |
| `device` | `-d` | `cuda` | `cuda` or `cpu` |

Cell-type columns available in `adata.obs` after loading:

| Column | Description |
|---|---|
| `parent` | Coarse cell type (e.g. `Mo+`, `Mo-`, `MdC`) |
| `sub.celltype` | Fine-grained cell type (e.g. `Mo1`, `Mo2`, `MdC2`) |
| `cell_type.v2` | `parent_tissue` (e.g. `Mo+_blood`, `MdC_brain`) |
| `cell_type.v3` | `sub.celltype_tissue` |

---

## Module structure

| File | Purpose |
|---|---|
| `__main__.py` | Human pipeline entry point — dispatches steps 1, 2, 3, or all |
| `step1_factorize.py` | `run_tensorcell2cell()` — LIANA + Tensor-cell2cell (human) |
| `step2_enrich.py` | `run_enrichment()` — LR filtering and Enrichr pathway enrichment (all modes) |
| `step3_ora.py` | `main()` — pseudo-bulk ORA via decoupler on MSigDB gene sets |
| `factorization_run_mouse.py` | `run_tensorcell2cell_mouse()` — LIANA mouse consensus + Tensor-cell2cell |
| `adata_processing_utils.py` | Data loading, filtering, and preprocessing |
| `lr_loadings_utils.py` | LR pair processing, DE filtering, pseudo-bulk product, enrichment utilities |
| `ora_analysis.py` | Legacy ORA module (superseded by `step3_ora.py`) |
| `plot_utils.py` | Interactive Bokeh embedding plots; `plot_pathway_network` for enrichment maps |

---

## Output structure

**Human:**
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
        Net_paths_mono_to_micro.png
        Net_paths_micro_to_mono.png
    [cohort_2 + brain_coarse only]
    lr_loads_neuro_to_micro_filtered_9.csv
    lr_loads_micro_to_neuro_filtered_10.csv
    Enriched_paths_neuro_to_micro_9.csv
    Enriched_paths_micro_to_neuro_10.csv
    Net_paths_neuro_to_micro_9.png
    Net_paths_micro_to_neuro_10.png
    [cohort_2 + brain_blood_coarse only]
    lr_loads_neuro_to_mono_filtered.csv
    lr_loads_mono_to_neuro_filtered.csv
    Enriched_paths_neuro_to_mono.csv
    Enriched_paths_mono_to_neuro.csv
    Net_paths_neuro_to_mono.png
    Net_paths_mono_to_neuro.png

<project_path>/<cohort>/ora_outputs/
    pseudobulk_ora_<db>_<pathway>.pdf
    [when -f is provided]
    lr_product_mono_to_micro.pdf
    lr_product_micro_to_mono.pdf
    [when -f brain_coarse + cohort_2]
    lr_product_neuro_to_micro_<factor>.pdf
    lr_product_micro_to_neuro_<factor>.pdf
```

**Mouse:**
```
<mouse_processed_path>/c2c_liana_outputs/<factorization_type>/rank_<r>/<merging_mode>/<groupby>/
    Loadings.xlsx
    Factorization_mouse.pdf
    clustermap_mouse.pdf
    Clustermap-LRs_mouse.pdf
    Gini_mouse.csv
    Networks_mouse.pdf
    liana_agg_mouse.pdf
```
