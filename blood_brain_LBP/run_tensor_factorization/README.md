**Cell-to-cell communications between Monocytes in blood and Microglia in brain - Batch 1 and Batch 2**

*Step 1*: Tensor factorization using Tensor-cell2cell factorization and LIANA ligand-receptor pairs

`factorization_c2c_liana.py`

inputs:
    * AnnData object created combined filtered blood and brain data (log1p_cell_counts>7 and cell_types with >100 cells) into a single object. 
    (For reproducibility, only common cell types in batch 1 and 2 are considered for factorization)
    * sample_key: Patient column in adata.obs to be used for liana ccc prediction by samples
    * groupby: Cell-type column in adata.obs to be used for liana ccc predictions

output:
    * Factor Loadings resulting from c2c tensor factorization


*Step 2*: Extract genes that summarize communication between Monocytes and Microglia 

`enrichment_mono_micro_factors.py`

*Step 2a*: Process loadings and extract LR pairs corresponding to Factors "From" and "To" Microglia
*Step 2b*: For each factor, extract genes within LR pairs that are DE on Microglia and Monocytes using original data (all genes and cell types). Specifically:
     - "From" Microglia factor: extract ligands that are differentially expressed (logfc=0.5) on Microglia; Extract receptors in LR pairs that interact with filtered ligands; Filter receptors that are differentially expressed on either CD16 or CD14 Monocytes. 
       - "To" Microglia factor: extract ligands that are differentially expressed (logfc=0.5) on either CD16 or CD14 Monocytes; Extract receptors in LR pairs that interact with filtered ligands; Filter receptors that are differentially expressed on Microglia.
*Step 2c*: Further filter LR pairs: as tensor factorization was computed on a simulated brain-blood tissue, we need to check that ligands in "From" factors and receptors in "To" factors are also not highly expressed on Microglia. Filter genes differentially expressed (logfc=0.5) on Microglia. Extract list of ligands and receptors for each factor.

*Step 2d*: Find pathways enriched that are enriched on both set of ligands and receptors. Enrichments is computed using as background genes:
  - "From" Microglia Ligands: genes expressed in at least 10% of Microglia cells in brain data
  - "From" Microglia Receptors: genes expressed in at least 10% of Monocytes cells in blood data
  - "To" Microglia Ligands: genes expressed in at least 10% of Monocytes cells in blood data
  - "To" Microglia Receptors: genes expressed in at least 10% of Microglia cells in brain data

**Cell-to-cell communications between Monocytes and Microglia in brain - Batch 1**

*Step 1*: Tensor factorization using Tensor-cell2cell factorization and LIANA ligand-receptor pairs

`factorization_c2c_liana.py`

inputs:
    * Filtered AnnData object (log1p_cell_counts>7 and cell_types with >100 cells). 
    * sample_key: Patient column in adata.obs to be used for liana ccc prediction by samples
    * groupby: Cell-type column in adata.obs to be used for liana ccc predictions

output:
    * Factor Loadings resulting from c2c tensor factorization


*Step 2*: Extract genes that summarize communication between Monocytes and Microglia 

`enrichment_mono_micro_factors.py`

*Step 2a*: Process loadings and extract LR pairs corresponding to Factors "From" and "To" Microglia
*Step 2b*: For each factor, extract genes within LR pairs that are DE on Microglia and Monocytes using original data (all genes and cell types). Specifically:
     - "From" Microglia factor: extract ligands that are differentially expressed (logfc=0.5) on Microglia; Extract receptors in LR pairs that interact with filtered ligands; Filter receptors that are differentially expressed on either CD16 or CD14 Monocytes. 
       - "To" Microglia factor: extract ligands that are differentially expressed (logfc=0.5) on either CD16 or CD14 Monocytes; Extract receptors in LR pairs that interact with filtered ligands; Filter receptors that are differentially expressed on Microglia. Extract list of ligands and receptors for each factor.
*Step 2c*: Find pathways enriched that are enriched on both set of ligands and receptors. 


**Cell-to-cell communications between blood and brain - Batch 1**


**Cell to cell communications between blood and brain - Batch 2**
