# blood-brain-LBP



## Quality control per sample

1. Load CellRanger raw count matrix (`raw_feature_bc_matrix.h5`)
2. Remove ambient RNA with `CellBender` <a name="cite_ref-1"></a>[<sup>[1]</sup>](#cite_note-1)
3. Basic filtering: \
        * Filter out genes expressed in < 3 cells \
        * Filter cells with <200 genes detected \
        * Remove `MALAT1` gene 
7. Compute QC metric (`pct_counts_mt`, `pct_counts_ribo`)
8. Detect doublets (`scDblFinder` and `scrublet`)
9. Feature selection: most deviant genes on raw counts <a name="cite_ref-2"></a>[<sup>[2]</sup>](#cite_note-2)
10. PCA and UMAP
11. Select thresholds for QC metric and doublets

<a name="cite_note-1"></a>1. [^](#cite_ref-1) Janssen P. The efect of background noise and its removal on the analysis of single-cell expression data, *Genome Biology* (June 2023) 
> we use our genotype-based estimates to evaluate the performance of three methods (CellBender, DecontX, SoupX) that are designed to quantify and remove background noise. We fnd that CellBender provides the most precise estimates of background noise levels and also yields the highest improvement for marker gene detection. By contrast, clustering and classifcation of cells are fairly robust towards background noise and only small improvements can be achieved by background removal that may come at the cost of distortions in fne structure.

<a name="cite_note-2"></a>2. [^](#cite_ref-2) Germain P.L., pipeComp, a general framework for the evaluation of computational pipelines, reveals performant single cell rna-seq preprocessing tools. *Genome Biology* (September 2020)