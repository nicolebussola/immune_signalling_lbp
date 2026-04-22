QC_LABELS_SAMPLE = [
    "log1p_total_counts",
    "log1p_n_genes_by_counts",
    "pct_counts_in_top_20_genes",
    "pct_counts_mt",
    "pct_counts_ribo",
    "cluster_labels",
    "passed_qc",
]

CHEMISTRY_V2_PATIENTS = ["182", "185"]

QC_LABELS_BATCH_2 = [
    "pt",
    "side",
    "log1p_total_counts",
    "log1p_n_genes_by_counts",
    "pct_counts_in_top_20_genes",
    "pct_counts_mt",
    "pct_counts_ribo",
]


QC_LABELS_BATCH_1 = QC_LABELS_BATCH_2 + ["chemistry"]

FAILED_QC_SAMPLES = ["PT-292-brain-L_QC.h5ad"]
