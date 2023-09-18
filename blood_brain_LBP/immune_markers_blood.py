t_cell_markers = {
    "Regulatory T cells (Tregs)": ["FOXP3"],
    "Naive CD4+ T cells": ["CCR7", "CD45RA"],
    "Central Memory CD4+ T cells": ["CCR7", "CD45RO"],
    "Effector Memory CD4+ T cells": ["CD45RO"],
    "Follicular helper T cells (Tfh)": ["CXCR5"],
    "Cytotoxic CD8+ T cells": ["GZMB", "PRF1", "TBX21", "CD45RO"],
    "Naive CD8+ T cells": ["CCR7", "CD45RA"],
    "Central Memory CD8+ T cells": ["CCR7", "CD45RO"],
    "Effector Memory CD8+ T cells": ["CD45RO"],
    "Regulatory T cells (Activated Tregs)": ["CTLA4", "CD25"],
    "Tissue-Resident Memory T cells (TRM)": ["CD69", "CD103"],
    "Th17 cells": ["IL17A", "RORC"],
    "Th1 cells": ["IFNG", "TBX21"],
    "Th2 cells": ["IL4", "IL5", "IL13", "GATA3"],
}

b_cell_markers = {
    "Naive B cells": {
        "CD19": "B cell surface marker",
        "CD27": "Not expressed on naive B cells",
        "CD24": "Marker for immature B cells",
    },
    "Memory B cells": {
        "CD19": "B cell surface marker",
        "CD27": "Marker for memory B cells",
        "CD20": "B cell surface marker",
    },
    "Plasma cells": {
        "CD138": "Plasma cell marker (Syndecan-1)",
        "IRF4": "Transcription factor involved in plasma cell differentiation",
        "XBP1": "Transcription factor essential for plasma cell development",
    },
    "Marginal zone B cells": {
        "CD19": "B cell surface marker",
        "CD27": "Expressed on some marginal zone B cells",
        "CD21": "High expression on marginal zone B cells",
    },
    "B1 B cells": {
        "CD19": "B cell surface marker",
        "CD5": "Expressed on B1 B cells",
        "CD11b": "Expressed on some B1 B cells",
    },
    "Follicular B cells": {
        "CD19": "B cell surface marker",
        "CD23": "Marker for follicular B cells",
        "CD21": "Expressed on some follicular B cells",
    },
}

immune_cell_markers = {
    "B cells": {
        "Naive B cells": ["CD19", "CD27", "CD24"],
        "Memory B cells": ["CD19", "CD27", "CD20"],
        "Plasma cells": ["CD138", "IRF4", "XBP1"],
    },
    "Natural Killer (NK) cells": {
        "CD56bright NK cells": ["CD56", "CD16"],
        "CD56dim NK cells": ["CD56", "CD16", "KLRB1"],
    },
    "Monocytes": {
        "Classical Monocytes": ["CD14", "CD16"],
        "Non-classical Monocytes": ["CD14", "CD16", "CX3CR1"],
        "Intermediate Monocytes": ["CD14", "CD16", "CX3CR1", "CD163"],
    },
    "Dendritic cells (DCs)": {
        "Conventional DCs (cDCs)": ["CD11c", "HLA-DR", "CD86"],
        "Plasmacytoid DCs (pDCs)": ["CD123", "CD303", "CD304"],
    },
    "Neutrophils": {"Neutrophils": ["CD66b", "CD16", "CD11b"]},
    "Eosinophils": {"Eosinophils": ["CD125", "CCR3", "SIGLEC8"]},
    "Basophils": {"Basophils": ["CD203c", "CD123", "CD193"]},
}
