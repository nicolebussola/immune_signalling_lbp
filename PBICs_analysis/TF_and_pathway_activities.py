from pathlib import Path

import decoupler as dc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

DATA_DIR = Path(".")
OUTPUT_DIR = Path("figures")
OUTPUT_DIR.mkdir(exist_ok=True)


CELL_TYPES = [
    ("cd14", "CD14+ Monocytes", "consensus"),
    ("cd16", "CD16+ Monocytes", "consensus"),
    ("cd8", "CD8+ T cells", "consensus"),
    ("cd4", "CD4+ T cells", "consensus"),
    ("nk", "CD16+ NK", "consensus"),
    ("bc", "B cells", "ulm"),
]


def load_mat(stem: str, label: str) -> pd.DataFrame:
    """Load and preprocess a dreamlet contrast matrix for one cell type."""
    mat = pd.read_csv(DATA_DIR / f"mat_{stem}_dreamlet.csv").dropna(axis="columns")
    mat.index = mat["Unnamed: 0"]
    mat = mat.iloc[:, 1:].drop(index=["(Intercept)", "chemistryv3"])
    mat.index = mat.index.str.replace("tissuebrain", label, regex=False)
    return mat[[c for c in mat.columns if not c.startswith("MT-")]].copy()


def run_activity(mat: pd.DataFrame, net: pd.DataFrame, method: str):
    """Dispatch to run_consensus or run_ulm based on method string."""
    if method == "consensus":
        return dc.run_consensus(mat=mat, net=net)
    return dc.run_ulm(mat=mat, net=net)


def pval_annotations(pval_df: pd.DataFrame) -> np.ndarray:
    """Return a significance-star annotation array from a p-value DataFrame."""
    ann = np.full(pval_df.shape, "", dtype=object)
    v = pval_df.values
    ann[v < 0.001] = "***"
    ann[(v >= 0.001) & (v < 0.01)] = "**"
    ann[(v >= 0.01) & (v < 0.05)] = "*"
    return ann


VOLCANO_CONFIG = {
    "cd14": {
        "highlight_genes": ["DNAJA1", "HSP90AA1", "JUNB", "JUND", "UBE2S"],
        "manual_offsets": {
            "DNAJA1": (0.5, +2.0),
            "HSP90AA1": (+1.0, +0.6),
            "JUNB": (+0.8, -0.8),
            "JUND": (-1.5, +0.05),
            "UBE2S": (+3.0, +0.2),
        },
    },
    "cd16": {"highlight_genes": [], "manual_offsets": {}},
    "cd8": {"highlight_genes": [], "manual_offsets": {}},
    "cd4": {"highlight_genes": [], "manual_offsets": {}},
    "nk": {"highlight_genes": [], "manual_offsets": {}},
    "bc": {"highlight_genes": [], "manual_offsets": {}},
}


def plot_volcano(result_df, title, out_path, highlight_genes=(), manual_offsets=None):
    """Volcano plot with optional manually positioned gene labels."""
    if manual_offsets is None:
        manual_offsets = {}

    dc.plot_volcano_df(
        result_df,
        x="logFC",
        y="adj.P.Val",
        top=0,
        figsize=(10, 10),
        color_pos="#F8083D",
        color_neg="#045FC4",
        color_null="#EFECF4",
        lFCs_thr=1,
    )
    ax = plt.gca()

    # Hide auto-generated labels for genes we'll redraw manually
    existing_labels = {t.get_text(): t for t in ax.texts}
    for gene in highlight_genes:
        if gene in existing_labels:
            existing_labels[gene].set_visible(False)

    for gene in highlight_genes:
        if gene not in result_df.index:
            continue
        xd = result_df.loc[gene, "logFC"]
        yd = -np.log10(result_df.loc[gene, "adj.P.Val"])
        dx, dy = manual_offsets.get(gene, (0.2, 0.2))
        ax.scatter(
            xd, yd, s=60, facecolor="none", edgecolor="black", linewidth=2.0, zorder=10
        )
        ax.annotate(
            gene,
            xy=(xd, yd),
            xytext=(xd + dx, yd + dy),
            fontsize=18,
            fontweight="bold",
            color="black",
            bbox=dict(
                boxstyle="round,pad=0.2", facecolor="white", edgecolor="red", alpha=0.5
            ),
            arrowprops=dict(arrowstyle="-", color="black", lw=1.2),
        )

    for text in ax.texts:
        text.set_fontsize(18)

    plt.title(title, fontsize=30)
    plt.xlabel("logFCs", fontsize=18)
    plt.ylabel("-log10(pvals)", fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()


results_dfs = {
    stem: pd.read_csv(DATA_DIR / f"df_{stem}_dreamlet.csv", index_col="Unnamed: 0")
    for stem, *_ in CELL_TYPES
}

mats = {stem: load_mat(stem, label) for stem, label, _ in CELL_TYPES}

for stem, label, _ in CELL_TYPES:
    cfg = VOLCANO_CONFIG[stem]
    plot_volcano(
        results_dfs[stem],
        title=label,
        out_path=OUTPUT_DIR / f"{stem.upper()}_volcano.pdf",
        highlight_genes=cfg["highlight_genes"],
        manual_offsets=cfg["manual_offsets"],
    )


#  TF activity (CollecTRI)
collectri = dc.get_collectri(organism="human", split_complexes=False)

tf_acts, tf_pvals = {}, {}
for stem, _, method in CELL_TYPES:
    tf_acts[stem], tf_pvals[stem] = run_activity(mats[stem], collectri, method)

# B cells excluded from heatmap (ulm vs consensus incompatibility)
tf_stems = ["cd4", "cd8", "nk", "cd16", "cd14"]


sig_acts = {s: tf_acts[s][tf_pvals[s] < 0.05] for s in tf_stems}
tf_matrix = pd.concat(sig_acts.values()).T
tf_matrix = tf_matrix.dropna(how="all")
pval_matrix = pd.concat([tf_pvals[s] for s in tf_stems]).T.loc[tf_matrix.index]

fig, ax = plt.subplots(figsize=(15, 6))
sns.heatmap(
    tf_matrix.T,
    ax=ax,
    cmap="icefire",
    center=0,
    linewidths=0.5,
    annot=pval_annotations(pval_matrix.T),
    fmt="",
    cbar_kws={"orientation": "vertical", "shrink": 1},
)
ax.collections[0].cmap.set_bad("0.8")
plt.xticks(rotation=90)
plt.xlabel("TF", fontsize=15)
plt.ylabel("Cell type", fontsize=15)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "collectri_heatmap_bq.png")
plt.close()


# Pathway activity (PROGENy)
progeny = dc.get_progeny(top=5000)
selected_pathways = ["Estrogen", "Hypoxia", "MAPK", "NFkB", "TNFa"]

pw_acts, pw_pvals = {}, {}
for stem, _, _ in CELL_TYPES:
    acts, pvals = dc.run_mlm(mat=mats[stem], net=progeny)
    pw_acts[stem] = acts[selected_pathways]
    pw_pvals[stem] = pvals[selected_pathways]

pathway_acts = pd.concat(pw_acts.values())
pval_df = pd.concat(pw_pvals.values())

fig, ax = plt.subplots(figsize=(10, 7))
sns.heatmap(
    pathway_acts,
    ax=ax,
    cmap="coolwarm",
    center=0,
    linewidths=4,
    annot=pval_annotations(pval_df),
    annot_kws={"fontsize": 15},
    fmt="",
    cbar_kws={"orientation": "vertical", "shrink": 1},
)
plt.xticks(fontsize=18, rotation=0)
plt.yticks(fontsize=20, rotation=0)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "progeny_heatmap_human_subset.png", dpi=600)
plt.close()
