import doubletdetection
import rpy2.robjects as ro
import scrublet as scr
import scvi
from rpy2.robjects.packages import importr

importr("Seurat")
importr("scater")
importr("scDblFinder")
importr("BiocParallel")
importr("scry")
importr("scds")


def scdblfinder(adata):
    data_mat = adata.X.T.copy()
    ro.r(
        """
        scDblFinder_score <- function(data_mat){
        set.seed(42)
        sce = scDblFinder(
            SingleCellExperiment(
                list(counts=data_mat),
            )
        )
        doublet_score = sce$scDblFinder.score
        return(doublet_score)}
            """
    )
    ro.r(
        """
        scDblFinder_class <- function(data_mat){
        set.seed(42)
        sce = scDblFinder(
            SingleCellExperiment(
                list(counts=data_mat),
            )
        )
        doublet_class = sce$scDblFinder.class
        return(doublet_class)}
            """
    )
    adata.obs["doublet_scores_scdblfinder"] = ro.globalenv["scDblFinder_score"](
        data_mat
    )
    adata.obs["predicted_doublets_scdblfinder"] = ro.globalenv["scDblFinder_class"](
        data_mat
    ).astype("object")
    return adata


def scrublet(adata):
    dbl_rate = adata.X.shape[0] / 1000 * 0.008  # from demuxify
    scrub = scr.Scrublet(adata.X, expected_doublet_rate=dbl_rate, sim_doublet_ratio=2)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(
        min_counts=3, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30
    )

    adata.obs["doublet_scores_scrublet"] = doublet_scores
    adata.obs["predicted_doublets_scrublet"] = predicted_doublets
    return adata


def scds(adata):
    data_mat = adata.X.T.copy()
    scds = ro.r(
        """
        scds <- function(data_mat){
        set.seed(42)
        sce = cxds(SingleCellExperiment(
                list(counts=data_mat),
            ),retRes = TRUE)
        sce = bcds(sce,retRes = TRUE,verb=TRUE)
        sce = cxds_bcds_hybrid(sce, estNdbl=TRUE)
        }
        """
    )
    doublet_results = ro.globalenv["scds"](data_mat)
    adata.obs["predicted_doublets_scds"] = list(doublet_results.obs["hybrid_call"])
    adata.obs["doublet_scores_scds"] = list(doublet_results.obs["hybrid_score"])
    return adata


def doubletdetection_method(adata, output_path_plot, name):
    clf = doubletdetection.BoostClassifier(
        n_iters=50,
        clustering_algorithm="louvain",
        standard_scaling=True,
        pseudocount=0.1,
        n_jobs=-1,
        random_state=42,
    )

    doublets_clf = clf.fit(adata.X).predict(p_thresh=1e-16, voter_thresh=0.5)
    doublet_score = clf.doublet_score()
    adata.obs["predicted_doublets_doubletdetection"] = doublets_clf
    adata.obs["doublet_scores_doubletdetection"] = doublet_score
    doubletdetection.plot.convergence(
        clf,
        show=False,
        save=f"{output_path_plot}/convergence_test_{name}.pdf",
        p_thresh=1e-16,
        voter_thresh=0.5,
    )
    return adata


def solo(adata):
    scvi.model.SCVI.setup_anndata(adata)
    vae = scvi.model.SCVI(adata)
    vae.train(early_stopping=True)
    solo = scvi.external.SOLO.from_scvi_model(vae)
    solo.train()
    doublets_solo = solo.predict(soft=False).tolist()
    adata.obs["predicted_doublets_solo"] = doublets_solo
    return adata
