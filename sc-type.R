
# load libraries
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)

# Load the blood dataset
blood.data <- Read10X_h5( "Desktop/Projects/LBP_MountSinai/LBP_brain_blood_pairs/data/narsad_cellRanger_outs/blood/PT-205-blood-L/PT-205-L-B_CellBender_filtered.h5")
# Initialize the Seurat object with the raw (non-normalized data).
blood <- CreateSeuratObject(counts = blood.data, min.cells = 3, min.features = 200)
#blood<- blood[ ! grepl("MALAT1", rownames(blood)), ]
library(scater)
library(scDblFinder)
library(BiocParallel)



# normalize data
blood[["percent.mt"]] <- PercentageFeatureSet(blood, pattern = "^MT-")
blood <- subset(blood, percent.mt < 20) # make some filtering based on QC metrics visualizations, see Seurat tutorial: https://satijalab.org/seurat/articles/blood3k_tutorial.html
blood <- NormalizeData(blood, normalization.method = "LogNormalize", scale.factor = 10000)
blood <- FindVariableFeatures(blood, selection.method = "vst", nfeatures = 5000)

# scale and run PCA
blood <- ScaleData(blood, features = rownames(blood))
blood <- RunPCA(blood, features = VariableFeatures(object = blood))

# Check number of PC components (we selected 10 PCs for downstream analysis, based on Elbow plot)
ElbowPlot(blood)

# cluster and visualize
blood <- FindNeighbors(blood, dims = 1:10)
blood <- FindClusters(blood, resolution = 0.8)
blood <- RunUMAP(blood, dims = 1:10)
DimPlot(blood, reduction = "umap")

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)


# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = blood[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either blood[["RNA"]]@scale.data (default), blood[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or blood[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(blood@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(blood@meta.data[blood@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(blood@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

blood@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  blood@meta.data$customclassif[blood@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(blood, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')        


# load auto-detection function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R")
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";

# guess a tissue type
tissue_guess = auto_detect_tissue_type(path_to_db_file = db_, seuratObject = blood, scaled = TRUE, assay = "RNA")  # if saled = TRUE, make sure the data is scaled, as seuratObject[[assay]]@scale.data is used. If you just created a Seurat object, without any scaling and normalization, set scaled = FALSE, seuratObject[[assay]]@counts will be used         
