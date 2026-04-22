library(zellkonverter)
library(dreamlet)
library(SingleCellExperiment)
library(reshape2)
library(ggplot2)
library(dplyr)
library(ggbeeswarm)
library(plotly)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(Orthology.eg.db)

MOUSE_H5AD <- "/path/to/paired_blood_brain_mouse.h5ad"
HUMAN_H5AD <- "/path/to/paired_blood_brain_cohort1.h5ad"


mapIt <- function(mouseids, horg, morg, orth) {
  mouseg <- mapIds(morg, mouseids, "ENTREZID", "SYMBOL")
  mapped <- select(orth, mouseg, "Homo_sapiens", "Mus_musculus")
  names(mapped) <- c("Mus_egid", "Homo_egid")
  husymb <- select(horg, as.character(mapped[, 2]), "SYMBOL", "ENTREZID")
  data.frame(Mus_symbol = mouseids, mapped, Homo_symbol = husymb[, 2])
}


mouse_h5 <- readH5AD(MOUSE_H5AD, verbose = TRUE, reader = "R")
human_h5  <- readH5AD(HUMAN_H5AD, verbose = TRUE)

gene_mapped <- mapIt(row.names(mouse_h5), org.Hs.eg.db, org.Mm.eg.db, Orthology.eg.db)


# HUMAN 

human_h5$id <- paste0(human_h5$pt, human_h5$tissue)

pb_human <- aggregateToPseudoBulk(
  human_h5,
  assay      = "counts",
  cluster_id = "cell_type.v2",
  sample_id  = "id",
  verbose    = FALSE
)

res.proc_human <- processAssays(pb_human, ~tissue + chemistry + (1 | pt), min.count = 5)

vp.lst_human <- fitVarPart(res.proc_human, ~tissue + chemistry + (1 | pt))
plotVarPart(vp.lst_human, label.angle = 60)

res.dl_human <- dreamlet(res.proc_human, ~tissue + chemistry + (1 | pt))

write.csv(topTable(res.dl_human[["CD14.mono"]], coef = "tissuebrain", number = Inf), "df_cd14_dreamlet.csv")
write.csv(topTable(res.dl_human[["CD16.mono"]], coef = "tissuebrain", number = Inf), "df_cd16_dreamlet.csv")
write.csv(topTable(res.dl_human[["CD4.T"]],    coef = "tissuebrain", number = Inf), "df_cd4_dreamlet.csv")
write.csv(topTable(res.dl_human[["CD8.T"]],    coef = "tissuebrain", number = Inf), "df_cd8_dreamlet.csv")
write.csv(topTable(res.dl_human[["CD16.NK"]],  coef = "tissuebrain", number = Inf), "df_nk_dreamlet.csv")
write.csv(topTable(res.dl_human[["B"]],        coef = "tissuebrain", number = Inf), "df_bc_dreamlet.csv")

write.csv(t(res.dl_human[["CD14.mono"]]$t), "mat_cd14_dreamlet.csv")
write.csv(t(res.dl_human[["CD16.mono"]]$t), "mat_cd16_dreamlet.csv")
write.csv(t(res.dl_human[["CD4.T"]]$t),     "mat_cd4_dreamlet.csv")
write.csv(t(res.dl_human[["CD8.T"]]$t),     "mat_cd8_dreamlet.csv")
write.csv(t(res.dl_human[["CD16.NK"]]$t),   "mat_nk_dreamlet.csv")
write.csv(t(res.dl_human[["B"]]$t),         "mat_bc_dreamlet.csv")


#  MOUSE
mouse_h5$cell_type.v2 <- as.character(mouse_h5$cell_type.v2)
mouse_h5$parent       <- as.character(mouse_h5$parent)

mouse_h5$parent[(mouse_h5$sub.celltype == "Tc2") & (mouse_h5$tissue == "PB")]     <- "CD4.T"
mouse_h5$parent[(mouse_h5$sub.celltype == "Tc1") & (mouse_h5$tissue == "PB")]     <- "CD8.T"
mouse_h5$parent[(mouse_h5$sub.celltype == "Tc2") & (mouse_h5$tissue == "brain")]  <- "CD8.T"
mouse_h5$parent[(mouse_h5$sub.celltype == "Tc1") & (mouse_h5$tissue == "brain")]  <- "CD4.T"

mouse_h5$parent[(mouse_h5$sub.celltype == "Mo1")  & (mouse_h5$tissue == "PB")]    <- "CD14.mono"
mouse_h5$parent[(mouse_h5$sub.celltype == "Mo2")  & (mouse_h5$tissue == "PB")]    <- "CD14.mono"
mouse_h5$parent[(mouse_h5$sub.celltype == "MdC2") & (mouse_h5$tissue == "brain")] <- "CD14.mono"
mouse_h5$parent[(mouse_h5$sub.celltype == "MdC3") & (mouse_h5$tissue == "brain")] <- "CD14.mono"
mouse_h5$parent[(mouse_h5$sub.celltype == "MdC4") & (mouse_h5$tissue == "brain")] <- "CD14.mono"

mouse_h5$cell_type.v2[mouse_h5$cell_type.v2 == "NK_brain"] <- "NK"
mouse_h5$cell_type.v2[mouse_h5$cell_type.v2 == "NK_PB"]    <- "NK"
mouse_h5$parent[mouse_h5$parent == "MdC"] <- "Mo"
mouse_h5$id <- paste0(mouse_h5$treatment, mouse_h5$sample, mouse_h5$tissue)

obs.mouse <- as.data.frame(colData(mouse_h5)) %>% select("treatment", "parent", "tissue")

mode(mouse_h5@assays@data$counts) <- "integer"

pb_mouse <- aggregateToPseudoBulk(
  mouse_h5,
  assay      = "counts",
  cluster_id = "parent",
  sample_id  = "id",
  verbose    = TRUE
)

res.proc_mouse <- processAssays(
  pb_mouse,
  assays      = c("Bc", "CD14.mono", "CD4.T", "CD8.T"),
  ~tissue + treatment + (1 | sample),
  min.cells   = 5,
  min.count   = 5,
  min.prop    = 0.1,
  min.samples = 1
)

res.dl_mouse <- dreamlet(res.proc_mouse, ~tissue + treatment + (1 | sample))

write.csv(topTable(res.dl_mouse[["CD14.mono"]], coef = "tissuebrain", number = Inf), "df_cd14_dreamlet_mouse.csv")
write.csv(topTable(res.dl_mouse[["CD4.T"]],    coef = "tissuebrain", number = Inf), "df_cd4_dreamlet_mouse.csv")
write.csv(topTable(res.dl_mouse[["CD8.T"]],    coef = "tissuebrain", number = Inf), "df_cd8_dreamlet_mouse.csv")
write.csv(topTable(res.dl_mouse[["Bc"]],       coef = "tissuebrain", number = Inf), "df_bc_dreamlet_mouse.csv")

write.csv(t(res.dl_mouse[["CD14.mono"]]$t), "mat_cd14_dreamlet_mouse.csv")
write.csv(t(res.dl_mouse[["CD4.T"]]$t),     "mat_cd4_dreamlet_mouse.csv")
write.csv(t(res.dl_mouse[["CD8.T"]]$t),     "mat_cd8_dreamlet_mouse.csv")



# Boxplot: Ccl5 in CD8 T cells 
df.ccl5 <- extractData(res.proc_mouse, "CD8.T", genes = "Ccl5")
thm <- theme_minimal() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 12))

ggplot(df.ccl5, aes(tissue, Ccl5, fill = factor(treatment))) +
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(aes(color = factor(treatment)), stroke = 0.2, color = "black",
                alpha = 0.6, size = 1, dodge.width = 0.8) +
  thm +
  ylab(bquote(Expression ~ (log[2] ~ Ccl5))) +
  xlab("") +
  ggtitle("CD8+ T cells - Ccl5")


