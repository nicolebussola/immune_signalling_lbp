library(data.table)
library(qvalue)
library(ggplot2)
library(ggthemes)
library(patchwork)
library(tidyverse)
library(readxl)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(Orthology.eg.db)

DATA_DIR    <- "/path/to/paired_blood_brain_cohort1.h5ad"
BULK_RDS    <- "/path/to/bulk_data.RDS"
CSF_XLSX    <- "/path/to/PMID31937773_SD2.xlsx"
REP_XLSX    <- "/path/to/DataResource5.xlsx"
COHORT2_DIR  <-  "/path/to/sigs_cohort2"

HUMAN_CELL_TYPES <- c("bc", "cd14", "cd16", "cd4", "cd8", "nk")
MOUSE_CELL_TYPES <- c("bc", "cd14", "cd4", "cd8", "nk")
DREAMLET_COLS    <- c("V1", "P.Value", "adj.P.Val", "logFC", "t")

check_required_file <- function(path) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path, call. = FALSE)
  }
  path
}

check_required_cols <- function(dt, cols, label) {
  missing_cols <- setdiff(cols, colnames(dt))
  if (length(missing_cols) > 0) {
    stop(
      sprintf("Missing required columns in %s: %s", label, paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }
  invisible(dt)
}

safe_frac <- function(num, denom) {
  ifelse(denom > 0, num / denom, NA_real_)
}

read_dreamlet_table <- function(path, label) {
  dt <- fread(check_required_file(path))
  check_required_cols(dt, DREAMLET_COLS, label)
  dt
}
map_mouse_to_human <- function(mouse_symbols) {
  mouse_entrez <- mapIds(org.Mm.eg.db, mouse_symbols, "ENTREZID", "SYMBOL")
  orth <- select(Orthology.eg.db, mouse_entrez, "Homo_sapiens", "Mus_musculus")
  names(orth) <- c("Mus_egid", "Homo_egid")
  human_symbols <- select(org.Hs.eg.db, as.character(orth[, 2]), "SYMBOL", "ENTREZID")
  data.table(
    mouse_symbol = mouse_symbols,
    Mus_egid = orth$Mus_egid,
    Homo_egid = orth$Homo_egid,
    human_symbol = human_symbols[, 2]
  )
}

mydat <- list("human" = list(), "mouse" = list())

for (i in HUMAN_CELL_TYPES) {
  mydat[["human"]][[i]] <- read_dreamlet_table(
    file.path(".", paste0("df_", i, "_dreamlet.csv")),
    paste0("human ", i)
  )
}

for (i in MOUSE_CELL_TYPES) {
  add <- read_dreamlet_table(
    file.path(".", paste0("df_", i, "_dreamlet_mouse.csv")),
    paste0("mouse ", i)
  )
  colnames(add)[1] <- "mouse_symbol"
  mapper           <- map_mouse_to_human(add$mouse_symbol)
  mapper           <- mapper[!is.na(human_symbol)]
  keep1            <- mapper[, .N, mouse_symbol][N == 1]$mouse_symbol
  keep2            <- mapper[, .N, human_symbol][N == 1]$human_symbol
  mapper           <- mapper[mouse_symbol %in% keep1 & human_symbol %in% keep2]
  mydat[["mouse"]][[i]] <- merge(add, mapper[, .(mouse_symbol, human_symbol)])
}

human_deg_summary <- c()
for (i in HUMAN_CELL_TYPES) {
  cur <- mydat[["human"]][[i]]
  add <- data.table(
    cell   = i,
    ngene  = nrow(cur),
    pi1    = 1 - pi0est(cur$P.Value)$pi0,
    ndeg   = nrow(cur[adj.P.Val < 0.05]),
    ndegbr = nrow(cur[adj.P.Val < 0.05 & logFC > 0]),
    ndegbl = nrow(cur[adj.P.Val < 0.05 & logFC < 0])
  )
  human_deg_summary <- rbind(human_deg_summary, add)
}
human_deg_summary[, pctdeg    := safe_frac(ndeg, ngene)]
human_deg_summary[, pctdegbr  := safe_frac(ndegbr, ngene)]
human_deg_summary[, pctdegbl  := safe_frac(ndegbl, ngene)]
human_deg_summary[, pctdegbr2 := safe_frac(ndegbr, ndeg)]
human_deg_summary[, pctdegbl2 := safe_frac(ndegbl, ndeg)]
human_deg_summary

shared_human_deg_genes <- c()
for (i in HUMAN_CELL_TYPES) {
  cur <- copy(mydat[["human"]][[i]])
  cur[, brain_deg := FALSE]
  cur[, blood_deg := FALSE]
  cur[adj.P.Val < 0.05 & logFC > 0, brain_deg := TRUE]
  cur[adj.P.Val < 0.05 & logFC < 0, blood_deg := TRUE]
  shared_human_deg_genes <- rbind(shared_human_deg_genes,
    cur[brain_deg == TRUE, .(cell = i, symbol = V1, deg = "brain")],
    cur[blood_deg == TRUE, .(cell = i, symbol = V1, deg = "blood")])
}
nrow(shared_human_deg_genes[, .N, list(symbol, deg)][N >= 2])
shared_human_deg_genes[, .N, list(symbol, deg)][N >= 2][, .N, deg]
shared_human_deg_genes[, .N, list(symbol, deg)][order(N, decreasing = TRUE)][N == 6]

mouse_deg_summary <- c()
for (i in MOUSE_CELL_TYPES) {
  cur <- mydat[["mouse"]][[i]]
  add <- data.table(
    cell   = i,
    ngene  = nrow(cur),
    ndeg   = nrow(cur[adj.P.Val < 0.05]),
    ndegbr = nrow(cur[adj.P.Val < 0.05 & logFC > 0]),
    ndegbl = nrow(cur[adj.P.Val < 0.05 & logFC < 0])
  )
  mouse_deg_summary <- rbind(mouse_deg_summary, add)
}
mouse_deg_summary[, pctdeg    := safe_frac(ndeg, ngene)]
mouse_deg_summary[, pctdegbr  := safe_frac(ndegbr, ngene)]
mouse_deg_summary[, pctdegbl  := safe_frac(ndegbl, ngene)]
mouse_deg_summary[, pctdegbr2 := safe_frac(ndegbr, ndeg)]
mouse_deg_summary[, pctdegbl2 := safe_frac(ndegbl, ndeg)]
mouse_deg_summary

# ── Human vs mouse ────────────────────────────────────────────────────────────
human_mouse_correlations <- c()
human_mouse_plot_data    <- c()
for (i in MOUSE_CELL_TYPES) {
  cur1 <- mydat[["human"]][[i]]
  cur2 <- mydat[["mouse"]][[i]]
  cur3 <- merge(
    cur1[, .(human_symbol = V1, human_logFC = logFC, human_t = t)],
    cur2[, .(human_symbol, mouse_logFC = logFC, mouse_t = t)]
  )
  curtest <- cor.test(cur3$human_logFC, cur3$mouse_logFC, method = "spearman")
  add <- data.table(
    cell        = i,
    ngene       = nrow(cur3),
    human_ngene = nrow(cur1),
    mouse_ngene = nrow(cur2),
    rho         = curtest$estimate,
    p           = curtest$p.value
  )
  human_mouse_correlations <- rbind(human_mouse_correlations, add)
  human_mouse_plot_data    <- rbind(human_mouse_plot_data, data.table(cell = i, cur3))
}
human_mouse_correlations[, padj := p.adjust(p, "fdr")]
human_mouse_correlations

human_mouse_deg_overlap <- c()
for (i in MOUSE_CELL_TYPES) {
  cur1 <- copy(mydat[["human"]][[i]])
  cur2 <- copy(mydat[["mouse"]][[i]])
  hupg <- cur1[adj.P.Val < 0.05 & logFC > 0]$V1
  hdwg <- cur1[adj.P.Val < 0.05 & logFC < 0]$V1
  mupg <- cur2[adj.P.Val < 0.05 & logFC > 0]$human_symbol
  mdwg <- cur2[adj.P.Val < 0.05 & logFC < 0]$human_symbol
  cur1[, hup := FALSE][, hdw := FALSE][, mup := FALSE][, mdw := FALSE]
  cur1[adj.P.Val < 0.05 & logFC > 0, hup := TRUE]
  cur1[adj.P.Val < 0.05 & logFC < 0, hdw := TRUE]
  cur1[V1 %in% mupg, mup := TRUE]
  cur1[V1 %in% mdwg, mdw := TRUE]
  cur1[, hup := factor(hup, levels = c("TRUE", "FALSE"))]
  cur1[, hdw := factor(hdw, levels = c("TRUE", "FALSE"))]
  cur1[, mup := factor(mup, levels = c("TRUE", "FALSE"))]
  cur1[, mdw := factor(mdw, levels = c("TRUE", "FALSE"))]
  ftres1 <- fisher.test(table(cur1$hup, cur1$mup))
  ftres2 <- fisher.test(table(cur1$hdw, cur1$mdw))
  add1 <- data.table(cell = i, test = "up", h = length(hupg), m = length(mupg),
                     hm = length(intersect(hupg, mupg)), bg = nrow(cur1),
                     or = ftres1$estimate, p = ftres1$p.value)
  add2 <- data.table(cell = i, test = "dw", h = length(hdwg), m = length(mdwg),
                     hm = length(intersect(hdwg, mdwg)), bg = nrow(cur1),
                     or = ftres2$estimate, p = ftres2$p.value)
  human_mouse_deg_overlap <- rbind(human_mouse_deg_overlap, add1, add2)
}
human_mouse_deg_overlap_filtered <- human_mouse_deg_overlap[h >= 5 & m >= 5]
human_mouse_deg_overlap_filtered[, padj := p.adjust(p, "fdr")]

# ── Plot: human vs mouse logFC ────────────────────────────────────────────────
ggplot(human_mouse_plot_data, aes(human_logFC, mouse_logFC)) +
  geom_point(alpha = 0.3) +
  geom_hline(yintercept = 0, col = "red", lty = "dotted") +
  geom_vline(xintercept = 0, col = "red", lty = "dotted") +
  facet_wrap(~cell, scale = "free")

# ── Replication: CD4 vs PMC7427333 ───────────────────────────────────────────
rep <- as.data.table(read_excel(check_required_file(REP_XLSX)))
cd4_replication_data <- merge(mydat[["human"]][["cd4"]], as.data.table(rep)[, .(V1 = Gene, rep_logFC = avg_logFC)])
plot(cd4_replication_data$logFC, cd4_replication_data$rep_logFC)
cor.test(cd4_replication_data$logFC, cd4_replication_data$rep_logFC, method = "spearman")$estimate
cor.test(cd4_replication_data$logFC, cd4_replication_data$rep_logFC, method = "spearman")$p.value

replication_correlations <- c()
for (i in MOUSE_CELL_TYPES) {
  cur     <- copy(mydat[["human"]][[i]])
  cur3    <- merge(cur, as.data.table(rep)[, .(V1 = Gene, rep_logFC = avg_logFC)])
  curtest <- cor.test(cur3$logFC, cur3$rep_logFC, method = "spearman")
  add     <- data.table(cell = i, ngene = nrow(cur3), ngene1 = nrow(cur),
                        ngene2 = nrow(rep), rho = curtest$estimate, p = curtest$p.value)
  replication_correlations <- rbind(replication_correlations, add)
}

# ── Cross-cell-type correlation (human) ──────────────────────────────────────
cross_celltype_correlations <- c()
checked <- c()
for (i in HUMAN_CELL_TYPES) {
  cur1 <- mydat[["human"]][[i]]
  for (j in HUMAN_CELL_TYPES) {
    if (i != j & !paste(i, j) %in% checked & !paste(j, i) %in% checked) {
      checked <- c(checked, paste(i, j))
      cur2 <- mydat[["human"]][[j]]
      cur3 <- merge(cur1, cur2, by = "V1")
      add  <- data.table(
        cell1 = i, cell2 = j, ngene = nrow(cur3),
        rho   = cor(cur3$logFC.x, cur3$logFC.y, method = "spearman"),
        p     = cor.test(cur3$logFC.x, cur3$logFC.y, method = "spearman")$p.value
      )
      cross_celltype_correlations <- rbind(cross_celltype_correlations, add)
    }
  }
}
cross_celltype_correlations[, padj := p.adjust(p, "fdr")]
cross_celltype_correlations


# ── Human cohort 1 vs bulk ────────────────────────────────────────────────────
blk <- as.data.table(readRDS(check_required_file(BULK_RDS))$de)
check_required_cols(blk, c("symbol", "logFC"), "bulk differential expression")
bulk_replication_correlations <- c()
for (i in HUMAN_CELL_TYPES) {
  cur    <- merge(blk, mydat[["human"]][[i]], by.x = "symbol", by.y = "V1")
  curres <- cor.test(cur$logFC.x, cur$logFC.y, method = "spearman")
  add    <- data.table(cell = i, ngene = nrow(cur), rho = curres$estimate, p = curres$p.value)
  bulk_replication_correlations <- rbind(bulk_replication_correlations, add)
}
bulk_replication_correlations[, padj := p.adjust(p, "fdr")]
bulk_replication_correlations


# ── Human cohort 1 vs cohort 2 ────────────────────────────────────────────────
cohort2_correlations <- c()
for (i in HUMAN_CELL_TYPES) {
  cur1 <- mydat[["human"]][[i]]
  for (j in c("microglia", "astrocytes")) {
    curfile <- file.path(COHORT2_DIR, paste0("de_", j, "_", i, ".csv"))
    if (file.exists(curfile)) {
      cur2   <- fread(curfile)
      cur3   <- merge(cur1, cur2, by = "V1")
      curres <- cor.test(cur3$logFC.x, cur3$logFC.y, method = "spearman")
      add    <- data.table(batch1_cell = i, batch2_bloodcell = i, batch2_braincell = j,
                           ngene = nrow(cur3), rho = curres$estimate, p = curres$p.value)
      cohort2_correlations <- rbind(cohort2_correlations, add)
    }
  }
}
cohort2_correlations <- cohort2_correlations[order(rho, decreasing = TRUE)]
cohort2_correlations[, padj := p.adjust(p, "fdr")]
cohort2_correlations

# ── Human cohort 1 vs blood-CSF pairs (PMID31937773) ─────────────────────────
csf    <- list()
csf$up <- as.data.table(read_excel(check_required_file(CSF_XLSX), sheet = 1))
csf$dw <- as.data.table(read_excel(check_required_file(CSF_XLSX), sheet = 2))

mycsf <- rbind(
  csf$up[clusters == "B2",    .(symbol = `...1`, cell = "B2",    mycell = "bc",    dir = "csf")],
  csf$up[clusters == "CD4",   .(symbol = `...1`, cell = "CD4",   mycell = "cd4",  dir = "csf")],
  csf$up[clusters == "CD8a",  .(symbol = `...1`, cell = "CD8a",  mycell = "cd8",  dir = "csf")],
  csf$up[clusters == "CD8na", .(symbol = `...1`, cell = "CD8na", mycell = "cd8",  dir = "csf")],
  csf$up[clusters == "Mono2", .(symbol = `...1`, cell = "Mono2", mycell = "cd14", dir = "csf")],
  csf$up[clusters == "Mono2", .(symbol = `...1`, cell = "Mono2", mycell = "cd16", dir = "csf")],
  csf$up[clusters == "NK1",   .(symbol = `...1`, cell = "NK1",   mycell = "nk",   dir = "csf")],
  csf$up[clusters == "NK2",   .(symbol = `...1`, cell = "NK2",   mycell = "nk",   dir = "csf")],
  csf$up[clusters == "Tdg",   .(symbol = `...1`, cell = "Tdg",   mycell = "cd4",  dir = "csf")],
  csf$up[clusters == "Tdg",   .(symbol = `...1`, cell = "Tdg",   mycell = "cd8",  dir = "csf")],
  csf$up[clusters == "Tregs", .(symbol = `...1`, cell = "Tregs", mycell = "cd4",  dir = "csf")],
  csf$up[clusters == "Tregs", .(symbol = `...1`, cell = "Tregs", mycell = "cd8",  dir = "csf")],
  csf$dw[clusters == "CD4",   .(symbol = `...1`, cell = "CD4",   mycell = "cd4",  dir = "blood")],
  csf$dw[clusters == "CD8a",  .(symbol = `...1`, cell = "CD8a",  mycell = "cd8",  dir = "blood")],
  csf$dw[clusters == "CD8na", .(symbol = `...1`, cell = "CD8na", mycell = "cd8",  dir = "blood")],
  csf$dw[clusters == "Mono2", .(symbol = `...1`, cell = "Mono2", mycell = "cd14", dir = "blood")],
  csf$dw[clusters == "Mono2", .(symbol = `...1`, cell = "Mono2", mycell = "cd16", dir = "blood")],
  csf$dw[clusters == "NK1",   .(symbol = `...1`, cell = "NK1",   mycell = "nk",   dir = "blood")],
  csf$dw[clusters == "NK2",   .(symbol = `...1`, cell = "NK2",   mycell = "nk",   dir = "blood")],
  csf$dw[clusters == "Tdg",   .(symbol = `...1`, cell = "Tdg",   mycell = "cd4",  dir = "blood")],
  csf$dw[clusters == "Tdg",   .(symbol = `...1`, cell = "Tdg",   mycell = "cd8",  dir = "blood")],
  csf$dw[clusters == "Tregs", .(symbol = `...1`, cell = "Tregs", mycell = "cd4",  dir = "blood")],
  csf$dw[clusters == "Tregs", .(symbol = `...1`, cell = "Tregs", mycell = "cd8",  dir = "blood")]
)

iter   <- unique(mycsf[, .(cell, mycell, dir)])
csf_signature_enrichment <- c()
for (i in 1:nrow(iter)) {
  them   <- iter[i]$cell
  us     <- iter[i]$mycell
  direc  <- iter[i]$dir
  curdat <- copy(mydat[["human"]][[us]])
  curgen <- mycsf[cell == them & mycell == us & dir == direc]$symbol
  curdat[, d1 := FALSE][, d2 := FALSE]
  curdat[V1 %in% curgen, d2 := TRUE]
  if (direc == "csf") {
    curdat[adj.P.Val < 0.05 & logFC > 0, d1 := TRUE]
  } else {
    curdat[adj.P.Val < 0.05 & logFC < 0, d1 := TRUE]
  }
  curdat[, d1 := factor(d1, levels = c("TRUE", "FALSE"))]
  curdat[, d2 := factor(d2, levels = c("TRUE", "FALSE"))]
  if (nrow(curdat[d1 == TRUE]) >= 10 & nrow(curdat[d2 == TRUE]) >= 10) {
    ftres  <- fisher.test(table(curdat$d1, curdat$d2))
    add    <- data.table(cell = them, mycell = us, dir = direc, or = ftres$estimate, p = ftres$p.value)
    csf_signature_enrichment <- rbind(csf_signature_enrichment, add)
  }
}
csf_signature_enrichment[, padj := p.adjust(p, "fdr")]
