#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript scripts/deseq2_analysis.R config.json")

required_cran <- c("jsonlite", "ggplot2", "pheatmap", "ggrepel")
required_bioc <- c("tximport", "DESeq2", "limma")
missing_pkgs <- c(required_cran[!vapply(required_cran, requireNamespace, logical(1), quietly=TRUE)],
                  required_bioc[!vapply(required_bioc, requireNamespace, logical(1), quietly=TRUE)])
if (length(missing_pkgs)) {
  stop("Missing R packages: ", paste(missing_pkgs, collapse=", "),
       ". Run: Rscript scripts/install_R_packages.R")
}

suppressPackageStartupMessages({
  library(jsonlite)
  library(tximport)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(ggrepel)
  library(limma)
})

config_path <- normalizePath(args[1], winslash="/", mustWork=TRUE)
config_dir <- dirname(config_path)
cfg <- fromJSON(config_path, simplifyVector=TRUE)
project_root <- normalizePath(file.path(config_dir, cfg$project_root), winslash="/", mustWork=FALSE)
resolve_path <- function(x) {
  if (grepl("^([A-Za-z]:[/\\\\]|/)", x)) normalizePath(x, winslash="/", mustWork=FALSE)
  else normalizePath(file.path(project_root, x), winslash="/", mustWork=FALSE)
}

manifest_path <- resolve_path(cfg$manifest)
contrast_path <- resolve_path(cfg$contrasts)
tx2gene_path <- resolve_path(cfg$tx2gene)
quant_root <- resolve_path(cfg$kallisto_output_dir)
out_root <- resolve_path(cfg$deseq2_output_dir)
sensitivity_root <- resolve_path(cfg$sensitivity_output_dir)
expr_root <- resolve_path(cfg$expression_output_dir)
plot_root <- resolve_path(cfg$plot_output_dir)
dir.create(out_root, recursive=TRUE, showWarnings=FALSE)
dir.create(sensitivity_root, recursive=TRUE, showWarnings=FALSE)
dir.create(expr_root, recursive=TRUE, showWarnings=FALSE)
dir.create(plot_root, recursive=TRUE, showWarnings=FALSE)

meta <- read.delim(manifest_path, check.names=FALSE, stringsAsFactors=FALSE)
meta <- meta[tolower(as.character(meta$include)) %in% c("true","1","yes","y"), , drop=FALSE]
if (anyDuplicated(meta$sample_id)) stop("sample_id is not unique")
meta$celltype <- factor(meta$celltype, levels=c("MPS","APS"))
meta$condition <- factor(meta$condition, levels=c("CTRL","DOX"))
meta$day <- factor(meta$day, levels=c("D1","D2"))
meta$replicate <- factor(meta$replicate, levels=c("R1","R2","R3"))
meta$group <- factor(paste(meta$celltype, meta$day, meta$condition, sep="_"))
rownames(meta) <- meta$sample_id

files <- file.path(quant_root, meta$sample_id, "abundance.tsv")
names(files) <- meta$sample_id
if (!all(file.exists(files))) {
  stop("Missing abundance.tsv files: ", paste(names(files)[!file.exists(files)], collapse=", "))
}

txmap <- read.delim(gzfile(tx2gene_path), stringsAsFactors=FALSE, check.names=FALSE)
if (!all(c("target_id","gene_id","gene_name") %in% names(txmap))) stop("tx2gene must contain target_id, gene_id, gene_name")
tx2gene <- unique(txmap[,c("target_id","gene_id")])

message("Importing ", length(files), " kallisto abundance tables with tximport...")
# Supply an explicit base-R importer. This avoids a tximport/readr importer
# lookup problem seen in some Windows R installations.
kallisto_tsv_importer <- function(path) {
  utils::read.delim(path, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                    check.names=FALSE, comment.char="", quote="")
}
txi <- tximport(files, type="none", txIn=TRUE, txOut=FALSE,
                tx2gene=tx2gene,
                txIdCol="target_id", abundanceCol="tpm",
                countsCol="est_counts", lengthCol="eff_length",
                importer=kallisto_tsv_importer)

# Preserve Ensembl IDs as the unique row identifier; symbols are annotation only.
gene_annot <- txmap[,c("gene_id","gene_name")]
gene_annot <- gene_annot[!duplicated(gene_annot$gene_id),]
rownames(gene_annot) <- gene_annot$gene_id

# One balanced model handles all simple pairwise comparisons between the 8 biological groups,
# while replicate blocks R1/R2/R3 are explicitly adjusted for.
dds <- DESeqDataSetFromTximport(txi, colData=meta, design=~ replicate + group)
keep <- rowSums(counts(dds) >= 10) >= 4
dds <- dds[keep,]
message("Genes retained after prefilter: ", sum(keep), " / ", length(keep))
dds <- DESeq(dds)
png("results/deseq2/DESeq2_dispersion_plot.png",
    width = 1800, height = 1500, res = 200)
plotDispEsts(dds)
dev.off()

saveRDS(dds, "results/deseq2/dds_full_model.rds")


vsd <- vst(dds, blind=FALSE)

saveRDS(dds, file.path(expr_root, "dds_group_model.rds"))
saveRDS(vsd, file.path(expr_root, "vst_group_model.rds"))
write.table(meta, file.path(expr_root, "metadata_used.tsv"), sep="\t", quote=FALSE, col.names=NA)

annotate_matrix <- function(mat) {
  ids <- rownames(mat)
  symbols <- gene_annot[ids,"gene_name"]
  data.frame(gene_id=ids, gene_name=symbols, mat, check.names=FALSE)
}
write.csv(annotate_matrix(counts(dds, normalized=TRUE)), gzfile(file.path(expr_root, "normalized_counts.csv.gz")), row.names=FALSE)
write.csv(annotate_matrix(assay(vsd)), gzfile(file.path(expr_root, "vst_counts.csv.gz")), row.names=FALSE)
write.csv(annotate_matrix(txi$abundance[rownames(dds),,drop=FALSE]), gzfile(file.path(expr_root, "gene_tpm.csv.gz")), row.names=FALSE)

# Genes of interest export
interest <- unlist(cfg$genes_of_interest)
interest_ids <- rownames(gene_annot)[gene_annot$gene_name %in% interest]
norm <- counts(dds, normalized=TRUE)
if (length(interest_ids)) {
  goi <- as.data.frame(as.table(norm[intersect(interest_ids, rownames(norm)),,drop=FALSE]))
  names(goi) <- c("gene_id","sample_id","normalized_count")
  goi$gene_name <- gene_annot[as.character(goi$gene_id),"gene_name"]
  meta_export <- meta[, setdiff(names(meta), "sample_id"), drop=FALSE]
  meta_export$sample_id <- rownames(meta_export)
  goi <- merge(goi, meta_export, by="sample_id", all.x=TRUE)
  write.csv(goi, file.path(expr_root, "genes_of_interest_normalized_counts.csv"), row.names=FALSE)
}

# Global QC plots
pca_data <- plotPCA(vsd, intgroup=c("celltype","condition","day","replicate"), returnData=TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))
p1 <- ggplot(pca_data, aes(PC1, PC2, color=celltype, shape=condition)) +
  geom_point(size=3) + facet_wrap(~day) + theme_classic(base_size=12) +
  labs(title="PCA: VST counts (replicate included in DESeq2 design)",
       x=paste0("PC1: ", percent_var[1], "%"), y=paste0("PC2: ", percent_var[2], "%"))
ggsave(file.path(plot_root,"PCA_VST_raw.png"), p1, width=8, height=5, dpi=300)
write.csv(pca_data, file.path(expr_root,"PCA_VST_raw_coordinates.csv"))

# Replicate-adjusted values are for visualization only, never for DE testing.
mm <- model.matrix(~ celltype * day * condition, data=meta)
adj <- removeBatchEffect(assay(vsd), batch=meta$replicate, design=mm)
pca <- prcomp(t(adj), scale.=FALSE)
meta_plot <- meta[, setdiff(names(meta), "sample_id"), drop=FALSE]
adj_df <- data.frame(sample_id=rownames(pca$x), PC1=pca$x[,1], PC2=pca$x[,2], meta_plot[rownames(pca$x),], check.names=FALSE)
adj_var <- 100 * (pca$sdev^2 / sum(pca$sdev^2))
p2 <- ggplot(adj_df, aes(PC1, PC2, color=celltype, shape=condition)) +
  geom_point(size=3) + facet_wrap(~day) + theme_classic(base_size=12) +
  labs(title="PCA: replicate-adjusted VST (visualization only)",
       x=sprintf("PC1: %.1f%%", adj_var[1]), y=sprintf("PC2: %.1f%%", adj_var[2]))
ggsave(file.path(plot_root,"PCA_VST_replicate_adjusted.png"), p2, width=8, height=5, dpi=300)
write.csv(adj_df, file.path(expr_root,"PCA_VST_replicate_adjusted_coordinates.csv"), row.names=FALSE)

png(file.path(plot_root,"sample_correlation_VST.png"), width=2400, height=2200, res=250)
pheatmap(cor(assay(vsd), method="pearson"), annotation_col=as.data.frame(meta[,c("celltype","condition","day","replicate")]), main="Sample correlation (VST)")
dev.off()

lib_df <- data.frame(sample_id=colnames(dds), library_size=colSums(counts(dds)), meta_plot[colnames(dds),], check.names=FALSE)
p_lib <- ggplot(lib_df, aes(x=reorder(sample_id, library_size), y=library_size, fill=celltype)) +
  geom_col() + coord_flip() + theme_classic(base_size=10) + labs(x=NULL, y="Imported estimated counts", title="Library size after tximport")
ggsave(file.path(plot_root,"library_sizes.png"), p_lib, width=8, height=8, dpi=300)
write.csv(lib_df, file.path(expr_root,"library_sizes.csv"), row.names=FALSE)

safe_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)
add_annotation <- function(res_df) {
  res_df$gene_id <- rownames(res_df)
  res_df$gene_name <- gene_annot[res_df$gene_id,"gene_name"]
  res_df <- res_df[,c("gene_id","gene_name",setdiff(names(res_df),c("gene_id","gene_name")))]
  res_df
}

plot_contrast <- function(res_df, contrast_name, dds_obj, numerator_samples=NULL, denominator_samples=NULL) {
  clean <- res_df
  clean$padj_plot <- ifelse(is.na(clean$padj), 1, pmax(clean$padj, 1e-300))
  clean$status <- "not significant"
  clean$status[!is.na(clean$padj) & clean$padj < 0.05 & clean$log2FoldChange >= 1] <- "up"
  clean$status[!is.na(clean$padj) & clean$padj < 0.05 & clean$log2FoldChange <= -1] <- "down"
  lab <- clean[order(clean$padj_plot),,drop=FALSE]
  lab <- lab[!is.na(lab$gene_name) & lab$gene_name != "",,drop=FALSE]
  lab <- head(lab, 12)
  p <- ggplot(clean, aes(log2FoldChange, -log10(padj_plot), color=status)) +
    geom_point(size=1, alpha=0.55) +
    geom_vline(xintercept=c(-1,1), linetype=2) + geom_hline(yintercept=-log10(0.05), linetype=2) +
    scale_color_manual(values=c("down"="#3B6FB6","not significant"="grey70","up"="#B54A4A")) +
    geom_text_repel(data=lab, aes(label=gene_name), size=3, max.overlaps=Inf) +
    theme_classic(base_size=11) + labs(title=contrast_name, color=NULL, x="log2 fold change", y="-log10 adjusted p-value")
  ggsave(file.path(plot_root,paste0("volcano__",safe_name(contrast_name),".png")), p, width=7, height=6, dpi=300)

  ma <- clean[is.finite(clean$baseMean) & is.finite(clean$log2FoldChange),,drop=FALSE]
  p_ma <- ggplot(ma, aes(log10(baseMean + 1), log2FoldChange, color=status)) +
    geom_point(size=1, alpha=0.55) + geom_hline(yintercept=0, linetype=2) +
    scale_color_manual(values=c("down"="#3B6FB6","not significant"="grey70","up"="#B54A4A")) +
    theme_classic(base_size=11) + labs(title=paste0("MA: ",contrast_name), color=NULL,
      x="log10(baseMean + 1)", y="log2 fold change")
  ggsave(file.path(plot_root,paste0("MA__",safe_name(contrast_name),".png")), p_ma, width=7, height=5.5, dpi=300)

  ranked <- clean[order(clean$padj_plot, clean$pvalue),,drop=FALSE]
  sig_ids <- ranked$gene_id[!is.na(ranked$padj) & ranked$padj < 0.05]
  selected <- head(if (length(sig_ids) >= 10) sig_ids else ranked$gene_id, 40)
  selected <- intersect(selected, rownames(vsd))
  if (length(selected) >= 2) {
    z <- assay(vsd)[selected,,drop=FALSE]
    z <- t(scale(t(z)))
    z[is.na(z)] <- 0
    labels <- gene_annot[selected,"gene_name"]
    labels[is.na(labels) | labels==""] <- selected[is.na(labels) | labels==""]
    rownames(z) <- make.unique(labels)
    ann <- as.data.frame(meta[,c("celltype","condition","day","replicate")])
    png(file.path(plot_root,paste0("heatmap__",safe_name(contrast_name),".png")), width=2200, height=2400, res=250)
    pheatmap(z, annotation_col=ann, main=paste0(contrast_name," — top 40"), show_colnames=TRUE, fontsize_row=7)
    dev.off()
  }
}

summary_rows <- list()
record_result <- function(res, contrast_name, interpretation, dds_obj=dds) {
  df <- add_annotation(as.data.frame(res))
  df <- df[order(df$padj, na.last=TRUE),]
  write.csv(df, file.path(out_root,paste0(safe_name(contrast_name),".csv")), row.names=FALSE)
  up <- sum(!is.na(df$padj) & df$padj < 0.05 & df$log2FoldChange >= 1)
  down <- sum(!is.na(df$padj) & df$padj < 0.05 & df$log2FoldChange <= -1)
  summary_rows[[length(summary_rows)+1]] <<- data.frame(contrast_name=contrast_name, interpretation=interpretation,
      tested_genes=nrow(df), FDR05_up_LFC1=up, FDR05_down_LFCminus1=down, stringsAsFactors=FALSE)
  plot_contrast(df, contrast_name, dds_obj)
  invisible(df)
}

# 12 primary pairwise comparisons from a single blocked group model.
main_results <- list()
contrast_table <- read.delim(contrast_path, stringsAsFactors=FALSE, check.names=FALSE)
for (i in seq_len(nrow(contrast_table))) {
  x <- contrast_table[i,]
  message("Contrast: ", x$contrast_name)
  res <- results(dds, contrast=c("group", x$numerator_group, x$denominator_group), alpha=0.05)
  main_results[[x$contrast_name]] <- record_result(res, x$contrast_name, x$interpretation, dds)
}


# Four difference-of-differences interaction tests.
run_interaction <- function(subset_idx, design_formula, factor1, factor2, contrast_name, interpretation) {
  sub_meta <- droplevels(meta[subset_idx,,drop=FALSE])
  sample_ids <- rownames(sub_meta)
  sub_txi <- txi
  for (nm in c("abundance","counts","length")) sub_txi[[nm]] <- txi[[nm]][,sample_ids,drop=FALSE]
  d <- DESeqDataSetFromTximport(sub_txi, colData=sub_meta, design=design_formula)
  d <- d[rowSums(counts(d) >= 10) >= 4,]
  d <- DESeq(d)
  rn <- resultsNames(d)
  candidates <- rn[grepl(factor1, rn) & grepl(factor2, rn)]
  if (length(candidates) != 1) stop("Could not uniquely identify interaction coefficient for ", contrast_name, ". Results names: ", paste(rn,collapse=", "))
  res <- results(d, name=candidates[1], alpha=0.05)
  record_result(res, contrast_name, interpretation, d)
}

for (dlevel in c("D1","D2")) {
  idx <- meta$day == dlevel
  run_interaction(idx, ~ replicate + celltype + condition + celltype:condition,
                  "celltype","condition",
                  paste0("interaction_DOX_effect__APS_vs_MPS_",dlevel),
                  paste0("Difference of DOX effects: (APS DOX-CTRL) minus (MPS DOX-CTRL) on ",dlevel))
}
for (ct in c("MPS","APS")) {
  idx <- meta$celltype == ct
  run_interaction(idx, ~ replicate + day + condition + day:condition,
                  "day","condition",
                  paste0("interaction_DOX_effect__D2_vs_D1_",ct),
                  paste0("Difference of DOX effects: (D2 DOX-CTRL) minus (D1 DOX-CTRL) within ",ct))
}

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(out_root,"contrast_summary.csv"), row.names=FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_root,"R_sessionInfo.txt"))
message("Finished. Differential-expression outputs: ", out_root)
message("QC and result plots: ", plot_root)
