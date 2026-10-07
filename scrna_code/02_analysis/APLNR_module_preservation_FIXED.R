#!/usr/bin/env Rscript

# =============================================================================
# APLNR Ctrl -> KD formal WGCNA module preservation
#
# Put this script in the SAME folder as:
#   modulePreservation_Ctrl_expression.csv
#   modulePreservation_KD_expression.csv
#   Ctrl_WGCNA_gene_module_colors.csv
#
# Run from Terminal:
#   Rscript APLNR_module_preservation.R
#
# Outputs:
#   ALL_modules_formal_preservation.csv
#   APLNR_black_module_formal_preservation.csv
#   APLNR_modulePreservation.RData
#   APLNR_modulePreservation_summary.txt
# =============================================================================

options(stringsAsFactors = FALSE)

# Work from the folder containing this script.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) == 1) {
  script_path <- normalizePath(sub("^--file=", "", file_arg))
  setwd(dirname(script_path))
}
cat("\nWorking directory:\n", getwd(), "\n\n")

# Install/load WGCNA if needed.
if (!requireNamespace("WGCNA", quietly = TRUE)) {
  cat("WGCNA is not installed. Installing required packages...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install(
    c("impute", "preprocessCore", "GO.db", "AnnotationDbi"),
    ask = FALSE,
    update = FALSE
  )
  install.packages("WGCNA", repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages(library(WGCNA))
try(allowWGCNAThreads(), silent = TRUE)

CTRL_FILE <- "modulePreservation_Ctrl_expression.csv"
KD_FILE   <- "modulePreservation_KD_expression.csv"
MOD_FILE  <- "Ctrl_WGCNA_gene_module_colors.csv"

required <- c(CTRL_FILE, KD_FILE, MOD_FILE)
missing_files <- required[!file.exists(required)]
if (length(missing_files) > 0) {
  stop(paste0(
    "\nMissing required file(s):\n  ",
    paste(missing_files, collapse = "\n  "),
    "\n\nPut this R script in the same folder as the three CSV files."
  ))
}

# Read data. Expression CSV rows = pseudocells, columns = genes.
ctrl <- read.csv(CTRL_FILE, row.names = 1, check.names = FALSE)
kd   <- read.csv(KD_FILE,   row.names = 1, check.names = FALSE)
mods <- read.csv(MOD_FILE, check.names = FALSE)

if (!all(c("gene", "module") %in% colnames(mods))) {
  stop("Ctrl_WGCNA_gene_module_colors.csv must contain columns named 'gene' and 'module'.")
}

ctrl[] <- lapply(ctrl, as.numeric)
kd[]   <- lapply(kd, as.numeric)

# Match genes across Ctrl, KD, and module assignments.
common <- Reduce(intersect, list(colnames(ctrl), colnames(kd), as.character(mods$gene)))
if (length(common) < 100) {
  stop(paste("Only", length(common), "common genes found. Check the input files."))
}

common <- colnames(ctrl)[colnames(ctrl) %in% common]
ctrl <- ctrl[, common, drop = FALSE]
kd   <- kd[, common, drop = FALSE]
mods <- mods[match(common, mods$gene), , drop = FALSE]

stopifnot(
  identical(colnames(ctrl), colnames(kd)),
  identical(colnames(ctrl), as.character(mods$gene))
)

module_colors <- as.character(mods$module)
names(module_colors) <- as.character(mods$gene)

cat("Input summary\n-------------\n")
cat("Genes:            ", length(common), "\n")
cat("Ctrl pseudocells: ", nrow(ctrl), "\n")
cat("KD pseudocells:   ", nrow(kd), "\n\n")
cat("Ctrl module sizes:\n")
print(sort(table(module_colors), decreasing = TRUE))

if (!"black" %in% module_colors) stop("No 'black' module found.")
cat("\nBlack module genes:", sum(module_colors == "black"), "\n\n")

# QC.
if (any(!is.finite(as.matrix(ctrl)))) stop("Ctrl matrix contains NA/Inf values.")
if (any(!is.finite(as.matrix(kd))))   stop("KD matrix contains NA/Inf values.")

sd_ctrl <- apply(ctrl, 2, sd)
sd_kd   <- apply(kd, 2, sd)
keep <- is.finite(sd_ctrl) & is.finite(sd_kd) & sd_ctrl > 0 & sd_kd > 0

if (any(!keep)) {
  cat("Removing", sum(!keep), "genes with zero/non-finite variance in Ctrl or KD.\n")
  ctrl <- ctrl[, keep, drop = FALSE]
  kd   <- kd[, keep, drop = FALSE]
  module_colors <- module_colors[keep]
}

cat("Genes entering modulePreservation:", ncol(ctrl), "\n")
cat("Black genes entering test:", sum(module_colors == "black"), "\n\n")

# Ctrl = reference; KD = test. We supply module labels only for Ctrl.
multiData <- list(
  Ctrl = list(data = as.data.frame(ctrl)),
  KD   = list(data = as.data.frame(kd))
)
multiColor <- list(Ctrl = module_colors)

check <- checkSets(multiData, checkStructure = TRUE)
if (!isTRUE(check$structureOK)) {
  print(check)
  stop("checkSets() reported a problem with the input data.")
}

cat("checkSets OK:",
    check$nSets, "sets,",
    check$nGenes, "genes; samples =",
    paste(check$nSamples, collapse = ", "),
    "\n\n")

# Formal preservation test. Original network was signed.
set.seed(17)
cat("\nRunning WGCNA modulePreservation with 200 permutations...\n")
cat("This can take a while.\n\n")

mp <- modulePreservation(
  multiData = multiData,
  multiColor = multiColor,
  dataIsExpr = TRUE,
  networkType = "signed",
  corFnc = "cor",
  corOptions = "use = 'p'",
  referenceNetworks = 1,
  testNetworks = 2,
  nPermutations = 200,
  includekMEallInSummary = FALSE,
  randomSeed = 17,
  quickCor = 0,
  verbose = 3
)

save(mp, file = "APLNR_modulePreservation.RData")

# Extract Ctrl-reference -> KD-test preservation results.
obs <- as.data.frame(mp$preservation$observed[[1]][[2]])
z   <- as.data.frame(mp$preservation$Z[[1]][[2]])
obs$module <- rownames(obs)
z$module   <- rownames(z)

results <- merge(
  obs, z,
  by = "module",
  suffixes = c("_observed", "_Z"),
  all = TRUE,
  sort = FALSE
)

write.csv(results, "ALL_modules_formal_preservation.csv", row.names = FALSE)
black <- results[results$module == "black", , drop = FALSE]
write.csv(black, "APLNR_black_module_formal_preservation.csv", row.names = FALSE)

# Pull headline values robustly.
get_value <- function(df, row_name, possible_cols) {
  if (!row_name %in% rownames(df)) return(NA_real_)
  found <- possible_cols[possible_cols %in% colnames(df)]
  if (length(found) == 0) return(NA_real_)
  as.numeric(df[row_name, found[1]])
}

zsummary <- get_value(z, "black", c("Zsummary.pres", "Zsummary"))
median_rank <- get_value(obs, "black", c("medianRank.pres", "medianRank"))
module_size <- get_value(z, "black", c("moduleSize"))
if (is.na(module_size)) module_size <- sum(module_colors == "black")

interpretation <- if (is.na(zsummary)) {
  "Zsummary column was not found automatically; inspect the black-module CSV."
} else if (zsummary < 2) {
  "Little evidence of preservation."
} else if (zsummary <= 10) {
  "Moderate preservation."
} else {
  "Strong preservation."
}

summary_lines <- c(
  "APLNR Ctrl-defined black module -> KD formal WGCNA preservation",
  "===============================================================",
  paste("Genes tested:", ncol(ctrl)),
  paste("Ctrl pseudocells:", nrow(ctrl)),
  paste("KD pseudocells:", nrow(kd)),
  paste("Black module size:", module_size),
  paste("Zsummary.pres:", signif(zsummary, 5)),
  paste("medianRank.pres:", signif(median_rank, 5)),
  paste("Interpretation:", interpretation),
  "",
  "Rule of thumb for Zsummary.pres:",
  "  < 2  : little/no evidence of preservation",
  "  2-10 : moderate preservation",
  "  > 10 : strong preservation",
  "",
  "Lower medianRank.pres = more preserved relative to the other Ctrl modules."
)

writeLines(summary_lines, "APLNR_modulePreservation_summary.txt")

cat("\n\n", paste(summary_lines, collapse = "\n"), "\n\n", sep = "")
cat("Saved:\n")
cat("  ALL_modules_formal_preservation.csv\n")
cat("  APLNR_black_module_formal_preservation.csv\n")
cat("  APLNR_modulePreservation.RData\n")
cat("  APLNR_modulePreservation_summary.txt\n")
