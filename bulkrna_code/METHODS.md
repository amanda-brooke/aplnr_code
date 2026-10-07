# Methods — bulk RNA-seq (n = 3 biological replicates)

*Matches the LaTeX methods pasted into the manuscript. Items in [brackets] need confirming/filling.*

## Samples and design

Paired-end bulk RNA-seq of hESC-derived anterior (APS) and mid (MPS)
primitive-streak populations at two differentiation days (24 and 48 h), each with
a doxycycline-inducible **APLNR knockdown (shAPLNR)** and a matched control, in
**n = 3 biological replicates** per condition (**R1–R3; 24 libraries**).

## Sequencing QC

Per-sample read quality was assessed with **FastQC** and aggregated with
**MultiQC** (sequencing core; reports in `results/qc/`), including an rRNA-content
check. For QC, reads were mapped to **hg38** — genome mapping was used for QC only;
transcript quantification used kallisto against the GENCODE transcriptome (below).

## Quantification

Reads were pseudoaligned and transcript abundances estimated with
**kallisto v0.46.1** against the **GENCODE v46 (GRCh38)** transcriptome
(`gencode.v46.transcripts.fa.gz`; 254,070 transcript targets); libraries were
**unstranded** (median pseudoalignment rate ≈ 87% [confirm]). Transcript-level
estimates were summarized to genes with **tximport v1.30.0** using the GENCODE v46
transcript-to-gene map.

## Differential expression

**DESeq2 v1.42.1** (R v4.3.3). Counts imported via `DESeqDataSetFromTximport`
with design **`~ replicate + group`** (`group` = celltype × day × condition);
genes pre-filtered to **≥ 10 counts in ≥ 4 samples**. The variance-stabilizing
transform (`vst`, `blind = FALSE`) was used for visualization and PCA; for PCA,
replicate batch was removed with limma `removeBatchEffect` under a
`~ celltype * day * condition` model. Contrasts:

- **Knockdown effect** — shAPLNR vs control within each celltype × day
- **Temporal** — 48 vs 24 h within each celltype × condition
- **Celltype** — APS vs MPS within each condition × day
- **Interactions** — difference of the knockdown effect across celltype or day

Differentially expressed at **padj < 0.05 and |log2FC| ≥ 1** (Benjamini–Hochberg).

## Gene-set enrichment

Preranked GSEA with **GSEApy** (`prerank`), ranking genes by the DESeq2 Wald
statistic, against **MSigDB v2026.1**: **Hallmark**, **GO Biological Process**
(`c5.go.bp`) and **Reactome** (`c2.cp`) collections (gene symbols). Gene sets
restricted to **15–500 genes**; **1000 permutations**; **seed 149**; weighted
(weight 1.0). Terms significant at **FDR < 0.25**. Developmental, patterning,
migration, adhesion and matrix terms were highlighted from the GO:BP results.

## Figure-level analysis (python)

Figures generated in python (numpy, pandas, matplotlib, scipy) from the
DESeq2/expression/GSEA outputs. QC for symbol-level and unbiased analyses:

- Duplicate gene symbols (multiple Ensembl IDs) collapsed to the
  **highest-expressed** isoform.
- **Ribosomal-protein, mitochondrial and rRNA** genes excluded from unbiased
  variable-gene selection (heatmaps).

**Trajectory attenuation:** for genes with a significant temporal program in
control, the per-gene 48-vs-24 h log2FC under shAPLNR was compared with control.
Concordance is reported as Pearson r; the magnitude (slope) is estimated by
**Deming regression** (errors-in-both-variables), with the error-variance ratio
set from the DESeq2 standard errors and a **95% CI from 500 gene-level
bootstraps**.

## Data / code availability

This repository contains the analysis code, sample metadata, and the figure/GSEA
notebooks. Raw reads and processed expression matrices are deposited at
[GEO/SRA accession ??]. Reference files (GENCODE v46 transcriptome, MSigDB v2026.1
gene sets) are obtained from their original sources.
