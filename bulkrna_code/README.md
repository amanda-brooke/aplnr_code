# bulk_rna_pipeline

Bulk RNA-seq analysis pipeline and code for the associated manuscript.
Transcript quantification with **kallisto**, QC (FastQC/MultiQC), and
differential expression with **DESeq2**, plus downstream GSEA and gene-set
scoring. Built for an MPS/APS knockdown study, and parameterized so a lab
member can point it at their own samples by editing three files.

## Analyzing your own data

1. **`metadata/samples.tsv`** — one row per sample: `sample_id`, grouping
   columns (e.g. `condition`, `day`, `replicate`), and an `include` flag.
   Add columns pointing to your FASTQ files.
2. **`metadata/contrasts.tsv`** — the DESeq2 comparisons you want
   (numerator vs. denominator group + a plain-language interpretation).
3. **`scripts/config.json`** — paths (raw FASTQ dir, reference, outputs),
   thread count, strandedness, and genes of interest.

Then run the steps below.

## Pipeline steps

1. Install dependencies: `Rscript scripts/install_R_packages.R`, plus Python
   tools on PATH (`kallisto`, `fastqc`, `multiqc`) and Python packages
   (`pandas`, `gseapy` for GSEA, `scanpy` if using the scoring notebooks).
2. Build a reference: a transcriptome FASTA + tx2gene table (e.g. GENCODE) and
   a kallisto index. Point `scripts/config.json` at them.
3. Quantify + QC + assemble counts: `python scripts/run_pipeline.py`
4. Differential expression: `Rscript scripts/deseq2_analysis.R`
   (`scripts/_no_R2_sensitivity/` holds a leave-one-replicate-out variant).
5. Downstream analysis: run the notebooks in `notebook/` (GSEA, gene-set scoring).

## Repository layout

```
scripts/            Pipeline driver + DESeq2 code + config
metadata/           samples.tsv (design) + contrasts.tsv (comparisons)
notebook/           Downstream analysis: GSEA preranked, gene-set scoring
notes/              Analysis notes
METHODS.md          Methods write-up
```

## Data availability & citation

Raw sequencing reads and processed count matrices for the study are deposited
under GEO accession **[GSE_XXXXXX — fill in]**. Reference annotations are
GENCODE v46 and MSigDB (obtain from their sources).

If you use this pipeline or data, please cite **[manuscript citation — fill in]**.

## Not included (by design)

Raw FASTQs, reference/genome/MSigDB files (large; MSigDB is
redistribution-restricted), and derived results. Regenerate results by running
the pipeline; obtain data from GEO. These are excluded via `.gitignore`.
