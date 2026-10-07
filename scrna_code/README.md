# APLNR knockdown in 2D human gastruloids — single-cell analysis code

Code used for the single-cell RNA-seq analyses and figures in:

> _[Authors]. [Title]. [Journal] ([Year]). DOI: [to be added]_

This repository documents how the analyses were done. It is not a software package.
Notebooks were run interactively; their outputs have been cleared to keep the repository
small.

## Data

Raw and processed data are available at **[GEO accession / Zenodo DOI — to be added]**.
External references used:

- IH human gastruloid timeseries (label-transfer reference) — [citation]
- Tyser et al. 2021, CS7 human gastrula (Figure 2 embryo panels) — [citation]

The code expects a data folder with this layout:

```
<data folder>/
  00_DATA_FROZEN/
    00_raw/amanda_liang_raw_unfiltered.h5ad
    01_source/sample_barcodes.csv            # final cell filter (barcode list)
    01_source/models/                        # trained scVI reference model
    02_annotated/adata_full_joint_OCM_labels_corrected.h5ad
    03_analysis_ready/adata_master_normalized[_magic].h5ad
  01_references/                             # IH timeseries, Tyser embryo
```

Tell the code where that folder is by setting an environment variable before running:

```bash
export SC_WORKFLOW_ROOT=/path/to/data_folder
```

All paths are resolved in `utils/config.py`. Outputs are written into the data folder
(`06_tables/`, `07_figures/`, `pseudobulk_ALL_<date>/`) or to `04_supplement/output/`.

## Repository layout and run order

### `01_preprocessing/`
| file | what it does |
|---|---|
| `2D_gastruloid_analysis_Amanda_JOINT_SCVI_newmodel.ipynb` | Loads the eight Cell Ranger matrices (2 sites × Ctrl/KD × 42 h/72 h), QC filtering, joint scVI integration with the IH gastruloid reference, kNN label transfer of cell states. |
| `07_canonical_master_build.ipynb` | Normalisation (`normalize_total` 1e4 + `log1p`), gene-symbol harmonisation, timepoint label correction → master object. |
| `08_figure2_panels.ipynb` | MAGIC imputation of APLNR and APELA only (used for UMAP colouring; never enters a statistical test). |

### `02_analysis/`
| file | what it does |
|---|---|
| `run_pseudobulk_ALL.py` | Pseudobulk DESeq2 (shAPLNR vs control, `~ site + condition`) for every cell state × timepoint. Source table for the KD figure. |
| `kd_pseudobulk.py` | Pseudobulk DE plus the pre-specified gene-panel test (BH within panel). |
| `kd_aplnr_correlation.py` | Genes co-varying with APLNR in control cells (input gene set for set enrichment). |
| `kd_set_enrichment.py`, `kd_set_enrichment_panels.py` | Expression-matched permutation tests of whether gene sets shift in the KD; panels. |
| `export_control_expression.py` | Control-only pseudobulk expression levels across states and timepoints. |
| `APLNR_module_preservation_FIXED.R` | Formal WGCNA module preservation (Ctrl → KD) of the APLNR-associated module. Reads the expression/module CSVs exported by `FIGURE2_current.ipynb`. |

### `03_figures/`
| file | what it does |
|---|---|
| `FIGURE2_current.ipynb` | Figure 2: control WGCNA module, Tyser embryo panels, control gastruloid UMAPs, dot plots, MAGIC panels, composition. |
| `FIGURE_KD_current.ipynb` | KD figure: knockdown validation, composition, UMAPs, pseudobulk heatmaps, module scores, GSEA. |

### `04_supplement/`
Run order: `qc_supplement.py` → `ih_validation.py` → `supp_figures.py`.
`supp_figures.py` assembles the supplementary figures (cell filtering, replicate
concordance, label validation against IH, KD controls). See `04_supplement/README.md`.

### `utils/`
`config.py` (paths), `stan_filter.py` (the barcode cell filter, applied everywhere),
`master_io.py` (reads the master object with the filter on by default), `palettes.py`
(figure colours and state groupings).

## Key analysis decisions

- **Cell filter.** A single barcode list (`sample_barcodes.csv`, 29,419 cells of 45,435
  post-QC) is applied at runtime by `utils/stan_filter.py` before any per-group statistic.
- **Replication.** Chicago and Michigan are two independent biological repeats; pseudobulk
  models use `~ site + condition`. Formal p-values are reported only for the pre-specified
  gene panel; genome-wide results are treated at the set level.
- **MAGIC** values are used for visualisation only.

## Environment

See `requirements.txt` (Python) — preprocessing used Python 3.12 / scanpy 1.12.2 /
anndata 0.13.2 / scvi-tools 1.5.0; downstream analyses used Python 3.11. The module
preservation script needs R with the `WGCNA` package.

## License

[To be chosen — e.g. MIT]

## Contact

Amanda Brooke — amandabrooke@uchicago.edu
