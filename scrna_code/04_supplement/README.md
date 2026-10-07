# Supplementary QC and validation figures

| file | what it does |
|---|---|
| `qc_style.py` | Shared figure style (fonts, sizes, colormap). |
| `qc_supplement.py` | QC measurements: filtering thresholds, batch structure, per-sample library-size gate, depth-standardised composition, label-transfer confidence. Writes one panel per question plus `QC_SOURCE_DATA.xlsx`. |
| `ih_validation.py` | Label validation: for each of our cells, the labels of its 25 nearest IH reference cells in the joint scVI latent space. Writes `IH_VALIDATION_SOURCE_DATA.xlsx`. |
| `supp_figures.py` | Makes no measurements; imports the two scripts above and lays out the supplementary figures. |

Run order: `qc_supplement.py` → `ih_validation.py` → `supp_figures.py`.
Outputs go to `./output/` (PNG, PDF, SVG); expensive intermediates are cached in `./cache/`.

Figures produced by `supp_figures.py`:

| figure | content |
|---|---|
| `SF2_1_cell_filtering` | genes / UMIs / mito with thresholds, library-size gate, cells per sample at each stage |
| `SF2_2_replicate_concordance` | UMAP and state composition per site |
| `SF2_3_label_validation` | neighbourhood concordance with IH labels; eNasM/lNasM breakdown |
| `SF5_1_KD_controls` | condition UMAP, depth standardisation, label-confidence sensitivity |

QC thresholds (min 500 genes, min 1,000 UMIs, max 20% mitochondrial, per-sample upper cut
at median + 5 × MAD on genes and counts) are those used in
`01_preprocessing/2D_gastruloid_analysis_Amanda_JOINT_SCVI_newmodel.ipynb`. No separate
doublet caller was run; the per-sample upper MAD cut serves as the doublet control.
