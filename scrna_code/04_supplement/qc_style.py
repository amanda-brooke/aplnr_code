"""
qc_style.py -- ONE place for figure style in the QC supplement.

Same reasoning as palettes.py: style that lives inside each script drifts, and
then panels that sit side by side in one supplementary figure don't match. Both
qc_supplement.py and ih_validation.py import this and call apply_style().

Sizes follow the figure2_clean notebooks (labels 11, titles 13). Arial is listed
first; Liberation Sans is metrically identical and is the fallback on machines
without Arial, so layout does not shift between them.

NOTE ON FONTS IN OUTPUT
  .svg  -- svg.fonttype "none" keeps text as <text> with font-family Arial first,
           so Illustrator picks up real Arial and sizes stay editable.
  .pdf  -- fonttype 42 EMBEDS whatever font actually rendered. Regenerate on a
           machine with Arial installed if the embedded PDF font matters.
"""
import matplotlib.pyplot as plt

FONT_STACK = ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"]

BASE, TITLE, TICK, LEGEND = 11, 13, 10, 10
SMALL, TINY = 9, 8          # in-panel annotations; use these, not raw numbers

RC = {
    "font.family": "sans-serif", "font.sans-serif": FONT_STACK,
    "font.size": BASE, "axes.labelsize": BASE, "axes.titlesize": TITLE,
    "xtick.labelsize": TICK, "ytick.labelsize": TICK, "legend.fontsize": LEGEND,
    "figure.titlesize": TITLE,
    "figure.dpi": 120, "savefig.dpi": 300, "savefig.bbox": "tight",
    "axes.grid": False, "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.9, "xtick.major.width": 0.9, "ytick.major.width": 0.9,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
}

# Sequential colourmap for the concordance matrix.
MATRIX_CMAP = "magma_r"


def apply_style():
    plt.rcParams.update(RC)
