"""
stan_filter.py -- ONE definition of Stan's barcode filter, imported everywhere.

WHY THIS EXISTS
The filter was applied in FIGURE2_CLEAN.ipynb, FIGURE2_SPATIAL.ipynb,
FIGURE_KD_FIGURE.ipynb and FIGURE_KD_CLEAN.ipynb -- and NOT in maturation_kd.ipynb,
kd_pseudobulk.py, chat_scores_maturation.py or kd_aplnr_correlation.py. So the figure
panels were drawn on 29,419 cells while every pseudobulk, score and correlation was
computed on 45,435. That is a 35% difference in the cell set between panels of the
same figure, and it is what made KD_pseudobulk_QC.csv report 53 NasM cells for
Chicago_Ctrl_72h while 2_composition.csv reported 36.

Same `predicted_leiden_sub_joint` labels in both. The discrepancy was ONLY the filter.

WHAT IT DOES
`sample_barcodes.csv` (29,519 rows) is matched against `obs_names` by pulling the first
run of >=14 ACGT characters out of each. Cells whose barcode is not in that set are
dropped. Retention is 64.7% overall and varies by group (50.6% to 71.9%), so this is
NOT a uniform thinning -- it changes group sizes unevenly and must be applied before
any per-group statistic, not after.

MIN_CELLS_PER_SAMPLE CONSEQUENCE
Filtering takes Chicago_Ctrl_72h NasM from 53 cells to 36, below the old floor of 50.
The floor is therefore 35, which keeps all four 72h NasM pseudobulks and preserves the
paired 2x2 design with 1 residual df. Do not raise it above 35 without re-checking
KD_pseudobulk_QC.csv -- at 40 the NasM fit loses Chicago and becomes unfittable.

USAGE
    from stan_filter import stan_mask
    keep = stan_mask(obs_names, ROOT)      # boolean array, same length as obs_names
    M, obs = M[keep], obs.loc[keep].reset_index(drop=True)
"""
import re
import numpy as np
import pandas as pd
from pathlib import Path

BARCODE_RE = re.compile(r"[ACGT]{14,}")
CSV_NAME = "sample_barcodes.csv"

# Canonical home is the frozen data folder. The project-root location is the
# pre-freeze path and is kept as a fallback so older checkouts still work.
CSV_CANDIDATES = (
    Path("00_DATA_FROZEN") / "01_source" / CSV_NAME,
    Path(CSV_NAME),
)


def parse_barcode(x):
    """First run of >=14 ACGT chars, or None. Matches FIGURE2_CLEAN's bc16()."""
    h = BARCODE_RE.findall(str(x))
    return h[0] if h else None


def stan_barcodes(root):
    p = next((Path(root) / c for c in CSV_CANDIDATES if (Path(root) / c).exists()), None)
    if p is None:
        tried = "\n  ".join(str(Path(root) / c) for c in CSV_CANDIDATES)
        raise FileNotFoundError(
            f"sample_barcodes.csv not found. Tried:\n  {tried}\n"
            f"Stan's barcode list is the single cell filter for this "
            f"project; without it the KD pipeline silently analyses 45,435 cells while "
            f"the figure panels use 29,419.")
    s = pd.read_csv(p, index_col=0)
    bcs = {parse_barcode(x) for x in s.index}
    bcs.discard(None)
    return bcs


def stan_mask(obs_names, root, verbose=True):
    """Boolean mask over obs_names for cells present in Stan's barcode list."""
    bcs = stan_barcodes(root)
    parsed = np.array([parse_barcode(x) for x in obs_names], dtype=object)
    keep = np.array([b is not None and b in bcs for b in parsed])
    n_unparsed = int(sum(b is None for b in parsed))
    if verbose:
        print(f"Stan filter: {keep.sum():,}/{len(keep):,} cells kept "
              f"({100*keep.mean():.1f}%)"
              + (f"   [{n_unparsed:,} obs_names had no parseable barcode]"
                 if n_unparsed else ""))
    if keep.sum() == 0:
        raise RuntimeError("Stan filter kept ZERO cells -- barcode formats do not match. "
                           "Check that obs_names still contain the raw 16bp barcode.")
    if keep.mean() < 0.4:
        print(f"  WARNING: retention {100*keep.mean():.1f}% is lower than the expected "
              f"~64.7%. Verify the barcode list is the current one.")
    return keep


def read_obs_names(h5_path):
    """obs_names from an .h5ad without loading the matrix."""
    import h5py
    with h5py.File(h5_path, "r") as h:
        og = h["obs"]
        ik = og.attrs.get("_index", "_index")
        ik = ik.decode() if isinstance(ik, bytes) else ik
        it = og[ik]
        raw = it["values"][:] if isinstance(it, h5py.Group) else it[:]
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in raw])

# =============================================================================
# WHAT THE FILTER ACTUALLY IS  -- resolved 2026-08-09
# =============================================================================
# Until now this module treated `sample_barcodes.csv` as an opaque list. It is not.
# Stan confirmed by email that he applied ONLY a per-sample library-size gate, with
# the windows in LIBRARY_SIZE_GATE below, to the RAW (pre-QC) object.
#
# VERIFIED: applying those windows to 00_DATA_FROZEN/00_raw/amanda_liang_raw_unfiltered.h5ad
# reproduces his barcode list EXACTLY -- 29,519 of 29,519 cells, zero mismatches, all
# eight samples. See 05_notebooks/qc_supplement/.
#
# Two consequences worth knowing:
#
# 1. THE GATE IS THE PRIMARY DEFINITION. The CSV is now a verification artefact, not
#    a dependency. `library_size_gate()` needs no external file.
#
# 2. `stan_mask()` MATCHES ON BARCODE SEQUENCE ALONE, IGNORING THE SAMPLE. 154 barcodes
#    occur in more than one of the eight samples, so stan_mask keeps 72 cells on the
#    strength of a barcode belonging to a DIFFERENT sample:
#         stan_mask()          29,419 cells
#         library_size_gate()  29,347 cells
#    The difference is 0.24% and does not move the results -- lNasM 72h is 521 cells
#    either way, eNasM 72h goes 1,695 -> 1,689.
#
#    Every frozen output in 01_DERIVED_LOCKED and every current figure was built with
#    stan_mask(). It is therefore LEFT UNCHANGED so those stay reproducible. Use
#    library_size_gate() for anything new, and see verify_gate_matches_list().
#
# The gate is ASYMMETRIC between arms of the same comparison -- Michigan 72h is
# Ctrl 5,000-30,000 against KD 4,000-20,000, giving controls 1.67x the median depth.
# Depth-standardised composition is in the QC supplement; lNasM, pluripotent, pPS-LC
# and AE-LC all hold or strengthen, eNasM retains 58% of its effect.
# =============================================================================

# (min, max) total_counts, inclusive. Keys are the CORRECTED sample labels.
LIBRARY_SIZE_GATE = {
    "Chicago_Ctrl_48h":  (15_000, 35_000),
    "Chicago_DOX_48h":   (10_000, 35_000),
    "Chicago_Ctrl_72h":  ( 4_500, 20_000),
    "Chicago_DOX_72h":   ( 5_000, 20_000),
    "Michigan_Ctrl_48h": ( 7_000, 25_000),
    "Michigan_DOX_48h":  ( 8_000, 25_000),
    "Michigan_Ctrl_72h": ( 5_000, 30_000),
    "Michigan_DOX_72h":  ( 4_000, 20_000),
}


def library_size_gate(sample_id, total_counts):
    """Stan's filter, computed rather than looked up. No external file needed.

    sample_id     corrected labels, e.g. "Michigan_Ctrl_72h"
    total_counts  UMIs per cell

    Returns a boolean array. On the QC-passed master object this keeps 29,347 cells.
    """
    sample_id = np.asarray(sample_id, dtype=object)
    total_counts = np.asarray(total_counts, dtype=float)
    unknown = sorted(set(sample_id) - set(LIBRARY_SIZE_GATE))
    if unknown:
        raise KeyError(f"no library-size gate defined for {unknown}. Keys must be the "
                       f"CORRECTED sample labels -- note Ctrl_72h and DOX_48h were "
                       f"exchanged at both sites by the OCM correction.")
    lo = np.array([LIBRARY_SIZE_GATE[s][0] for s in sample_id], float)
    hi = np.array([LIBRARY_SIZE_GATE[s][1] for s in sample_id], float)
    return (total_counts >= lo) & (total_counts <= hi)


def verify_gate_matches_list(sample_id, total_counts, obs_names, root):
    """Prove the gate and the supplied barcode list agree. Returns a dict."""
    gate = library_size_gate(sample_id, total_counts)
    listed = stan_mask(obs_names, root, verbose=False)
    return {"gate": int(gate.sum()), "barcode_list": int(listed.sum()),
            "both": int((gate & listed).sum()),
            "gate_only": int((gate & ~listed).sum()),
            "list_only": int((~gate & listed).sum())}
