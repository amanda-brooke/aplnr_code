"""
kd_pseudobulk.py — CANONICAL pseudobulk DE for the shAPLNR knockdown.

Single source of truth. Replaces and consolidates:
    kd_pseudobulk_deseq2.py                  (superseded by this file)
    FIGURE_KD_FIGURE.ipynb  cell 36          (duplicate pydeseq2, different params,
                                              wrote the SAME filename -- delete it)
    FIGURE_KD_FIGURE.ipynb  concordance cell (raw-CPM pseudobulk -- now reads this output)
    19_pseudobulk_de.ipynb                   (exploratory; produced 06_tables/pseudobulk/)

Run in the pseudobulk env:  python kd_pseudobulk.py
Needs: numpy, pandas, h5py, scipy, pydeseq2.  Reads the master via h5py, so it is
immune to the anndata uns/log1p read bug.

-------------------------------------------------------------------------------
DESIGN NOTES -- read before changing the parameters
-------------------------------------------------------------------------------
Chicago and Michigan are independent BIOLOGICAL REPEATS. The design is therefore
n = 2, paired by site, and `~site + cond` is the correct blocked model.

This script does NOT report a genome-wide padj. The reason was stated wrongly in an
earlier version of this docstring, which claimed genome-wide FDR was simply
unreachable. Measured: it is reachable, and that is the problem.

  tissue        n_pb   genes    best p     padj<0.05   best gene
  AE-LC            7  15,673   6.7e-06            0    CDKN1C
  pluripotent      8  19,367   9.5e-08            2    CDKN1C
  NasM             4  17,185   1.9e-09          223    LRRC75A
  aPS-LC           4  16,044   9.4e-16          351    H4C3
  pPS-LC           4  16,008   1.8e-34          378    TUBB
  all_mesoderm     4  20,124   3.3e-16          524    ANKRD1
  pluri_AE         4  19,690  2.7e-112          173    HNRNPDL

Two things are wrong with the bottom five rows. First, p = 2.7e-112 from FOUR samples
is not a finding, it is a broken variance estimate: those groups have 1 residual df,
so the Wald test leans entirely on the dispersion prior, and when that prior is too
tight the p-values collapse. The 7-8 pseudobulk groups, which have real residual df,
behave sanely (0 and 2 significant genes). The split is by DESIGN, not by biology.

Second, look at WHAT is significant. Of the top 50 genome-wide genes, 21/50 (NasM),
27/50 (pPS-LC) and 18/50 (pluri_AE) are ribosomal, mitochondrial, snRNP, histone or
tubulin -- TUBB, H4C3, COMMD6, SNRPG, FTH1, RPL39, RPS13. That is the compositional
shift from KD cells being 25-38% less complex, not APLNR biology.

Sign agreement does NOT rescue these: the top-50 agree 50/50 across repeats, but the
compositional shift is present at BOTH sites, so an artefact agrees in sign too.

So: no genome-wide padj. Use the pre-specified panel, and the permutation-based set
tests in kd_set_enrichment.py, both of which are immune to this failure mode.
APLNR itself ranks #3 of 17,185 in NasM, which is meaningful even though the absolute
p-value is not.

What n = 2 DOES support, and what this script reports instead:
  * shrunk effect sizes (LFC shrinkage stabilises log2FC when dispersion is poorly
    estimated, which is exactly the n = 2 situation)
  * per-repeat log2FC computed on DESeq2 median-of-ratios normalised counts, then
    MEDIAN-CENTRED.  Centring matters: raw CPM log2FCs on this data are shifted
    -0.285 genome-wide because KD cells are 25-38% less complex, which makes 76% of
    all genes negative in both repeats and sign agreement nearly automatic. Centring
    cuts "reproducible" genes from 4,164 to 490 and makes the survivors meaningful.
  * a reproducibility flag: same sign in both repeats AND |log2FC| above a floor.

    *** DO NOT USE THIS FLAG AS A DISCOVERY FILTER. ***
    Measured by permuting condition labels across cells within each site (12 permutations,
    NasM, 72h) the flag has a high and irreducible false-discovery rate:

        |log2FC| floor   real genes   null median   empirical FDR   key genes retained
             0.50            1,266          364          29%              12/12
             0.75              338           88          26%               3/12
             1.00              108           18          17%               1/12
             1.25               35            4          11%               0/12

    There is no threshold that is both clean and useful. The noise floor is set by
    Chicago_Ctrl_72h having only 53 NasM cells against Michigan's 858, and no downstream
    filter can repair that. A related symptom: of the 1,266 genes at floor 0.5, 89 of the
    863 "up" are ribosomal/OXPHOS/translation machinery against 6 of the 403 "down" -- the
    compositional shift surviving the median-centring.

    The flag IS fine for genes named in advance on biology, because the artefact does not
    preferentially target them. It is NOT fine for "here are the N genes that changed",
    a volcano, or GSEA input. For formal inference use KD_pseudobulk_panel_tests.csv,
    which BH-corrects within a pre-specified 23-gene panel -- a different and valid frame.

    ** DO NOT compare the reproducible COUNT between groups with different
    n_pseudobulks. ** Groups fitted over both timepoints (7-8 pseudobulks) average two
    samples per site x condition when the per-repeat log2FC is formed, which halves the
    noise and mechanically lowers the count. Measured here: 1.3-2.1% reproducible for the
    7-8 pseudobulk groups vs 6.7-9.9% for every 4-pseudobulk group, mesoderm and
    non-mesoderm alike. The count tracks the design, not the biology. Compare it only
    within groups that share n_pseudobulks, and use the panel tests for specificity.
  * targeted p-values on a small PRE-SPECIFIED mechanistic panel, BH-corrected within
    the panel. This is the only formal test here and it is labelled as such.

MIN_CELLS_PER_SAMPLE = 35 is chosen from the data, not by convention. The binding
constraint is Chicago_Ctrl_72h NasM, which has 53 cells BEFORE Stan's barcode filter
and 36 AFTER it. The filter is now applied (see stan_filter.py), so a floor of 50 would
drop that sample, leaving 3 pseudobulks and zero residual df -- unfittable. At >=35 all
four 72h NasM pseudobulks survive and the paired 2x2 design is preserved.
Do not raise this without re-checking KD_pseudobulk_QC.csv.

Runtime: ~25 s to load the master, then ~35 s per tissue (most of it LFC shrinkage).
About 4-5 minutes end to end for the seven groups in STATE_PLAN.

Expected warnings
-----------------
"overflow encountered in exp"  -- benign, pydeseq2's IRLS inner loop. Suppressed below.
"residual degrees of freedom is less than 3 ... poorly estimated by the MAD"
    -- REAL, and deliberately left visible. This is the n = 2 design speaking: with 4
    pseudobulks and 1 residual df in the 72h mesoderm groups, the dispersion prior is
    shaky. It is why this script emits no genome-wide padj, and it means the panel
    p-values are approximate too. The per-repeat concordance columns do not depend on
    the dispersion estimate and are the more robust evidence.
"""
import os

# pydeseq2 fits genes in joblib subprocesses, and warning filters set in this process do
# NOT propagate to them -- which is why filterwarnings() alone leaves the console full of
# "overflow encountered in exp". Setting PYTHONWARNINGS before joblib spawns its workers
# does propagate. RuntimeWarning only: this deliberately leaves pydeseq2's UserWarnings
# visible, in particular the "residual degrees of freedom is less than 3" one, which is a
# genuine signal about the n = 2 design and should not be hidden.
os.environ.setdefault("PYTHONWARNINGS", "ignore::RuntimeWarning")

import numpy as np, pandas as pd, h5py, scipy.sparse as sp, time, warnings
from pathlib import Path
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
from config import ROOT
from stan_filter import stan_mask, read_obs_names

warnings.filterwarnings("ignore", message="overflow encountered in exp")
np.seterr(over="ignore")

# ----------------------------------------------------------------------------- config
MASTER = ROOT / "00_DATA_FROZEN" / "03_analysis_ready" / "adata_master_normalized_magic.h5ad"
TAB    = ROOT / "06_tables" / "FIGURE_KD_clean"; TAB.mkdir(parents=True, exist_ok=True)

MIN_CELLS_PER_SAMPLE = 35      # see design notes; 50 breaks NasM AFTER the Stan filter
MIN_COUNT            = 10      # a gene needs >= this many counts ...
MIN_FRAC_SAMPLES     = 0.5     # ... in at least this fraction of pseudobulks
REPRO_LFC_FLOOR      = 0.5     # |centred log2FC| required in BOTH repeats

# Which states to test, and over which timepoints.
# Mesoderm is 72h-only: at 42h Chicago_DOX has ZERO nascent-mesoderm cells and
# Michigan_DOX has 10 (recount 2026-08-07, AFTER the Stan filter: eNasM 9 + lNasM 1;
# the previous figure of 30 was pre-filter), so a timepoint term there fits noise.
#
# aPS-LC and AE-LC stay 72h-only. Not a cell-count question: at 42h Chicago_DOX has
# 1,012 pluripotent cells against 1,024 control -- i.e. the cells are there -- but only
# 33 aPS-LC (vs 364) and 15 AE-LC (vs 133). That sample under-differentiated. Lowering
# MIN_CELLS_PER_SAMPLE to admit it would fit a model on a failed differentiation.
# If aPS-LC at 42h is ever needed, run it Michigan-only (1,166 vs 1,186) and state the
# Chicago exclusion.
STATE_PLAN = {
    "pluripotent":  {"states": ["pluripotent"],                              "timepoints": ["42h", "72h"]},
    "AE-LC":        {"states": ["AE-LC"],                                    "timepoints": ["42h", "72h"]},
    "aPS-LC":       {"states": ["aPS-LC"],                                   "timepoints": ["72h"]},
    # pPS-LC gets 42h too (added 2026-08-07): post-Stan-filter all four 42h groups
    # clear MIN_CELLS_PER_SAMPLE -- Chi 237/105, Mich 539/580. It is the only
    # differentiated state besides pluripotent where Chicago_DOX_42h is usable.
    "pPS-LC":       {"states": ["pPS-LC"],                                   "timepoints": ["42h", "72h"]},
    "NasM":         {"states": ["eNasM", "lNasM"],                           "timepoints": ["72h"]},
    "all_mesoderm": {"states": ["aPS-LC", "pPS-LC", "eNasM", "lNasM"],       "timepoints": ["72h"]},
    # pooled non-mesoderm reference at 72h, used as the contrast for
    # mesoderm_specific_down.csv in FIGURE_KD_CLEAN.ipynb. Kept 72h-only so it is
    # like-for-like with all_mesoderm.
    "pluri_AE":     {"states": ["pluripotent", "AE-LC"],                     "timepoints": ["72h"]},
}

# Pre-specified on biology, before seeing any p-value. The only formal test here.
#
# AMENDMENT 2026-07-25: SNAI1 and SERPINE2 added (23 -> 25 genes).
# These two carry the EMT limb of the stated hypothesis and were named in the plot
# plan's headline sentence, but had been omitted from the panel by oversight -- so
# the headline was asserting an effect for two genes that had no p-value at all.
# Adding them is a correction, not an expansion, and it costs the whole panel a
# slightly stricter BH threshold, which is the honest trade.
#
# *** DO NOT ADD GENES THAT CAME OUT OF THE APLNR CORRELATION ANALYSIS. ***
# BMP2, GATA6, VIM, HAND1, GJA1 and IRX3 all fall in the KD and all look tempting.
# They were identified BY looking at the data, so putting them in a "pre-specified"
# panel would be circular and would destroy the one clean inferential frame in this
# figure. They are reported instead as a SET-LEVEL enrichment (arm 3), which is a
# valid test precisely because it does not depend on picking individual winners.
PANEL = ["APLNR", "APELA", "WNT8A", "WNT3", "WNT3A", "FZD7", "LRP5",
         "AXIN2", "SP5", "LEF1", "NKD1", "DKK1",
         "TBXT", "EOMES", "MIXL1", "TBX6", "MSGN1", "DLL3", "MESP1", "MESP2",
         "CDX2", "FOXC1", "FOXC2", "SNAI1", "SERPINE2",
         # 2026-07-26: RSPO3 and MSGN1 added (25 -> 27). Both are posterior-PSM genes and
         # both belong in the score's own audit trail: RSPO3 is IN the paraxial marker
         # score, and MSGN1 is the canonical readout of RSPO3-driven posterior PSM fate
         # (Chal et al. 2015). They were missing, so neither had a p-value.
         # THIS IS THE FINAL PANEL. Adding genes shifts the BH denominator and therefore
         # every panel_padj in the file. Do not add more without re-freezing.
         # MSGN1 is NOT repeated here -- it is already in the fate block above. Adding it
         # twice silently duplicated its row in every downstream merge.
         "RSPO3"]
assert len(PANEL) == len(set(PANEL)), f"duplicate panel genes: {[g for g in set(PANEL) if PANEL.count(g) > 1]}"


# ------------------------------------------------------------------------------ io
def _read_col(it):
    if isinstance(it, h5py.Group):
        ks = set(it.keys())
        if {"categories", "codes"} <= ks:
            return pd.Categorical.from_codes(it["codes"][:],
                                             categories=np.asarray(_read_col(it["categories"])))
        if "values" in ks:
            v = it["values"][:]
            return (np.array([x.decode() if isinstance(x, bytes) else x for x in v], dtype=object)
                    if v.dtype == object or v.dtype.kind in "SU" else v)
    a = it[:]
    return a.astype(str) if a.dtype.kind == "S" else a


def read_master(path):
    """Raw integer counts + the obs columns we need. Never touches X (log-normalised)."""
    with h5py.File(path, "r") as h:
        X = h["layers"]["counts"]                       # RAW COUNTS -- required by DESeq2
        M = sp.csr_matrix((X["data"][:], X["indices"][:], X["indptr"][:]),
                          shape=tuple(int(s) for s in X.attrs["shape"]))
        v = h["var"]; ik = v.attrs.get("_index", "_index")
        ik = ik.decode() if isinstance(ik, bytes) else ik
        vit = v[ik]; vraw = vit["values"][:] if isinstance(vit, h5py.Group) else vit[:]
        var = np.array([(x.decode() if isinstance(x, bytes) else str(x)).split(" (")[0]
                        for x in vraw])
        og = h["obs"]
        obs = pd.DataFrame({c: np.asarray(_read_col(og[c])).astype(str) for c in
                            ["condition", "site", "timepoint_true", "predicted_leiden_sub_joint"]})
    if len(set(var)) != len(var):
        warnings.warn("duplicate gene symbols after parsing -- check the var index")

    # STAN BARCODE FILTER -- the single cell filter for this project. Applied here so
    # every downstream pseudobulk uses the SAME cell set as the figure panels. Without
    # it this script silently analysed 45,435 cells while FIGURE2_CLEAN used 29,419.
    keep = stan_mask(read_obs_names(path), ROOT)
    M = M[keep]
    obs = obs.loc[keep].reset_index(drop=True)
    return M, var, obs


# ---------------------------------------------------------------------- normalisation
def size_factors(counts):
    """DESeq2 median-of-ratios. Computed here rather than pulled off the dds object so
    the per-repeat log2FCs do not depend on pydeseq2's internal attribute names."""
    c = counts.to_numpy(dtype=float)                    # samples x genes
    with np.errstate(divide="ignore"):
        logc = np.log(c)
    keep = np.all(np.isfinite(logc), axis=0)            # genes detected in every sample
    if keep.sum() < 100:
        warnings.warn(f"only {keep.sum()} genes present in all samples; size factors may be unstable")
    ref = np.exp(logc[:, keep].mean(axis=0))            # geometric mean per gene
    return np.median(c[:, keep] / ref, axis=1)


# --------------------------------------------------------------- pseudobulk assembly
def build_pseudobulks(M, var, obs):
    """Sum raw counts per (sample x tissue). Returns {tissue: (cmat, mdf)} and a QC log."""
    obs = obs.copy()
    obs["cond"] = obs["condition"].map({"Ctrl": "Ctrl", "DOX": "KD"})
    if obs["cond"].isna().any():
        raise ValueError(f"unmapped condition values: "
                         f"{sorted(obs.loc[obs['cond'].isna(),'condition'].unique())}")
    obs["sample"] = obs["site"] + "_" + obs["cond"] + "_" + obs["timepoint_true"]

    built, qc_rows = {}, []
    for tissue, plan in STATE_PLAN.items():
        smask = (obs["predicted_leiden_sub_joint"].isin(plan["states"])
                 & obs["timepoint_true"].isin(plan["timepoints"])).to_numpy()
        rows, meta = [], []
        for s in sorted(pd.unique(obs.loc[smask, "sample"])):
            mask = smask & (obs["sample"] == s).to_numpy()
            n    = int(mask.sum())
            used = n >= MIN_CELLS_PER_SAMPLE
            qc_rows.append({"tissue": tissue, "sample": s, "n_cells": n, "included": used,
                            "reason": "" if used else
                                      f"< MIN_CELLS_PER_SAMPLE ({MIN_CELLS_PER_SAMPLE})"})
            if not used:
                continue
            rows.append(np.asarray(M[mask].sum(0)).ravel())
            g = obs[mask].iloc[0]
            meta.append({"sample": s, "site": g["site"],
                         "timepoint": g["timepoint_true"], "cond": g["cond"]})
        if not meta:
            built[tissue] = (None, pd.DataFrame()); continue
        mdf  = pd.DataFrame(meta).set_index("sample")
        cmat = pd.DataFrame(np.vstack(rows).astype(int), index=mdf.index, columns=var)
        built[tissue] = (cmat, mdf)
    return built, qc_rows


def fit_tissue(tissue, cmat, mdf):
    """DESeq2 + shrinkage + centred per-repeat log2FCs for one tissue.
    Returns (results_df, panel_df) or (None, None) if the design is not fittable."""
    if cmat is None or len(mdf) < 4 or mdf["cond"].nunique() < 2:
        print(f"SKIP {tissue}: {len(mdf)} usable pseudobulks", flush=True)
        return None, None
    paired = mdf.groupby("site")["cond"].nunique()
    if (paired < 2).any():
        print(f"SKIP {tissue}: sites without both conditions -> {paired.to_dict()}", flush=True)
        return None, None

    keep = (cmat >= MIN_COUNT).sum(0) >= max(2, int(np.ceil(MIN_FRAC_SAMPLES * len(cmat))))
    cmat = cmat.loc[:, keep]
    mdf  = mdf.copy()
    mdf["cond"] = pd.Categorical(mdf["cond"], categories=["Ctrl", "KD"])
    terms  = [f for f in ["site", "timepoint"] if mdf[f].nunique() > 1] + ["cond"]
    design = "~" + "+".join(terms)
    resid  = len(mdf) - (1 + sum(mdf[f].nunique() - 1 for f in terms))

    dds = DeseqDataSet(counts=cmat, metadata=mdf, design=design, quiet=True)
    dds.deseq2()
    ds = DeseqStats(dds, contrast=["cond", "KD", "Ctrl"], quiet=True)
    ds.summary()
    r = ds.results_df.copy().rename(columns={"log2FoldChange": "log2FC_raw"})

    # --- LFC shrinkage: the right move at n = 2, stabilises the effect estimates ---
    shrunk = None
    for coeff in ("cond[T.KD]", "cond_KD_vs_Ctrl"):
        try:
            ds.lfc_shrink(coeff=coeff)
            shrunk = ds.results_df["log2FoldChange"].copy(); break
        except Exception:
            continue
    if shrunk is None:
        warnings.warn(f"{tissue}: LFC shrinkage unavailable, using unshrunk log2FC")
        shrunk = r["log2FC_raw"]
    r["log2FC_shrunk"] = shrunk

    # --- per-repeat log2FC on normalised counts, MEDIAN-CENTRED (see design notes) ---
    sf   = size_factors(cmat)
    norm = cmat.div(sf, axis=0)
    for site in sorted(mdf["site"].unique()):
        c = norm[(mdf["site"] == site) & (mdf["cond"] == "Ctrl")].mean(0)
        k = norm[(mdf["site"] == site) & (mdf["cond"] == "KD")].mean(0)
        lfc = np.log2((k + 1) / (c + 1))
        r[f"log2FC_{site}"] = (lfc - lfc.median()).reindex(r.index)      # <-- centring

    sites = sorted(mdf["site"].unique())
    if len(sites) == 2:
        a, b = (r[f"log2FC_{s}"] for s in sites)
        r["sign_agree"]   = np.sign(a) == np.sign(b)
        r["reproducible"] = r["sign_agree"] & (np.minimum(a.abs(), b.abs()) > REPRO_LFC_FLOOR)
    r["tissue"], r["n_pseudobulks"], r["design"] = tissue, len(mdf), design
    r = r.drop(columns=[c for c in ["padj"] if c in r])       # no genome-wide FDR at n=2
    r = r.sort_values("pvalue")

    nrep = int(r["reproducible"].sum()) if "reproducible" in r else -1
    print(f"{tissue:13s} {design:24s} {len(mdf)} pseudobulks, {resid} resid df | "
          f"{len(r):5d} genes | {nrep:4d} reproducible", flush=True)

    # --- targeted panel: BH WITHIN the pre-specified panel only ---
    panel = r.reindex([g for g in PANEL if g in r.index])[["log2FC_shrunk", "pvalue"]].dropna()
    if len(panel):
        p = panel["pvalue"].to_numpy(float); n = len(p); o = np.argsort(p)
        q = np.empty(n); q[o] = p[o] * n / (np.arange(n) + 1)
        q[o] = np.minimum.accumulate(q[o][::-1])[::-1]
        panel = panel.assign(panel_padj=np.clip(q, 0, 1), tissue=tissue)
        sig = panel.index[panel["panel_padj"] < 0.05].tolist()
        print(f"              panel-FDR<0.05: {sig if sig else 'none'}", flush=True)
    return r, (panel if len(panel) else None)

# ------------------------------------------------------- identity pseudobulk export
IDENTITY_STATES = [
    "pluripotent",
    "AE-LC",
    "pPS-LC",
    "aPS-LC",
    "eNasM",
    "lNasM",
    "aPM-LC",
    "LPM-LC",
]

IDENTITY_TIMEPOINTS = ["42h", "72h"]


def export_identity_pseudobulks(M, var, obs, outdir):
    """
    Export exact-state raw-count pseudobulks for comparison with directed
    differentiation. This is separate from the KD differential-expression design.
    """
    obs = obs.copy()

    obs["condition_export"] = obs["condition"].map({
        "Ctrl": "CTRL",
        "DOX": "KD",
    })

    if obs["condition_export"].isna().any():
        bad = sorted(
            obs.loc[obs["condition_export"].isna(), "condition"].unique()
        )
        raise ValueError(f"Unmapped conditions: {bad}")

    state_col = "predicted_leiden_sub_joint"
    time_col = "timepoint_true"

    keep = (
        obs[state_col].isin(IDENTITY_STATES)
        & obs[time_col].isin(IDENTITY_TIMEPOINTS)
    )

    rows = []
    metadata = []

    grouping_columns = [
        "site",
        "condition_export",
        time_col,
        state_col,
    ]

    for keys, indices in obs.loc[keep].groupby(
        grouping_columns,
        observed=True
    ).groups.items():

        site, condition, timepoint, state = keys
        indices = np.asarray(list(indices), dtype=int)
        n_cells = len(indices)

        pseudobulk_id = (
            f"{site}_{condition}_{state}_{timepoint}"
        )

        summed_counts = np.asarray(
            M[indices].sum(axis=0)
        ).ravel()

        rows.append(
            pd.Series(summed_counts, index=var, name=pseudobulk_id)
        )

        metadata.append({
            "pseudobulk_id": pseudobulk_id,
            "site": site,
            "condition": condition,
            "timepoint": timepoint,
            "state": state,
            "n_cells": n_cells,
            "low_cell_count": n_cells < 20,
        })

    # Samples × genes initially
    count_matrix = pd.DataFrame(rows)

    # Export as genes × pseudobulk samples
    count_matrix.T.rename_axis("gene").to_csv(
        outdir / "gastruloid_identity_pseudobulk_counts.csv"
    )

    pd.DataFrame(metadata).to_csv(
        outdir / "gastruloid_identity_pseudobulk_metadata.csv",
        index=False,
    )

    print(
        "\nIdentity pseudobulk export:"
        f"\n  {outdir / 'gastruloid_identity_pseudobulk_counts.csv'}"
        f"\n  {outdir / 'gastruloid_identity_pseudobulk_metadata.csv'}"
    )


# ---------------------------------------------------------------------------- main
def main():
    t0 = time.time()
    M, var, obs = read_master(MASTER)
    print(f"loaded {M.shape} in {round(time.time()-t0)}s", flush=True)
    export_identity_pseudobulks(M, var, obs, TAB)
    
    built, qc_rows = build_pseudobulks(M, var, obs)

    results, panel_rows = {}, []
    for tissue, (cmat, mdf) in built.items():
        r, panel = fit_tissue(tissue, cmat, mdf)
        if r is not None:
            results[tissue] = r
            if panel is not None:
                panel_rows.append(panel)

    if not results:
        print("no tissue produced a fittable model"); return

    out = pd.concat(results.values())
    out.index.name = "gene"
    out.to_csv(TAB / "KD_pseudobulk_results.csv")
    pd.concat(panel_rows).rename_axis("gene").to_csv(TAB / "KD_pseudobulk_panel_tests.csv")
    pd.DataFrame(qc_rows).to_csv(TAB / "KD_pseudobulk_QC.csv", index=False)

    print(f"\nwrote:\n  {TAB/'KD_pseudobulk_results.csv'}"
          f"\n  {TAB/'KD_pseudobulk_panel_tests.csv'}"
          f"\n  {TAB/'KD_pseudobulk_QC.csv'}")
    print("\nNOTE: there is deliberately no genome-wide padj column. At n = 2 the best gene "
          "\nper tissue falls 10-100x short of what BH over ~20,000 genes requires, so a padj "
          "\ncolumn here would only ever read 1.0 and be misinterpreted. Use log2FC_shrunk, "
          "\nthe per-repeat centred log2FCs and `reproducible`; the panel file holds the only "
          "\nformal tests.")


if __name__ == "__main__":
    main()
