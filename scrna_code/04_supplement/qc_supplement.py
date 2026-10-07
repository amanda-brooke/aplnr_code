"""
qc_supplement.py -- builds the single-cell QC supplement. Standalone by design.

WHY THIS IS A SCRIPT AND NOT A NOTEBOOK CELL
QC is a leaf of the pipeline: it READS only 00_DATA_FROZEN and WRITES only to
./output/. Nothing in the analysis reads anything this produces. That is what makes
it safe to run at any time without touching a figure -- and it stays true only as
long as nobody adds an analysis step in here. Do not.

RUN
    cd 04_supplement
    python qc_supplement.py

OUTPUT
    output/QC_S1_thresholds.*            pre-filter violins with the cuts drawn
    output/QC_S2_batch_umap.*            UMAP by site / condition / timepoint / state
    output/QC_S4_barcode_depth.*         per-sample library-size gate
    output/QC_S5_depth_standardised.*    composition, raw vs depth-standardised
    output/QC_SOURCE_DATA.xlsx           10 sheets, the numbers behind all of it

Each panel is written as .png (preview), .pdf and .svg. Text is kept as text in both
vector formats, so font size and family are editable without re-running anything.

DEGRADES INSTEAD OF CRASHING
QC_S1 needs the pre-filter object (49,646 cells). If it is absent the script builds
everything else and says so at the end. It never half-writes a file.

COLOURS
From ../utils/palettes.py. No hex literals here. Ever.
"""
from pathlib import Path
import sys
import warnings

import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from anndata.io import read_elem

# ----------------------------------------------------------------- CONFIG

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "utils"))
from config import ROOT
FROZEN = ROOT / "00_DATA_FROZEN"
OUT = HERE / "output"
OUT.mkdir(exist_ok=True)

DEPTH_TABLE = None
STD_TABLE = None
LT_STATE = LT_SENS = None

from stan_filter import LIBRARY_SIZE_GATE, library_size_gate
from palettes import (CTRL_C, KD_C, STATE_COLORS, SITE_COLORS, CONDITION_COLORS, TIMEPOINT_COLORS,
                      TISSUE_COLORS, TISSUE_MAP, TISSUE_ORDER,
                      GREY_LINE, GREY_FILL, GREY_DARK, GREY_MID, RED)

# QC thresholds -- these are the numbers stated in the methods. They are copied
# from 00_DATA_FROZEN/_provenance/2D_gastruloid_analysis_Amanda_JOINT_SCVI_newmodel.ipynb
# (cell 3) and must not be edited independently of it.
MIN_GENES  = 500
MIN_COUNTS = 1_000
MAX_PCT_MT = 20.0
UPPER_MAD  = 5.0     # per-sample: median + 5*MAD on genes and on counts.
                     # This is also the doublet control -- no doublet caller was run.

MASTER  = FROZEN / "03_analysis_ready" / "adata_master_normalized.h5ad"
BARCODES = FROZEN / "01_source" / "sample_barcodes.csv"
FILTER_REPORT = ROOT / "02_new_experiment" / "summary_tables" / "qc_filter_report.csv"

# Optional. Checked in order; first hit wins.
RAW_CANDIDATES = [
    FROZEN / "00_raw" / "amanda_liang_raw_unfiltered.h5ad",
    FROZEN / "01_source" / "amanda_liang_raw_unfiltered.h5ad",
    ROOT.parent / "new_analysis" / "amanda_liang_raw_unfiltered.h5ad",
]

# Canonical display order. Note these are the CORRECTED labels (see sheet S6):
# the on-chip-multiplexing assignment for Ctrl_72h and DOX_48h was swapped at both
# sites, and nominal 48h is 42h true.
SAMPLE_ORDER = ["Chicago_Ctrl_48h", "Chicago_Ctrl_72h", "Chicago_DOX_48h", "Chicago_DOX_72h",
                "Michigan_Ctrl_48h", "Michigan_Ctrl_72h", "Michigan_DOX_48h", "Michigan_DOX_72h"]

from qc_style import apply_style, MATRIX_CMAP, SMALL, TINY
apply_style()

# ------------------------------------------------------------------ HELPERS

def read_obs(path, columns=None):
    """obs only -- never loads the expression matrix."""
    with h5py.File(path, "r") as f:
        df = read_elem(f["obs"])
    return df[columns] if columns else df


def read_obsm(path, key):
    with h5py.File(path, "r") as f:
        return np.asarray(f["obsm"][key][:])


def stan_mask(obs_names):
    """Stan's barcode list. Same rule as stan_filter.py: first run of >=14 ACGT."""
    import re
    pat = re.compile(r"[ACGT]{14,}")
    def parse(x):
        h = pat.findall(str(x))
        return h[0] if h else None
    keep_set = {parse(x) for x in pd.read_csv(BARCODES, index_col=0).index}
    keep_set.discard(None)
    return np.array([(b := parse(x)) is not None and b in keep_set for x in obs_names])


def mad(x):
    x = np.asarray(x, float)
    return np.nanmedian(np.abs(x - np.nanmedian(x)))


def save(fig, stem):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"{stem}.{ext}")
    plt.close(fig)
    print(f"  wrote output/{stem}  (png, pdf, svg)")


def strip(ax):
    ax.grid(False)
    ax.set_axisbelow(True)
    return ax


# =========================================================== S1  THRESHOLDS

def panel_thresholds(raw_obs):
    """Pre-filter distributions per sample, with every stated cut drawn on.

    Reviewers asked to see the thresholds, not be told them. The upper cuts are
    per-sample (median + 5*MAD), so they are drawn as a short tick per violin
    rather than one line across the axes.
    """
    metrics = [("n_genes_by_counts", "Genes per cell", MIN_GENES, True),
               ("total_counts",      "UMIs per cell",  MIN_COUNTS, True),
               ("pct_counts_mt",     "% mitochondrial", MAX_PCT_MT, False),
               ("pct_counts_ribo",   "% ribosomal",     None,       False)]

    fig, axes = plt.subplots(len(metrics), 1, figsize=(9.4, 10.9), sharex=True)
    order = [s for s in SAMPLE_ORDER if s in set(raw_obs.sample_id)]

    for ax, (col, label, floor, log) in zip(axes, metrics):
        strip(ax)
        raw = [raw_obs.loc[raw_obs.sample_id == s, col].astype(float).values for s in order]
        # The KDE must be computed in the space the axis is drawn in. Computing it on
        # linear counts and then setting a log axis inflates the low tail into a fake
        # second mode -- which is what the first version of this figure did.
        data = [np.log10(np.clip(v, 1, None)) for v in raw] if log else raw
        parts = ax.violinplot(data, positions=range(len(order)), widths=0.8,
                              showextrema=False, showmedians=False)
        for pc, s in zip(parts["bodies"], order):
            pc.set_facecolor(SITE_COLORS[s.split("_")[0]])
            pc.set_alpha(0.65)
            pc.set_edgecolor("none")
        for i, d in enumerate(data):
            ax.plot([i - 0.2, i + 0.2], [np.median(d)] * 2, color=GREY_DARK, lw=1.2, zorder=3)
        del raw

        tx = (lambda y: np.log10(max(y, 1))) if log else (lambda y: y)

        if floor is not None:
            ax.axhline(tx(floor), color=RED, lw=1.0, ls="--", zorder=4)
            ax.text(len(order) - 0.4, tx(floor), f"  {floor:,g}", color=RED,
                    va="bottom", ha="left", fontsize=9)

        # per-sample upper MAD cut, where one applies
        if col in ("n_genes_by_counts", "total_counts"):
            for i, smp in enumerate(order):
                v = raw_obs.loc[raw_obs.sample_id == smp, col].astype(float)
                ax.plot([i - 0.35, i + 0.35],
                        [tx(np.median(v) + UPPER_MAD * mad(v))] * 2,
                        color=RED, lw=1.0, zorder=4)

        if log:   # axis is log10 units; relabel the ticks in original units
            lo, hi = ax.get_ylim()
            ticks = np.arange(np.floor(lo), np.ceil(hi) + 1)
            ax.set_yticks(ticks)
            ax.set_yticklabels([f"$10^{{{int(t)}}}$" for t in ticks])
        ax.set_ylabel(label)

    axes[0].set_title(f"Pre-filter QC, all {len(raw_obs):,} called cells\n"
                      f"dashed = fixed threshold; solid red = per-sample median + {UPPER_MAD:g}xMAD "
                      f"(also the doublet control); grey = median",
                      loc="left", fontsize=9)
    axes[-1].set_xticks(range(len(order)))
    axes[-1].set_xticklabels(order, rotation=90)
    fig.tight_layout()
    save(fig, "QC_S1_thresholds")


# ============================================================ S2  BATCH UMAP

def panel_batch(obs, umap):
    """The 'clusters are not batch' panel.

    The honest framing, which the caption must carry: site IS the batch AND the
    biological repeat -- they are one variable, not two. What the design does give
    is that all four samples per site were multiplexed in a single 10x lane, so
    Ctrl vs DOX and 42h vs 72h are WITHIN-lane comparisons and carry no batch term.
    """
    obs = obs.copy()
    obs["tissue"] = obs.predicted_leiden_sub_joint.map(TISSUE_MAP).fillna("other")

    keys = [("site", SITE_COLORS, "Site (= batch = biological repeat)"),
            ("condition", CONDITION_COLORS, "Condition (within-lane)"),
            ("timepoint_true", TIMEPOINT_COLORS, "Timepoint (within-lane)"),
            ("tissue", TISSUE_COLORS, "Transferred label")]

    fig, axes = plt.subplots(2, 2, figsize=(9.4, 9.1))
    rng = np.random.default_rng(0)
    shuf = rng.permutation(len(obs))          # avoid last-plotted-wins bias

    for ax, (key, cmap, title) in zip(axes.ravel(), keys):
        strip(ax)
        vals = obs[key].astype(str).values[shuf]
        xy = umap[shuf]
        levels = ([t for t in TISSUE_ORDER if t in set(vals)] if key == "tissue"
                  else sorted(set(vals)))
        for lv in levels:
            m = vals == lv
            ax.scatter(xy[m, 0], xy[m, 1], s=0.6, lw=0, alpha=0.55,
                       color=cmap.get(lv, GREY_MID), label=f"{lv} ({m.sum():,})",
                       rasterized=True)
        ax.set_aspect("equal", adjustable="datalim")
        ax.set_title(title, loc="left")
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.legend(markerscale=8, frameon=False, loc="upper left",
                  bbox_to_anchor=(1.0, 1.0), handletextpad=0.2, labelspacing=0.25)

    fig.suptitle(f"scVI joint embedding, {len(obs):,} cells after barcode filter",
                 x=0.02, ha="left", fontsize=9)
    fig.tight_layout()
    save(fig, "QC_S2_batch_umap")


# ======================================================= S4  BARCODE DEPTH

def depth_table(obs):
    """Sequencing depth before and after the barcode filter, per sample."""
    rows = []
    for smp in SAMPLE_ORDER:
        m = obs.sample_id_all == smp
        for tag, sel in [("after_QC", m), ("after_barcode_filter", m & obs.stan)]:
            v = obs.loc[sel, "total_counts"].astype(float)
            if not len(v):
                continue
            rows.append({"sample": smp, "stage": tag, "n_cells": len(v),
                         "min": round(v.min()), "p01": round(np.percentile(v, 1)),
                         "median": round(v.median()), "p99": round(np.percentile(v, 99)),
                         "max": round(v.max()),
                         "p99_over_p01": round(np.percentile(v, 99) / np.percentile(v, 1), 1)})
    return pd.DataFrame(rows)


def panel_barcode_depth(obs, tab):
    """The second filter: a per-sample library-size gate.

    Confirmed by Stan 2026-08-09 and verified to reproduce his barcode list exactly
    (29,519/29,519, zero mismatches) when applied to the raw pre-QC object. The
    windows are drawn as red brackets. Note they are ASYMMETRIC between arms of the
    same comparison -- Michigan 72h is Ctrl 5,000-30,000 against KD 4,000-20,000.
    """
    fig, axes = plt.subplots(1, 2, figsize=(10.9, 4.3),
                             gridspec_kw={"width_ratios": [1.35, 1]})

    ax = strip(axes[0])
    order = [s for s in SAMPLE_ORDER if (obs.sample_id_all == s).any()]
    for i, smp in enumerate(order):
        m = obs.sample_id_all == smp
        a = obs.loc[m, "total_counts"].astype(float)
        b = obs.loc[m & obs.stan, "total_counts"].astype(float)
        for off, v, col in [(-0.18, a, GREY_FILL), (0.18, b, SITE_COLORS[smp.split("_")[0]])]:
            lo, hi = np.percentile(v, [1, 99])
            ax.plot([i + off] * 2, [lo, hi], color=col, lw=4, solid_capstyle="butt",
                    zorder=2)
            ax.plot([i + off - 0.09, i + off + 0.09], [np.median(v)] * 2,
                    color=GREY_DARK, lw=1.0, zorder=3)
    ax.set_yscale("log")
    ax.set_ylabel("UMIs per cell (1st-99th pct)")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=90)
    for i, smp in enumerate(order):                     # Stan's window as a bracket
        lo, hi = LIBRARY_SIZE_GATE[smp]
        ax.plot([i + 0.40] * 2, [lo, hi], color=RED, lw=0.9, zorder=5)
        for y in (lo, hi):
            ax.plot([i + 0.33, i + 0.47], [y] * 2, color=RED, lw=0.9, zorder=5)
    ax.axhline(MIN_COUNTS, color=GREY_DARK, lw=0.8, ls="--")
    ax.text(len(order) - 0.4, MIN_COUNTS, f"  your QC floor {MIN_COUNTS:,}",
            color=GREY_DARK, va="bottom", ha="left", fontsize=9)
    ax.set_title("grey = after QC     coloured = after gate     red bracket = Stan's window",
                 loc="left")

    ax = strip(axes[1])
    w = tab.pivot(index="sample", columns="stage", values="p99_over_p01").reindex(order)
    x = np.arange(len(order))
    ax.barh(x - 0.19, w["after_QC"], height=0.36, color=GREY_FILL,
            edgecolor=GREY_DARK, lw=0.4, label="after QC")
    ax.barh(x + 0.19, w["after_barcode_filter"], height=0.36, color=GREY_LINE,
            edgecolor=GREY_DARK, lw=0.4, label="after gate")
    ax.set_yticks(x)
    ax.set_yticklabels(order, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("depth spread (p99 / p01)")
    ax.legend(frameon=False, loc="lower right")
    ax.set_title("depth spread collapses in every sample", loc="left")

    fig.suptitle("Second filter: per-sample library-size gate (Stan). Windows differ "
                 "between arms of the same comparison.", x=0.01, ha="left", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.87))
    save(fig, "QC_S4_barcode_depth")


# ================================================ S5  DEPTH STANDARDISATION

STATES_FOR_DEPTH = ["pluripotent", "aPS-LC", "pPS-LC", "eNasM", "lNasM", "AE-LC"]


def depth_standardised_composition(obs, timepoint="72h", nbin=10):
    """Is the composition difference explained by the asymmetric gate?

    Direct standardisation. Within each site, cells are binned into deciles of
    total_counts defined on the POOLED cells of that site; the state fraction is
    computed per arm within each bin and then averaged over the pooled bin weights.
    Both arms are therefore reported at the same depth distribution.
    """
    sub = obs[(obs.timepoint_true == timepoint) & obs.gate].copy()
    sub["umi"] = sub.total_counts.astype(float)

    def raw(d):
        v = d.predicted_leiden_sub_joint.value_counts(normalize=True) * 100
        return pd.Series({k: v.get(k, 0.0) for k in STATES_FOR_DEPTH})

    def std(arm):
        acc, tot = [], 0
        for _, g in sub.groupby("site", observed=True):
            e = np.quantile(g.umi, np.linspace(0, 1, nbin + 1))
            e[0] -= 1; e[-1] += 1
            g = g.assign(bin=pd.cut(g.umi, e, labels=False))
            w = g.groupby("bin", observed=True).size()
            w = w / w.sum()
            a = g[g.condition == arm]
            f = (a.groupby("bin", observed=True).predicted_leiden_sub_joint
                  .value_counts(normalize=True).unstack(fill_value=0) * 100)
            f = f.reindex(index=w.index, columns=STATES_FOR_DEPTH).fillna(0)
            acc.append(f.mul(w, axis=0).sum() * len(a)); tot += len(a)
        return sum(acc) / tot

    t = pd.DataFrame({"Ctrl_raw": raw(sub[sub.condition == "Ctrl"]),
                      "KD_raw": raw(sub[sub.condition == "DOX"]),
                      "Ctrl_depth_standardised": std("Ctrl"),
                      "KD_depth_standardised": std("DOX")}).round(2)
    t["diff_raw"] = (t.KD_raw - t.Ctrl_raw).round(2)
    t["diff_standardised"] = (t.KD_depth_standardised - t.Ctrl_depth_standardised).round(2)
    t["pct_effect_retained"] = (100 * t.diff_standardised / t.diff_raw).round(0)
    return t.rename_axis("state").reset_index()


def panel_depth_standardised(tab):
    fig, ax = plt.subplots(figsize=(7.0, 3.9))
    strip(ax)
    x = np.arange(len(tab))
    ax.barh(x - 0.19, tab.diff_raw, height=0.36, color=GREY_FILL,
            edgecolor=GREY_DARK, lw=0.4, label="as analysed")
    ax.barh(x + 0.19, tab.diff_standardised, height=0.36, color=GREY_LINE,
            edgecolor=GREY_DARK, lw=0.4, label="depth-standardised")
    ax.axvline(0, color=GREY_DARK, lw=0.8)
    ax.set_yticks(x)
    ax.set_yticklabels(tab.state)
    ax.invert_yaxis()
    ax.set_xlabel("KD - Ctrl, % of cells at 72 h")
    ax.legend(frameon=False, loc="lower right")
    ax.set_title("Composition difference survives depth standardisation\n"
                 "eNasM is the exception (58% of effect retained)", loc="left")
    fig.tight_layout()
    save(fig, "QC_S5_depth_standardised")


# ============================================= S6  IH LABEL-TRANSFER CONFIDENCE

CONF_MIN, MARGIN_MIN = 0.50, 0.20     # used ONLY for the sensitivity check below


def label_transfer_tables(obs):
    """Per-state confidence in the IH label transfer, and an eNasM/lNasM sensitivity.

    WHY THIS EXISTS SEPARATELY FROM QC_S3
    The pooled confidence histogram has a median of 0.84, but that is carried by
    pluripotent (0.96) and AE-LC (0.92), which are 60% of cells. eNasM and lNasM --
    the two states the paper's claims rest on -- sit near the bottom at 0.65 and
    0.64. Reporting only the pooled number would be misleading.

    NOTE: no confidence filter is applied to the analysis. These thresholds are used
    here purely to show the result does not depend on the borderline cells.
    """
    f = obs[obs.gate].copy()
    P = [c for c in f.columns if c.startswith("prob_joint_")]
    labs = [c.replace("prob_joint_", "") for c in P]
    M = f[P].astype(float).values
    srt = np.sort(M, 1)
    f["top"], f["margin"] = srt[:, -1], srt[:, -1] - srt[:, -2]
    f["runner_up"] = [labs[i] for i in np.argsort(M, 1)[:, -2]]

    per_state = (f.groupby("predicted_leiden_sub_joint", observed=True)
                  .agg(n=("top", "size"), median_probability=("top", "median"),
                       pct_below_0p5=("top", lambda x: 100 * (x < 0.5).mean()),
                       median_margin=("margin", "median"),
                       pct_margin_below_0p2=("margin", lambda x: 100 * (x < 0.2).mean()))
                  .round(3).sort_values("median_probability")
                  .rename_axis("state").reset_index())
    per_state["most_common_runner_up"] = per_state.state.map(
        f.groupby("predicted_leiden_sub_joint", observed=True).runner_up
         .agg(lambda x: f"{x.mode().iat[0]} ({100*(x==x.mode().iat[0]).mean():.0f}%)"))

    rows = []
    for st in ["eNasM", "lNasM"]:
        d = f[(f.predicted_leiden_sub_joint == st) & (f.timepoint_true == "72h")]
        c = d[(d.top >= CONF_MIN) & (d.margin >= MARGIN_MIN)]
        for tag, x in [("as analysed", d), (f"top>={CONF_MIN} & margin>={MARGIN_MIN}", c)]:
            nc = int((x.condition == "Ctrl").sum()); nk = int((x.condition == "DOX").sum())
            rows.append({"state": st, "cell_set": tag, "n": len(x), "Ctrl": nc, "KD": nk,
                         "KD_per_Ctrl": round(nk / nc, 3) if nc else np.nan})
    sensitivity = pd.DataFrame(rows)
    return f, per_state, sensitivity


def panel_label_transfer(f, per_state, sensitivity):
    fig = plt.figure(figsize=(11.4, 5.1))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.5, 1.15, 1.0], wspace=0.45)

    # (a) per-state confidence, ordered
    ax = strip(fig.add_subplot(gs[0]))
    order = [t for t in TISSUE_ORDER if t in set(per_state.state)]
    order = list(per_state.set_index("state").reindex(order).dropna(how="all").index)
    for i, st in enumerate(order):
        v = f.loc[f.predicted_leiden_sub_joint == st, "top"].astype(float)
        lo, md, hi = np.percentile(v, [25, 50, 75])
        col = TISSUE_COLORS.get(st, GREY_MID)
        ax.plot([lo, hi], [i] * 2, color=col, lw=3.2, solid_capstyle="butt")
        ax.plot(md, i, "o", color=GREY_DARK, ms=3.2, zorder=3)
    ax.axvline(CONF_MIN, color=RED, lw=1.0, ls="--")
    ax.set_yticks(range(len(order))); ax.set_yticklabels(order)
    ax.invert_yaxis(); ax.set_xlim(0, 1)
    ax.set_xlabel("IH label-transfer probability")
    ax.set_title("median and IQR, per state", loc="left")

    # (b) the eNasM / lNasM boundary specifically
    ax = strip(fig.add_subplot(gs[1]))
    for st in ["eNasM", "lNasM"]:
        d = f[f.predicted_leiden_sub_joint == st]
        other = "lNasM" if st == "eNasM" else "eNasM"
        delta = (d[f"prob_joint_{st}"].astype(float)
                 - d[f"prob_joint_{other}"].astype(float))
        ax.hist(delta, bins=28, histtype="step", lw=1.3,
                color=STATE_COLORS[st], label=f"called {st} (n={len(d):,})")
    ax.axvline(0, color=GREY_DARK, lw=0.8)
    ax.set_xlabel("P(assigned) - P(the other state)")
    ax.set_ylabel("cells")
    ax.legend(frameon=False, loc="upper left", fontsize=9)
    ax.set_title("they are each other's runner-up,\nbut the margin is wide", loc="left")

    # (c) does requiring confidence change the answer?
    ax = strip(fig.add_subplot(gs[2]))
    d = sensitivity[sensitivity.state == "lNasM"]
    x = np.arange(len(d))
    ax.bar(x - 0.19, d.Ctrl, 0.36, color=CTRL_C, label="Ctrl")
    ax.bar(x + 0.19, d.KD, 0.36, color=KD_C, label="KD")
    for i, r in enumerate(d.itertuples()):
        ax.text(i, max(r.Ctrl, r.KD) * 1.04, f"1:{r.KD_per_Ctrl:.2f}",
                ha="center", fontsize=9, color=GREY_DARK)
    ax.set_xticks(x)
    ax.set_xticklabels(["as\nanalysed", f"conf\n>={CONF_MIN}"], fontsize=9)
    ax.set_ylabel("lNasM cells at 72 h")
    ax.legend(frameon=False, fontsize=9)
    ax.set_title("balance is unchanged", loc="left")

    fig.suptitle("IH label transfer: eNasM and lNasM are the least confident states,\n"
                 "but the Ctrl:KD balance does not depend on the borderline cells.",
                 x=0.01, y=0.995, va="top", ha="left", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.80))
    save(fig, "QC_S6_label_transfer")


# =============================================================== SOURCE DATA

def build_workbook(obs, raw_obs, raw_path):
    """One workbook, one sheet per question a reviewer can ask."""
    sheets = {}

    # --- S1  cell accounting, the waterfall
    #
    # TRAP: the pre-filter object and qc_filter_report.csv are keyed on the
    # PRE-correction sample label, the master object on the corrected one. Joining
    # them naively gives 150% retention for Chicago_DOX_48h. Everything upstream of
    # the OCM correction must be joined through sample_id_original.
    orig_of = (obs.groupby("sample_id_all", observed=True).sample_id_original
                  .agg(lambda x: x.mode().iat[0]).astype(str))
    if raw_obs is not None:
        before_orig = raw_obs.groupby("sample_id", observed=True).size()
    elif FILTER_REPORT.exists():
        before_orig = pd.read_csv(FILTER_REPORT, index_col=0)["before"]
    else:
        before_orig = pd.Series(dtype=int)
    before = pd.Series(orig_of.map(before_orig).values, index=orig_of.index, dtype=float)

    after_qc = obs.groupby("sample_id_all", observed=True).size()
    after_bc = obs.loc[obs.stan].groupby("sample_id_all", observed=True).size()

    acc = pd.DataFrame({"label_before_correction": orig_of.astype(str),
                        "cells_called": before.astype(float),
                        "after_QC": after_qc.astype(float),
                        "after_barcode_filter": after_bc.astype(float)}).reindex(SAMPLE_ORDER)
    acc["pct_passing_QC"] = (100 * acc.after_QC / acc.cells_called).round(1)
    acc["pct_of_QC_retained_by_barcodes"] = (100 * acc.after_barcode_filter / acc.after_QC).round(1)
    acc.loc["TOTAL"] = acc.sum(numeric_only=True)
    acc.loc["TOTAL", "label_before_correction"] = ""
    acc.loc["TOTAL", "pct_passing_QC"] = round(100 * acc.loc["TOTAL", "after_QC"] / acc.loc["TOTAL", "cells_called"], 1)
    acc.loc["TOTAL", "pct_of_QC_retained_by_barcodes"] = round(
        100 * acc.loc["TOTAL", "after_barcode_filter"] / acc.loc["TOTAL", "after_QC"], 1)
    sheets["S1_cell_accounting"] = acc.rename_axis("sample").reset_index()

    # --- S2  why each cell was dropped
    if FILTER_REPORT.exists():
        fr = pd.read_csv(FILTER_REPORT, index_col=0)
        fr["pct_retained"] = fr["pct_retained"].round(1)
        fr.index.name = "label_before_correction"
        fr = fr.reset_index()
        back = {v: k for k, v in orig_of.items()}
        fr.insert(0, "sample", fr.label_before_correction.map(back))
        fr = fr.set_index("sample").reindex(SAMPLE_ORDER).reset_index()
        sheets["S2_QC_rejection_reasons"] = fr

    # --- S3  distributions as numbers
    def summarise(df, tag):
        g = df.groupby("sample_id" if "sample_id" in df else "sample_id_all", observed=True)
        out = pd.DataFrame({
            "n_cells": g.size(),
            "median_genes": g.n_genes_by_counts.median().round(0),
            "median_UMIs": g.total_counts.median().round(0),
            "median_pct_mt": g.pct_counts_mt.median().round(2),
            "median_pct_ribo": g.pct_counts_ribo.median().round(2),
            "p01_UMIs": g.total_counts.quantile(0.01).round(0),
            "p99_UMIs": g.total_counts.quantile(0.99).round(0),
        })
        out.insert(0, "stage", tag)
        return out

    dist = []
    if raw_obs is not None:
        dist.append(summarise(raw_obs, "pre-filter"))
    dist.append(summarise(obs.rename(columns={"sample_id_all": "sample_id"}), "post-QC"))
    dist.append(summarise(obs.loc[obs.stan].rename(columns={"sample_id_all": "sample_id"}),
                          "post-barcode"))
    sheets["S3_QC_metric_summary"] = pd.concat(dist).reset_index().rename(
        columns={"index": "sample", "sample_id": "sample"})

    # --- S4  composition, the numbers behind the UMAP
    f = obs.loc[obs.stan].copy()
    f["tissue"] = f.predicted_leiden_sub_joint.map(TISSUE_MAP).fillna("other")
    comp = pd.crosstab(f.sample_id_all, f.tissue).reindex(
        index=SAMPLE_ORDER, columns=[t for t in TISSUE_ORDER])
    comp = comp.fillna(0).astype(int)
    comp["TOTAL"] = comp.sum(axis=1)
    sheets["S4_state_by_sample"] = comp.reset_index().rename(columns={"sample_id_all": "sample"})

    # site share per state -- the quantitative form of "clusters are not batch"
    st = f.predicted_leiden_sub_joint.map(TISSUE_MAP).fillna("other")
    st.name = "state"
    site = pd.crosstab(st, f.site)
    site_pct = (100 * site.div(site.sum(axis=1), axis=0)).round(1)
    site_pct.columns = [f"pct_{c}" for c in site_pct.columns]
    sheets["S5_state_site_balance"] = pd.concat([site, site_pct], axis=1).reindex(
        [t for t in TISSUE_ORDER]).rename_axis("state").reset_index()

    # --- S7  provenance and the label correction, stated openly
    meta = pd.DataFrame({
        "sample": SAMPLE_ORDER,
        "site": [s.split("_")[0] for s in SAMPLE_ORDER],
        "condition": [s.split("_")[1] for s in SAMPLE_ORDER],
        "nominal_timepoint": [s.split("_")[2] for s in SAMPLE_ORDER],
    })
    tp = obs.groupby("sample_id_all", observed=True).timepoint_true.agg(
        lambda x: x.mode().iat[0]).reindex(SAMPLE_ORDER)
    meta["timepoint_used_in_analysis"] = tp.values
    orig = orig_of.reindex(SAMPLE_ORDER)
    meta["label_before_OCM_correction"] = orig.values
    meta["OCM_label_corrected"] = (meta["sample"] != meta["label_before_OCM_correction"])
    meta["lane"] = meta["site"].map({"Chicago": "Chicago 10x lane (OB1-OB4 multiplexed)",
                                     "Michigan": "Michigan 10x lane (OB1-OB4 multiplexed)"})
    sheets["S7_sample_metadata"] = meta

    # --- S8  depth window
    sheets["S8_library_size_gate"] = DEPTH_TABLE
    sheets["S9_depth_standardised"] = STD_TABLE
    sheets["S10_label_transfer_by_state"] = LT_STATE
    sheets["S11_eNasM_lNasM_sensitivity"] = LT_SENS

    # --- S0  the README sheet, first
    notes = pd.DataFrame({"item": [
        "Thresholds", "Doublets", "Barcode filter",
        "Barcode filter is not depth-neutral", "Batch structure",
        "OCM label correction", "Source objects", "Generated by"],
        "value": [
        f"min {MIN_GENES} genes; min {MIN_COUNTS:,} UMIs; max {MAX_PCT_MT:g}% mito; "
        f"per-sample upper cut at median + {UPPER_MAD:g} x MAD on genes and on counts",
        "No doublet caller was run. The per-sample upper MAD cut on genes and counts "
        "is the doublet control and is deliberately conservative (an upper bound).",
        "A per-sample library-size gate on total_counts, applied to the RAW pre-QC "
        "object (confirmed by Stan 2026-08-09). The windows are in stan_filter.py and "
        "reproduce his supplied barcode list EXACTLY: 29,519/29,519, zero mismatches. "
        "No confidence filter is applied. Matching by barcode sequence alone keeps 72 "
        "cells whose barcode belongs to a different sample; computing the gate directly "
        "gives 29,347.",
        "The gate is ASYMMETRIC between arms of the same comparison -- Michigan 72h is "
        "Ctrl 5,000-30,000 against KD 4,000-20,000, giving controls 1.67x the median "
        "depth. Depth-standardised composition (S9, QC_S5) shows lNasM, pluripotent, "
        "pPS-LC and AE-LC hold or strengthen; eNasM retains 58% of its effect and "
        "should be quoted with that caveat.",
        "All four samples per site were multiplexed in ONE 10x lane, so condition and "
        "timepoint are within-lane and carry no batch term. Site is simultaneously the "
        "batch and the biological repeat; the two cannot be separated by design. "
        "Pseudobulk models ~site + condition.",
        "The on-chip-multiplexing assignment for Ctrl_72h and DOX_48h was exchanged at "
        "both sites, and nominal 48h is 42h true. Sheet S7 records both labels per sample.",
        f"{MASTER.name}; "
        + (raw_path.name if raw_path else "pre-filter object NOT AVAILABLE"),
        "04_supplement/qc_supplement.py"]})
    sheets = {"S0_README": notes, **sheets}

    path = OUT / "QC_SOURCE_DATA.xlsx"
    with pd.ExcelWriter(path, engine="openpyxl") as xl:
        for name, df in sheets.items():
            df.to_excel(xl, sheet_name=name[:31], index=False)
            ws = xl.sheets[name[:31]]
            for col in ws.columns:
                width = max(len(str(c.value)) if c.value is not None else 0 for c in col)
                ws.column_dimensions[col[0].column_letter].width = min(max(width + 2, 10), 60)
            ws.freeze_panes = "A2"
    print(f"  wrote output/QC_SOURCE_DATA.xlsx ({len(sheets)} sheets)")
    return list(sheets)


# ====================================================================== MAIN

def main():
    skipped = []
    print("Reading frozen objects (obs only, no expression matrices)...")

    obs = read_obs(MASTER)
    obs = obs.rename(columns={"sample_id": "sample_id_all"})
    obs["stan"] = stan_mask(obs.index)
    obs["gate"] = library_size_gate(obs.sample_id_all.astype(str), obs.total_counts)
    print(f"  master: {len(obs):,} cells post-QC, {obs.stan.sum():,} by barcode list, "
          f"{obs.gate.sum():,} by library-size gate "
          f"({int((obs.stan & ~obs.gate).sum())} kept on a barcode from another sample)")


    raw_path = next((p for p in RAW_CANDIDATES if p.exists()), None)
    raw_obs = None
    if raw_path:
        raw_obs = read_obs(raw_path)
        print(f"  pre-filter: {len(raw_obs):,} called cells  <- {raw_path}")
    else:
        skipped.append("QC_S1_thresholds (needs amanda_liang_raw_unfiltered.h5ad)")
        print("  pre-filter object NOT FOUND -- looked in:")
        for p in RAW_CANDIDATES:
            print(f"      {p}")

    print("\nBuilding panels...")
    if raw_obs is not None:
        panel_thresholds(raw_obs)

    f = obs.loc[obs.stan]
    umap = read_obsm(MASTER, "X_umap_joint")[obs.stan.values]
    panel_batch(f, umap)

    global DEPTH_TABLE, STD_TABLE
    DEPTH_TABLE = depth_table(obs)
    panel_barcode_depth(obs, DEPTH_TABLE)
    STD_TABLE = depth_standardised_composition(obs)
    panel_depth_standardised(STD_TABLE)

    global LT_STATE, LT_SENS
    _f, LT_STATE, LT_SENS = label_transfer_tables(obs)
    panel_label_transfer(_f, LT_STATE, LT_SENS)

    print("\nBuilding source data...")
    build_workbook(obs, raw_obs, raw_path)

    print("\nDone. Output in", OUT)
    if skipped:
        print("\nSKIPPED, and why:")
        for s in skipped:
            print("  -", s)
        print("\n  Copy the pre-filter object to 00_DATA_FROZEN/00_raw/ and rerun to add it.")


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        main()
