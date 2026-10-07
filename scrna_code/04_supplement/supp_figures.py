"""
supp_figures.py -- assembles the loose QC panels into the four supplementary
                   figures that actually go in the paper.

WHY THIS IS SEPARATE FROM THE ANALYSIS SCRIPTS
qc_supplement.py and ih_validation.py compute things and emit one panel per
question, which is the right shape while you are still deciding what is true.
This file makes no measurements at all -- it imports those functions and only
decides layout. So a change to the science happens in one place and a change to
the figure happens here, and they never argue.

WHAT WAS DROPPED, AND WHY
  %ribosomal      never used as a filter -- decoration, not evidence
  QC_S3           was the label-transfer histogram; QC_S9 measures the same thing
                  against IH's real annotations instead of the model's posterior
  QC_S6a,b        superseded by QC_S9 for the same reason; S6c survives in SF5.1
  S4 right panel  the p99/p01 bars restated the brackets; numbers live in sheet S8

OUTPUT   output/SF2_1_cell_filtering.*   SF2_2_replicate_concordance.*
         output/SF2_3_label_validation.* SF5_1_KD_controls.*
"""
from pathlib import Path
import warnings

import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec

import qc_supplement as qc
import ih_validation as ihv
from qc_style import apply_style, MATRIX_CMAP, SMALL, TINY, BASE as BASE_FS
from palettes import (TISSUE_COLORS, TISSUE_MAP, TISSUE_ORDER, SITE_COLORS,
                      PROGRAM_COLORS,
                      TIMEPOINT_COLORS, CONDITION_COLORS, CTRL_C, KD_C,
                      STATE_COLORS, POSITIVE_COLOR, GREY_DARK, GREY_FILL,
                      GREY_MID, GREY_LINE, RED, TREND_COLOR)

apply_style()
OUT = Path(__file__).resolve().parent / "output"

DOT = 5           # UMAP marker size -- matches the scatter sizes in figure2_clean
DOT_ALPHA = 0.55
PANEL_LETTER = dict(fontsize=15, fontweight="bold", va="top", ha="left")


def letter(fig, ch, x, y):
    fig.text(x, y, ch, **PANEL_LETTER)


def save(fig, stem):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"{stem}.{ext}")
    plt.close(fig)
    print(f"  wrote output/{stem}  (png, pdf, svg)")


def load():
    """Everything the four figures need, computed once."""
    obs = qc.read_obs(qc.MASTER).rename(columns={"sample_id": "sample_id_all"})
    obs["stan"] = qc.stan_mask(obs.index)
    obs["gate"] = qc.library_size_gate(obs.sample_id_all.astype(str), obs.total_counts)
    raw_path = next((p for p in qc.RAW_CANDIDATES if p.exists()), None)
    raw_obs = qc.read_obs(raw_path) if raw_path else None
    umap = qc.read_obsm(qc.MASTER, "X_umap_joint")
    return obs, raw_obs, umap


# =============================================== SF2.1  CELL FILTERING

def sf2_1(raw_obs, obs):
    """Every filter that was applied, on one shared sample axis."""
    rows = [("n_genes_by_counts", "Genes per cell", qc.MIN_GENES, True),
            ("total_counts", "UMIs per cell", qc.MIN_COUNTS, True),
            ("pct_counts_mt", "% mitochondrial", qc.MAX_PCT_MT, False)]
    order = [s for s in qc.SAMPLE_ORDER if s in set(raw_obs.sample_id)]

    fig, axes = plt.subplots(5, 1, figsize=(9.6, 13.6), sharex=True,
                             gridspec_kw={"height_ratios": [1, 1, 1, 1, 0.9]})

    for ax, (col, label, floor, log) in zip(axes[:3], rows):
        ax.grid(False)
        raw = [raw_obs.loc[raw_obs.sample_id == s, col].astype(float).values for s in order]
        data = [np.log10(np.clip(v, 1, None)) for v in raw] if log else raw
        parts = ax.violinplot(data, positions=range(len(order)), widths=0.8,
                              showextrema=False, showmedians=False)
        for pc, s in zip(parts["bodies"], order):
            pc.set_facecolor(SITE_COLORS[s.split("_")[0]])
            pc.set_alpha(0.65); pc.set_edgecolor("none")
        for i, d in enumerate(data):
            ax.plot([i - 0.2, i + 0.2], [np.median(d)] * 2, color=GREY_DARK, lw=1.3, zorder=3)
        tx = (lambda y: np.log10(max(y, 1))) if log else (lambda y: y)
        ax.axhline(tx(floor), color=RED, lw=1.1, ls="--", zorder=4)
        ax.text(len(order) - 0.45, tx(floor), f"  {floor:,g}", color=RED,
                va="bottom", ha="left", fontsize=SMALL)
        if col in ("n_genes_by_counts", "total_counts"):
            for i, smp in enumerate(order):
                v = raw_obs.loc[raw_obs.sample_id == smp, col].astype(float)
                ax.plot([i - 0.35, i + 0.35],
                        [tx(np.median(v) + qc.UPPER_MAD * qc.mad(v))] * 2,
                        color=RED, lw=1.1, zorder=4)
        if log:
            lo, hi = ax.get_ylim()
            t = np.arange(np.floor(lo), np.ceil(hi) + 1)
            ax.set_yticks(t); ax.set_yticklabels([f"$10^{{{int(x)}}}$" for x in t])
        ax.set_ylabel(label)

    # row 4 -- the library-size gate, same axis units as row 2
    ax = axes[3]; ax.grid(False)
    for i, smp in enumerate(order):
        m = obs.sample_id_all == smp
        a = obs.loc[m, "total_counts"].astype(float)
        b = obs.loc[m & obs.gate, "total_counts"].astype(float)
        for off, v, c in [(-0.16, a, GREY_FILL),
                          (0.16, b, SITE_COLORS[smp.split("_")[0]])]:
            lo, hi = np.percentile(v, [1, 99])
            ax.plot([i + off] * 2, [lo, hi], color=c, lw=7, solid_capstyle="butt", zorder=2)
            ax.plot([i + off - 0.1, i + off + 0.1], [np.median(v)] * 2,
                    color=GREY_DARK, lw=1.2, zorder=3)
        lo, hi = qc.LIBRARY_SIZE_GATE[smp]
        ax.plot([i + 0.40] * 2, [lo, hi], color=RED, lw=1.1, zorder=5)
        for y in (lo, hi):
            ax.plot([i + 0.33, i + 0.47], [y] * 2, color=RED, lw=1.1, zorder=5)
    ax.set_yscale("log"); ax.set_ylabel("UMIs per cell\nafter QC (grey) and gate")

    # row 5 -- the counts, which is what "cells per sample" actually asks for
    ax = axes[4]; ax.grid(False)
    orig = (obs.groupby("sample_id_all", observed=True).sample_id_original
               .agg(lambda x: x.mode().iat[0]).astype(str))
    called = raw_obs.groupby("sample_id", observed=True).size()
    stages = pd.DataFrame({
        "called":  orig.map(called),
        "after QC": obs.groupby("sample_id_all", observed=True).size(),
        "after gate": obs[obs.gate].groupby("sample_id_all", observed=True).size(),
    }).reindex(order).fillna(0)
    x = np.arange(len(order))
    for k, (col, c) in enumerate([("called", GREY_FILL), ("after QC", GREY_MID),
                                  ("after gate", GREY_DARK)]):
        ax.bar(x + (k - 1) * 0.27, stages[col], 0.26, color=c,
               edgecolor=GREY_DARK, lw=0.4, label=col)
    for i, smp in enumerate(order):
        ax.text(i, stages.loc[smp, "called"] * 1.06, f"{int(stages.loc[smp,'after gate']):,}",
                ha="center", fontsize=TINY, color=GREY_DARK)
    ax.set_ylabel("cells")
    ax.set_yscale("log")
    ax.legend(frameon=False, fontsize=SMALL, ncol=3, loc="lower left")
    ax.set_title(f"Cells retained at each stage. Total "
                 f"{int(stages['called'].sum()):,} called \u2192 "
                 f"{int(stages['after QC'].sum()):,} after QC \u2192 "
                 f"{int(stages['after gate'].sum()):,} after gate.",
                 loc="left", fontsize=SMALL)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=90)

    axes[0].set_title(f"Pre-filter QC, all {len(raw_obs):,} called cells.  dashed = fixed "
                      f"threshold, solid red = per-sample median + {qc.UPPER_MAD:g}xMAD "
                      f"(also the doublet control), grey = median", loc="left", fontsize=SMALL)
    axes[3].set_title("Library-size gate (red brackets). Windows differ between arms of the "
                      "same comparison.", loc="left", fontsize=SMALL)
    fig.tight_layout()
    for ch, y in zip("abcde", [0.995, 0.800, 0.605, 0.410, 0.205]):
        letter(fig, ch, 0.005, y)
    save(fig, "SF2_1_cell_filtering")


# ============================================ SF2.2  DATASET STRUCTURE

def umap_panel(ax, xy, vals, cmap, title, levels=None, ncol=1):
    """One scatter call in shuffled order.

    Plotting one level at a time means the last group drawn sits on top of every
    other -- in the condition panel that made KD look like it covered the embedding
    simply because DOX was drawn second. Colours are assigned per point instead and
    the legend is built from proxy handles.
    """
    from matplotlib.lines import Line2D
    ax.grid(False)
    levels = levels or sorted(set(vals))
    counts = {lv: int((vals == lv).sum()) for lv in levels}
    colours = np.array([cmap.get(v, GREY_MID) for v in vals], dtype=object)
    ax.scatter(xy[:, 0], xy[:, 1], s=DOT, lw=0, alpha=DOT_ALPHA,
               c=list(colours), rasterized=True)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", ms=6, color=cmap.get(lv, GREY_MID),
                              label=f"{lv} ({counts[lv]:,})") for lv in levels],
              frameon=False, loc="upper left", handletextpad=0.1, labelspacing=0.3,
              borderpad=0.2, fontsize=SMALL, ncol=ncol)
    ax.set_title(title, loc="left")
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal", adjustable="datalim")
    for sp in ax.spines.values():
        sp.set_visible(False)


def sf2_2(obs, umap):
    """One UMAP per biological replicate, same colouring.

    The main figure already shows the embedding by timepoint and by label, so the
    only job left for the supplement is: do the two repeats give the same thing?
    An overlaid two-colour UMAP answers that badly -- whichever group is drawn
    second looks dominant, which is exactly the artefact that made the condition
    panel misleading before. Facetting by replicate, with the full dataset in grey
    behind, makes co-occupancy something you can see rather than infer.
    """
    f = obs[obs.gate].copy()
    f["tissue"] = f.predicted_leiden_sub_joint.map(TISSUE_MAP).fillna("other")
    xy = umap[obs.gate.values]
    sites = ["Chicago", "Michigan"]
    levels = [t for t in TISSUE_ORDER if t in set(f.tissue)]

    fig, axes = plt.subplots(1, 3, figsize=(16.4, 5.4),
                             gridspec_kw={"width_ratios": [1, 1, 0.95]})
    for ax, site in zip(axes[:2], sites):
        ax.grid(False)
        m = (f.site.astype(str) == site).values
        ax.scatter(xy[~m, 0], xy[~m, 1], s=DOT * 0.7, lw=0, alpha=0.30,
                   color=GREY_FILL, rasterized=True)
        sub = f[m]
        cols = np.array([TISSUE_COLORS.get(t, GREY_MID) for t in sub.tissue], dtype=object)
        ax.scatter(xy[m, 0], xy[m, 1], s=DOT, lw=0, alpha=DOT_ALPHA,
                   c=list(cols), rasterized=True)
        ax.set_title(f"replicate {sites.index(site)+1} \u2014 {site}  (n = {m.sum():,})",
                     loc="left")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_aspect("equal", adjustable="datalim")
        for sp in ax.spines.values():
            sp.set_visible(False)

    # (c) the quantitative version -- "concordant across repeats" as numbers
    ax = axes[2]; ax.grid(False)
    frac = {}
    for site in sites:
        sub = f[f.site.astype(str) == site]
        frac[site] = (sub.tissue.value_counts(normalize=True)
                        .reindex(levels).fillna(0) * 100)
    y = np.arange(len(levels))
    ax.barh(y - 0.2, frac[sites[0]].values, height=0.38, color=SITE_COLORS[sites[0]],
            edgecolor=GREY_DARK, lw=0.4, label=f"{sites[0]} (rep 1)")
    ax.barh(y + 0.2, frac[sites[1]].values, height=0.38, color=SITE_COLORS[sites[1]],
            edgecolor=GREY_DARK, lw=0.4, label=f"{sites[1]} (rep 2)")
    ax.set_yticks(y); ax.set_yticklabels(levels, fontsize=SMALL); ax.invert_yaxis()
    ax.set_xlabel("% of that replicate's cells")
    ax.legend(frameon=False, fontsize=SMALL, loc="lower right")
    rho = np.corrcoef(frac[sites[0]].values, frac[sites[1]].values)[0, 1]
    ax.set_title(f"same states, different proportions\n"
                 f"(rank agreement, Pearson r = {rho:.2f})", loc="left")

    from matplotlib.lines import Line2D
    fig.legend(handles=[Line2D([], [], marker="o", ls="", ms=6,
                               color=TISSUE_COLORS.get(t, GREY_MID), label=t)
                        for t in levels],
               frameon=False, loc="center left", bbox_to_anchor=(0.995, 0.5),
               fontsize=SMALL)
    fig.suptitle(f"scVI joint embedding, {len(f):,} cells. Grey = all cells; coloured = "
                 f"that replicate. Each site is one 10x lane carrying all four samples, "
                 f"so site is simultaneously batch and biological repeat.",
                 x=0.005, y=0.995, va="top", ha="left", fontsize=SMALL)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    for ch, x in zip("abc", [0.005, 0.335, 0.665]):
        letter(fig, ch, x, 0.88)
    save(fig, "SF2_2_replicate_concordance")


# ============================================== SF2.3  VALIDATION vs IH

def sf2_3():
    """Do our eNasM / lNasM labels mean what IH's mean?"""
    comp, mylab = ihv.concordance()
    mat, purity = ihv.concordance_tables(comp, mylab)
    m = mat.set_index("our_state")

    fig, axes = plt.subplots(1, 3, figsize=(16.8, 6.2),
                             gridspec_kw={"width_ratios": [1.45, 0.95, 0.85]})

    ax = axes[0]; ax.grid(False)
    im = ax.imshow(m.values, cmap=MATRIX_CMAP, vmin=0, vmax=100, aspect="equal")
    ax.set_xticks(range(m.shape[1])); ax.set_xticklabels(m.columns, rotation=90, fontsize=SMALL)
    ax.set_yticks(range(m.shape[0])); ax.set_yticklabels(m.index, fontsize=SMALL)
    ax.set_xlabel("IH label of the 25 nearest reference cells")
    ax.set_ylabel("our transferred label")
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            v = m.values[i, j]
            if v >= 5:
                ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=TINY,
                        color="white" if v > 45 else GREY_DARK)
    cb = fig.colorbar(im, ax=ax, fraction=0.042, pad=0.03, shrink=0.8)
    cb.set_label("% of neighbours", size=SMALL); cb.ax.tick_params(labelsize=TINY)
    ax.set_title("measured concordance, not a model posterior", loc="left")

    ax = axes[1]; ax.grid(False)
    d = purity[purity.our_state.isin(["eNasM", "lNasM"])]
    for i, r in enumerate(d.itertuples()):
        other = r.pct_neighbours_nascent_mesoderm - r.pct_neighbours_same_label
        ax.barh(i, r.pct_neighbours_same_label, color=POSITIVE_COLOR, height=0.42,
                label="same label" if i == 0 else None)
        ax.barh(i, other, left=r.pct_neighbours_same_label, color=GREY_MID, height=0.42,
                label="the other nascent state" if i == 0 else None)
        ax.barh(i, 100 - r.pct_neighbours_nascent_mesoderm,
                left=r.pct_neighbours_nascent_mesoderm, color=GREY_FILL, height=0.42,
                label="anything else" if i == 0 else None)
        ax.text(103, i, f"{r.pct_neighbours_nascent_mesoderm:.0f}% nascent mesoderm",
                va="center", ha="left", fontsize=SMALL, color=GREY_DARK)
    ax.set_yticks(range(len(d))); ax.set_yticklabels(d.our_state)
    ax.set_ylim(1.7, -0.7); ax.set_xlim(0, 152)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("% of 25 nearest IH neighbours")
    ax.legend(frameon=False, fontsize=TINY, loc="upper center",
              bbox_to_anchor=(0.42, -0.16), ncol=3, columnspacing=1.2)
    ax.set_title("identity is solid; the sub-stage boundary is not", loc="left")

    # (c) the transfer posterior itself, per state -- part of the method
    ax = axes[2]; ax.grid(False)
    obs_all = qc.read_obs(qc.MASTER).rename(columns={"sample_id": "sample_id_all"})
    obs_all["gate"] = qc.library_size_gate(obs_all.sample_id_all.astype(str),
                                           obs_all.total_counts)
    ff = obs_all[obs_all.gate]
    st_order = [t for t in TISSUE_ORDER if t in set(ff.predicted_leiden_sub_joint)]
    for i, stt in enumerate(st_order):
        v = ff.loc[ff.predicted_leiden_sub_joint == stt,
                   "prediction_confidence_joint"].astype(float)
        lo, md, hi = np.percentile(v, [25, 50, 75])
        ax.plot([lo, hi], [i] * 2, color=TISSUE_COLORS.get(stt, GREY_MID), lw=4,
                solid_capstyle="butt")
        ax.plot(md, i, "o", color=GREY_DARK, ms=4, zorder=3)
    ax.axvline(0.5, color=RED, lw=1.0, ls="--")
    ax.set_yticks(range(len(st_order))); ax.set_yticklabels(st_order, fontsize=SMALL)
    ax.invert_yaxis(); ax.set_xlim(0, 1)
    ax.set_xlabel("IH label-transfer probability")
    ax.set_title("median and IQR per state\n(no confidence filter is applied)", loc="left")

    fig.suptitle("Our cells sit among IH cells of the same tissue. eNasM and lNasM "
                 "exchange with each other and with almost nothing else, and are the "
                 "least confident of the well-populated states.",
                 x=0.005, y=0.995, va="top", ha="left", fontsize=SMALL)
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    for ch, x in zip("abc", [0.005, 0.47, 0.755]):
        letter(fig, ch, x, 0.91)
    save(fig, "SF2_3_label_validation")


# ================================================ SF5.1  KD CONTROLS

def sf5_1(obs, umap):
    """Controls for the Ctrl-vs-KD comparison itself: batch, depth, label noise."""
    f = obs[obs.gate].copy()
    xy = umap[obs.gate.values]
    rng = np.random.default_rng(0); sh = rng.permutation(len(f))
    std = qc.depth_standardised_composition(obs)
    _, _, sens = qc.label_transfer_tables(obs)

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.8))

    umap_panel(axes[0], xy[sh], f.condition.astype(str).values[sh], CONDITION_COLORS,
               "Condition \u2014 the composition result, not a batch panel")

    ax = axes[1]; ax.grid(False)
    y = np.arange(len(std))
    ax.barh(y - 0.2, std.diff_raw, height=0.38, color=GREY_FILL, edgecolor=GREY_DARK,
            lw=0.5, label="as analysed")
    ax.barh(y + 0.2, std.diff_standardised, height=0.38, color=GREY_LINE,
            edgecolor=GREY_DARK, lw=0.5, label="depth-standardised")
    ax.axvline(0, color=GREY_DARK, lw=0.9)
    ax.set_yticks(y); ax.set_yticklabels(std.state); ax.invert_yaxis()
    ax.set_xlabel("KD - Ctrl, % of cells at 72 h")
    ax.legend(frameon=False, loc="lower right")
    ax.set_title("Not sequencing depth\n(eNasM retains 58%, the rest hold)", loc="left")

    ax = axes[2]; ax.grid(False)
    d = sens[sens.state == "lNasM"]
    x = np.arange(len(d))
    ax.bar(x - 0.19, d.Ctrl, 0.36, color=CTRL_C, label="Ctrl")
    ax.bar(x + 0.19, d.KD, 0.36, color=KD_C, label="KD")
    for i, r in enumerate(d.itertuples()):
        ax.text(i, max(r.Ctrl, r.KD) * 1.03, f"1 : {r.KD_per_Ctrl:.2f}", ha="center",
                fontsize=SMALL, color=GREY_DARK)
    ax.set_xticks(x); ax.set_xticklabels(["as analysed", "confidence\nfiltered"])
    ax.set_ylabel("lNasM cells at 72 h")
    ax.legend(frameon=False)
    ax.set_title("Not label noise\n(ratio unchanged)", loc="left")

    fig.tight_layout()
    for ch, x_ in zip("abc", [0.005, 0.345, 0.675]):
        letter(fig, ch, x_, 0.995)
    save(fig, "SF5_1_KD_controls")


def main():
    print("loading...")
    obs, raw_obs, umap = load()
    print(f"  {len(obs):,} cells, {obs.gate.sum():,} after the gate")
    if raw_obs is not None:
        sf2_1(raw_obs, obs)
    else:
        print("  SKIPPED SF2.1 -- pre-filter object not found")
    sf2_2(obs, umap)
    sf2_3()
    sf5_1(obs, umap)
    print("\nDone.")


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        main()
