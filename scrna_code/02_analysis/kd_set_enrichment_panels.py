"""
kd_set_enrichment_panels.py -- draws the two panels for ARMS 3 and 3.5.

Panel G3  "does the APLNR-correlated set depend on APLNR?"
          Null distribution of the set's mean log2FC under 5,000 expression-matched
          random gene sets, with the observed value marked, in TWO compartments:
            NasM      -- receptor compartment, APLNR in 75-83% of cells
            pluri_AE  -- ligand compartment, the specificity control
          The point of showing both is the SIGN FLIP. The set falls where the receptor
          is and rises where it isn't. A null result in the control would be weaker,
          because a null result is also what low power looks like.

Panel G4  "which programmes fall, APLNR-linked or not?"
          z per pre-specified category in NasM. Controls drawn in grey and labelled as
          controls, so a reader can verify the artefact check without reading methods.

Both read from kd_set_enrichment.py (same functions, same seed) so the figures and
KD_set_enrichment.csv cannot drift apart.

Run: python kd_set_enrichment_panels.py     (~40 s; recomputes the nulls to draw them)
Writes SVG + PDF into 07_figures/FIGURE_KD_maturation/.
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.patches import Patch

import kd_set_enrichment as SE

ROOT = SE.ROOT
TAB  = SE.TAB
FIG  = ROOT / "07_figures" / "FIGURE_KD_maturation"; FIG.mkdir(parents=True, exist_ok=True)

SAVE_FORMATS = ["svg", "pdf"]
CTRL_C, KD_C = "#BE6455", "#D99A90"     # house palette, matches maturation_kd.ipynb
NULL_C       = "#C9C4BD"
GREY         = "#8A8780"
INK          = "#2B2926"

APLNR_SET = "APLNR-correlated (up_both, top 200)"
CONTROLS  = {"pluripotency (negative control)", "ribosome (artefact control)"}


def save(fig, stem):
    for e in SAVE_FORMATS:
        fig.savefig(FIG / f"{stem}.{e}", bbox_inches="tight")
    print(f"  wrote {stem}.{'/'.join(SAVE_FORMATS)}")


def build(tissue, rng, res, sets):
    """Recompute matched nulls for one tissue, keeping the full null draws."""
    gcol = "gene" if "gene" in res.columns else res.columns[0]
    d = res[res.tissue == tissue].set_index(gcol)[["log2FC_shrunk", "baseMean"]].dropna()
    d = d[~d.index.duplicated()]
    d["bin"] = pd.qcut(np.log10(d["baseMean"] + 1), SE.N_BINS, labels=False,
                       duplicates="drop")
    out = {}
    for name, genes in sets.items():
        r = SE.matched_null(d, genes, SE.N_PERM, rng)
        if r:
            out[name] = r
    return out


def panel_G3(by_tissue, show=False):
    """Null distributions for the APLNR-correlated set, receptor vs ligand compartment."""
    tissues = [("NasM", "NasM (72h) — receptor compartment\nAPLNR in 75–83% of cells"),
               ("pluri_AE", "pluri_AE (72h) — ligand compartment\nspecificity control")]
    fig, axes = plt.subplots(2, 1, figsize=(7.4, 6.4), sharex=True)

    for ax, (t, title) in zip(axes, tissues):
        r = by_tissue.get(t, {}).get(APLNR_SET)
        if r is None:
            ax.text(.5, .5, f"{APLNR_SET}\nnot fitted in {t}", ha="center",
                    va="center", transform=ax.transAxes, color=GREY)
            continue
        null = r["_null_means"]; obs = r["mean_log2FC"]
        ax.hist(null, bins=60, color=NULL_C, edgecolor="white", linewidth=.4)
        ax.axvline(null.mean(), color=GREY, lw=1.2, ls="--", zorder=3)
        col = CTRL_C if obs < null.mean() else "#3B7A9E"
        ax.axvline(obs, color=col, lw=2.6, zorder=4)

        # arrow from the null centre to the observed value makes the direction explicit
        y = ax.get_ylim()[1]
        ax.annotate("", xy=(obs, y * .78), xytext=(null.mean(), y * .78),
                    arrowprops=dict(arrowstyle="-|>", color=col, lw=1.6,
                                    shrinkA=0, shrinkB=0))
        # Put the label in the EMPTY GAP between the observed line and the histogram,
        # never on the outboard side -- the observed value sits near the axis limit and
        # an outboard label runs off the figure.
        pad = 0.02 * (ax.get_xlim()[1] - ax.get_xlim()[0])
        left_of_null = obs < null.mean()
        ax.text(obs + pad if left_of_null else obs - pad, y * .92,
                f"observed {obs:+.3f}\nz = {r['z']:+.2f}",
                ha="left" if left_of_null else "right", va="top",
                color=col, fontsize=10, fontweight="bold")
        # white bbox: this label sits over the histogram bars and is unreadable without it
        ax.text(null.mean(), y * .30, f"null {null.mean():+.3f}",
                ha="right" if left_of_null else "left", va="center",
                color=GREY, fontsize=8.5, zorder=5,
                bbox=dict(facecolor="white", edgecolor="none", alpha=.85,
                          boxstyle="round,pad=0.22"))
        ax.set_title(title, fontsize=10.5, loc="left", color=INK)
        ax.set_ylabel("random gene sets")
        for sp in ["top", "right"]:
            ax.spines[sp].set_visible(False)

    axes[-1].set_xlabel("mean pseudobulk log$_2$FC (KD / Ctrl) of the gene set")
    fig.suptitle("The APLNR-correlated programme falls where the receptor is,\n"
                 "and rises where it is not",
                 fontsize=12.5, y=1.005, ha="left", x=.02, color=INK)
    fig.text(.02, -.035,
             f"{SE.N_PERM:,} random gene sets per panel, matched to the observed set's "
             f"expression (baseMean) decile distribution.\n"
             f"Set = {APLNR_SET.replace('(', '').replace(')', '')}, defined by "
             f"correlation in CONTROL cells only, then tested on the KD contrast.",
             fontsize=7.6, color=GREY, ha="left")
    fig.tight_layout()
    save(fig, "G3_aplnr_set_enrichment_null")
    if not show:
        plt.close(fig)
    return fig


def panel_G4(by_tissue, enr, show=False):
    """z per pre-specified category, NasM."""
    rows = []
    for name, r in by_tissue.get("NasM", {}).items():
        rows.append({"set": name, "z": r["z"], "n": r["n_genes"]})
    d = pd.DataFrame(rows)
    padj = (enr[enr.tissue == "NasM"].set_index("set")["padj"])
    d["padj"] = d["set"].map(padj)
    d = d.sort_values("z")

    fig, ax = plt.subplots(figsize=(8.2, 5.6))
    ax.axvspan(-2, 2, color="#F2F0ED", zorder=0)
    ax.axvline(0, color=GREY, lw=1.1, zorder=1)

    for i, r in enumerate(d.itertuples()):
        ctrl = r.set in CONTROLS
        sig = (r.padj < 0.05) and not ctrl
        col = GREY if ctrl else (CTRL_C if sig else KD_C)
        ax.plot([0, r.z], [i, i], color=col, lw=2.2, alpha=.9, zorder=2)
        ax.scatter(r.z, i, s=74, color=col, zorder=3,
                   edgecolor="white", linewidth=1.1)
        lab = f"{r.set}  (n={r.n})"
        ax.text(r.z - .28 if r.z < 0 else r.z + .28, i,
                f"padj {r.padj:.3f}" if r.padj >= .001 else "padj <0.001",
                va="center", ha="right" if r.z < 0 else "left",
                fontsize=7.8, color=col)

    # Italicise the control rows via label properties, NOT mathtext: $\it{a b}$ drops
    # the spaces and renders "artefactcontrol".
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels(list(d["set"]), fontsize=9.5)
    for lbl, name in zip(ax.get_yticklabels(), d["set"]):
        if name in CONTROLS:
            lbl.set_style("italic"); lbl.set_color(GREY)
    ax.set_xlabel("z vs expression-matched random gene sets   "
                  "(negative = programme falls in KD)")
    ax.set_title("Which programmes fall in the APLNR knockdown?   NasM, 72h",
                 fontsize=12, loc="left", color=INK)
    # generous left margin: the padj label for the largest effect otherwise collides
    # with the category name in the y tick labels
    ax.set_xlim(min(-14.5, d.z.min() - 3.4), max(3.4, d.z.max() + 2.2))
    for sp in ["top", "right", "left"]:
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis="y", length=0)

    # anchored outside the axes: "lower left" sat on top of the two largest effects
    ax.legend(handles=[Patch(color=CTRL_C, label="falls, padj < 0.05"),
                       Patch(color=KD_C,  label="not significant"),
                       Patch(color=GREY,  label="control set (italic)")],
              frameon=False, fontsize=8.4, loc="upper left",
              bbox_to_anchor=(0.005, -0.105), ncol=3, handlelength=1.2,
              columnspacing=1.4)
    fig.text(.01, -.155,
             "Grey band is |z| < 2. Controls behave as required: the ribosome set is UP, "
             "not down, so the\ncompositional shift (KD cells are 25–38% less complex) is "
             "not driving these results; pluripotency is flat.\n"
             "Categories are pre-specified from textbook membership — see kd_set_enrichment.py.",
             fontsize=7.6, color=GREY, ha="left")
    fig.tight_layout()
    save(fig, "G4_programme_enrichment_NasM")
    if not show:
        plt.close(fig)
    return fig


def main(show=False):
    plt.rcParams.update({"font.size": 10, "svg.fonttype": "none",
                         "pdf.fonttype": 42, "axes.edgecolor": INK,
                         "axes.labelcolor": INK, "text.color": INK,
                         "xtick.color": INK, "ytick.color": INK})
    res = pd.read_csv(TAB / "KD_pseudobulk_results.csv")
    enr = pd.read_csv(TAB / "KD_set_enrichment.csv")

    sets = dict(SE.CATEGORIES)
    p40 = TAB / "genes_rising_with_APLNR.csv"
    if p40.exists():
        sets["APLNR-correlated (top 40)"] = [x for x in pd.read_csv(p40)["gene"]
                                             if x != "APLNR"]
    pbr = TAB / "aplnr_percell_correlation_byrepeat.csv"
    if pbr.exists():
        c = pd.read_csv(pbr)
        rho = "rho_raw_mean" if "rho_raw_mean" in c.columns else "mean"
        up = c[c.up_both & ~c.gene.astype(str).str.startswith("ENSG")]
        sets[APLNR_SET] = [x for x in up.sort_values(rho, ascending=False).gene
                           if x != "APLNR"][:200]

    by_tissue = {}
    for t in ["NasM", "pluri_AE"]:
        by_tissue[t] = build(t, np.random.default_rng(SE.SEED), res, sets)

    panel_G3(by_tissue, show=show)
    panel_G4(by_tissue, enr, show=show)
    print(f"\nfigures -> {FIG}")


if __name__ == "__main__":
    main()
