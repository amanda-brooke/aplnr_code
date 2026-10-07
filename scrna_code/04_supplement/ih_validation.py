"""
ih_validation.py -- label-transfer validation against the IH gastruloid timeseries.

Do our cells actually sit among IH reference cells carrying the same label?
The joint scVI object holds IH's annotated cells and ours in one latent space, so for
each of our cells we look up its 25 nearest IH neighbours and tabulate their labels.
This is a measured concordance, not the model's own posterior.

Used by supp_figures.py (SF2_3_label_validation).

READS   the joint scVI reference object + the master object (see config.py)
WRITES  ./output/ and ./cache/
RUN     python ih_validation.py
"""
from pathlib import Path
import sys, time, warnings

import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from anndata.io import read_elem

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "utils"))
from config import ROOT
from stan_filter import stan_mask
from palettes import POSITIVE_COLOR, GREY_DARK, GREY_FILL, GREY_MID

OUT = HERE / "output"; OUT.mkdir(exist_ok=True)
CACHE = HERE / "cache"; CACHE.mkdir(exist_ok=True)

from qc_style import apply_style, MATRIX_CMAP, SMALL, TINY
apply_style()


def save(fig, stem):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"{stem}.{ext}")
    plt.close(fig)
    print(f"  wrote output/{stem}  (png, pdf, svg)")



# ============================================ 3. NEIGHBOURHOOD CONCORDANCE

JOINT = (ROOT / "00_DATA_FROZEN" / "01_source" / "models" / "gastruloid_reference_models"
         / "reference_amanda_liang_joint_scvi_model" / "adata.h5ad")
K_NEIGHBOURS = 25
STATE_ORDER = ["pluripotent", "aPS-LC", "pPS-LC", "eNasM", "lNasM", "AE-LC",
               "DE", "PGC-LC", "aPM-LC", "LPM-LC", "CS-LC", "stressed"]


def concordance():
    """Do our cells actually sit among IH cells carrying the same label?

    Every other confidence number in this supplement is the model's own posterior.
    This one is measured: the joint scVI object holds IH's 46,507 annotated cells
    and our 45,435 in ONE 30-dimensional latent, so for each of our cells we can
    look up its nearest IH neighbours and ask what IH calls them.

    Note this is close to re-deriving the kNN transfer that produced the labels, so
    the diagonal is expected to be high. The informative part is the OFF-diagonal --
    where the leakage goes, and whether a state's neighbourhood is the right tissue
    even when the exact sub-stage is ambiguous.
    """
    cache = CACHE / "ih_concordance.npz"
    if cache.exists():
        z = np.load(cache, allow_pickle=True)
        return pd.DataFrame(z["comp"], columns=list(z["labs"])), z["mylab"]

    from stan_filter import library_size_gate
    with h5py.File(JOINT, "r") as f:
        oj = read_elem(f["obs"])
        Z = np.asarray(f["obsm"]["X_scVI_joint"][:], dtype=np.float32)
    with h5py.File(ROOT / "00_DATA_FROZEN" / "03_analysis_ready"
                   / "adata_master_normalized.h5ad", "r") as f:
        om = read_elem(f["obs"])
    om["gate"] = library_size_gate(om.sample_id.astype(str), om.total_counts)

    isref = (oj.dataset == "Reference").values
    R = Z[isref]
    reflab = oj.known_leiden_sub[isref].astype(str).values
    mine = om.reindex(oj.index[~isref].str.replace("^QUERY__", "", regex=True))
    sel = mine.gate.values
    Q = Z[~isref][sel]
    mylab = mine.predicted_leiden_sub_joint.astype(str).values[sel]
    print(f"  {len(Q):,} gated cells vs {len(R):,} annotated IH cells, k={K_NEIGHBOURS}")

    labs = sorted(set(reflab))
    code = np.array([labs.index(l) for l in reflab])
    comp = np.zeros((len(Q), len(labs)), dtype=np.float32)
    Rn = (R ** 2).sum(1)
    for s0 in range(0, len(Q), 1500):                 # chunked brute force; a KD-tree
        q = Q[s0:s0 + 1500]                           # is no faster in 30 dimensions
        d = Rn[None, :] - 2 * q @ R.T                 # monotone in euclidean distance
        nn = np.argpartition(d, K_NEIGHBOURS, axis=1)[:, :K_NEIGHBOURS]
        for i in range(len(q)):
            np.add.at(comp[s0 + i], code[nn[i]], 1.0)
    comp /= K_NEIGHBOURS
    np.savez(cache, comp=comp, mylab=mylab, labs=np.array(labs, dtype=object))
    return pd.DataFrame(comp, columns=labs), mylab


def concordance_tables(comp, mylab):
    C = comp.copy(); C["ours"] = mylab
    mat = C.groupby("ours").mean() * 100
    order = [s for s in STATE_ORDER if s in mat.index]
    mat = mat.reindex(order)[[s for s in STATE_ORDER if s in mat.columns]]

    rows = []
    for st in order:
        d = comp[np.asarray(mylab) == st]
        own = mat.loc[st, st] if st in mat.columns else np.nan
        nas = (d[["eNasM", "lNasM"]].sum(1) * 100).mean()
        rows.append({"our_state": st, "n_cells": int((np.asarray(mylab) == st).sum()),
                     "pct_neighbours_same_label": round(own, 1),
                     "pct_neighbours_nascent_mesoderm": round(nas, 1),
                     "top_confusion": mat.loc[st].drop(st, errors="ignore").idxmax(),
                     "top_confusion_pct": round(mat.loc[st].drop(st, errors="ignore").max(), 1)})
    return mat.round(1).rename_axis("our_state").reset_index(), pd.DataFrame(rows)


def panel_concordance(mat, purity):
    m = mat.set_index("our_state")
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 6.6),
                             gridspec_kw={"width_ratios": [1.55, 1]})

    ax = axes[0]; ax.grid(False)
    im = ax.imshow(m.values, cmap=MATRIX_CMAP, vmin=0, vmax=100, aspect="equal")
    ax.set_xticks(range(m.shape[1])); ax.set_xticklabels(m.columns, rotation=90, fontsize=9)
    ax.set_yticks(range(m.shape[0])); ax.set_yticklabels(m.index, fontsize=9)
    ax.set_xlabel("IH label of the 25 nearest reference cells")
    ax.set_ylabel("our transferred label")
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            v = m.values[i, j]
            if v >= 5:
                ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=TINY,
                        color="white" if v > 45 else GREY_DARK)
    cb = fig.colorbar(im, ax=ax, fraction=0.042, pad=0.03, shrink=0.82)
    cb.set_label("% of neighbours", size=SMALL)
    cb.ax.tick_params(labelsize=TINY)
    ax.set_title("measured concordance, not a model posterior", loc="left")

    ax = axes[1]; ax.grid(False)
    d = purity[purity.our_state.isin(["eNasM", "lNasM"])]
    y = np.arange(len(d))
    for i, r in enumerate(d.itertuples()):
        other = r.pct_neighbours_nascent_mesoderm - r.pct_neighbours_same_label
        rest = 100 - r.pct_neighbours_nascent_mesoderm
        ax.barh(i, r.pct_neighbours_same_label, color=POSITIVE_COLOR, height=0.42,
                label="same label" if i == 0 else None)
        ax.barh(i, other, left=r.pct_neighbours_same_label, color=GREY_MID, height=0.42,
                label="the other nascent state" if i == 0 else None)
        ax.barh(i, rest, left=r.pct_neighbours_nascent_mesoderm, color=GREY_FILL,
                height=0.42, label="anything else" if i == 0 else None)
        ax.text(103, i, f"{r.pct_neighbours_nascent_mesoderm:.0f}%\nnascent\nmesoderm",
                va="center", ha="left", fontsize=9, color=GREY_DARK)
    ax.set_yticks(y); ax.set_yticklabels(d.our_state); ax.invert_yaxis()
    ax.set_ylim(1.75, -0.75)
    ax.set_xlim(0, 132); ax.set_xlabel("% of 25 nearest IH neighbours")
    ax.legend(frameon=False, fontsize=SMALL, loc="upper center",
              bbox_to_anchor=(0.45, -0.14), ncol=3, columnspacing=1.2)
    ax.set_title("identity is solid; the sub-stage boundary\nis where the uncertainty lives",
                 loc="left")

    fig.suptitle("Our cells sit among IH cells of the same tissue. eNasM and lNasM "
                 "exchange with each other and with almost nothing else.",
                 x=0.01, y=0.995, va="top", ha="left", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.83))
    save(fig, "QC_S9_IH_neighbourhood_concordance")


def main():
    t0 = time.time()
    print("Neighbourhood concordance with IH")
    comp, mylab = concordance()
    mat, purity = concordance_tables(comp, mylab)
    print(purity.to_string(index=False))
    panel_concordance(mat, purity)
    path = OUT / "IH_VALIDATION_SOURCE_DATA.xlsx"
    with pd.ExcelWriter(path, engine="openpyxl") as xl:
        mat.to_excel(xl, sheet_name="IH_concordance_matrix", index=False)
        purity.to_excel(xl, sheet_name="IH_neighbourhood_purity", index=False)
    print(f"  wrote output/{path.name}")
    print(f"\nDone in {time.time()-t0:.0f}s.")


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        main()
