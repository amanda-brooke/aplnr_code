"""
kd_aplnr_correlation.py -- ARM 2 producer.

Writes `aplnr_percell_correlation_byrepeat.csv`, which several cells of
maturation_kd.ipynb read but which NO script produced. It existed only as a file
on disk, so its provenance -- in particular whether it was computed on MAGIC-imputed
or raw values -- could not be checked. That is the same class of problem as the
bare-filename reads: the analysis was not reconstructible from the repo.

WHAT IT COMPUTES
Per-cell Spearman correlation between APLNR and every other gene, in CONTROL
mesoderm at 72h, separately for each biological repeat (Chicago, Michigan).

Three flavours are written:

  rho_raw_{site}    Spearman on log-normalised RAW counts.  <-- USE THIS ONE
  rho_mgc_{site}    MAGIC-imputed APLNR (the obs column) vs RAW partner gene.
  rho_adj_{site}    raw, after regressing log1p(library_size) and
                    log1p(genes_per_cell) out of BOTH gene and APLNR.

SETTLED: THERE IS NO GENOME-WIDE MAGIC IN THIS DATASET
The master is called `adata_master_normalized_magic.h5ad`, which invites the
assumption that `.X` is MAGIC-imputed. It is not. Checked directly: `.X` and
`layers["counts"]` have identical nnz (163,204,404) and identical sparsity (9.3%),
and `.X` values run 0.14-5.12. That is log1p(CP10K) on the same sparsity pattern as
the raw counts. Genome-wide MAGIC output would be dense with essentially no zeros.

MAGIC in this project exists ONLY as two obs columns, `APLNR_magic` and
`APELA_magic`. So any genome-wide correlation is already raw, and the earlier worry
that arm 2 rested on imputed values was unfounded.

What DOES remain is an asymmetry: several notebook cells correlate the *smoothed*
APLNR (obs column) against *unsmoothed* partner genes. Denoising one side of a
correlation inflates it. `rho_mgc_` quantifies exactly that, so the size of the
effect is visible rather than argued about. Measured previously at the all-mesoderm
level it was negligible (0.391 vs 0.389), but check it per site rather than trusting
that one number.

DO NOT rename the master to fix this. Notebooks reference it by name. Fix the
misleading impression in documentation instead, which ANALYSIS_AUDIT.md now does.

WHY THE DEPTH-ADJUSTED COLUMN EXISTS
Two genes both detected more often in deeply sequenced cells will correlate even if
they are biologically unrelated. KD cells are 25-38% less complex than controls in
this dataset, so depth is not a nuisance here -- it is structured. rho_adj is the
honest version. Where rho_raw and rho_adj disagree, believe rho_adj.

REPRODUCIBILITY FLAGS (`up_both` / `down_both`) are computed on rho_raw and require
the SAME SIGN in both repeats AND |rho| >= MIN_RHO in both, so nothing enters the
gene list on the strength of one site.

*** A CAVEAT THAT MUST TRAVEL WITH ANY FIGURE BUILT ON THIS FILE ***
Chicago control mesoderm at 72h is ~267 cells with APLNR detected in ~22% of them.
A Spearman on 267 cells that are 78% zero is driven mostly by detected-vs-not, and
its confidence interval is wide. Michigan (~2,130 cells) carries the estimate. The
both-repeats requirement is therefore a check on SIGN, not an independent
replication of magnitude. Do not quote Chicago rho values to two decimals.

Run in the pseudobulk env:  python kd_aplnr_correlation.py
Needs: numpy, pandas, h5py, scipy.  Runtime ~2 min, dominated by loading the master.
"""
import numpy as np, pandas as pd, h5py, scipy.sparse as sp, time
from pathlib import Path
from scipy.stats import rankdata

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
from config import ROOT
from stan_filter import stan_mask, read_obs_names

MASTER = ROOT / "00_DATA_FROZEN" / "03_analysis_ready" / "adata_master_normalized_magic.h5ad"
TAB    = ROOT / "06_tables" / "FIGURE_KD_clean"; TAB.mkdir(parents=True, exist_ok=True)

ANCHOR    = "APLNR"
MESO      = ["aPS-LC", "pPS-LC", "eNasM", "lNasM"]
TIMEPOINT = "72h"
MIN_PCT   = 0.02      # gene must be detected in >= 2% of cells at that site
MIN_RHO   = 0.10      # |rho| required in BOTH repeats for up_both / down_both


def _read_col(it):
    if isinstance(it, h5py.Group):
        ks = set(it.keys())
        if {"categories", "codes"} <= ks:
            return pd.Categorical.from_codes(
                it["codes"][:], categories=np.asarray(_read_col(it["categories"])))
        if "values" in ks:
            v = it["values"][:]
            return (np.array([x.decode() if isinstance(x, bytes) else x for x in v], dtype=object)
                    if v.dtype == object or v.dtype.kind in "SU" else v)
    a = it[:]
    return a.astype(str) if a.dtype.kind == "S" else a


def _csr(g):
    return sp.csr_matrix((g["data"][:], g["indices"][:], g["indptr"][:]),
                         shape=tuple(int(s) for s in g.attrs["shape"]))


def corr_vs_anchor(dense, anchor_vec, covars=None):
    """Spearman of every column of `dense` against `anchor_vec`.

    Rank first, then Pearson -- that IS Spearman, and it lets the whole thing be one
    matrix product instead of 25k scipy calls. If `covars` is given, the ranks are
    residualised on it first, which turns this into a partial Spearman.
    """
    # rankdata(..., axis=0) is vectorised; apply_along_axis would make one
    # python-level call per gene and is ~100x slower on a 15k-gene matrix.
    R = rankdata(dense, axis=0).astype(np.float32)
    a = rankdata(anchor_vec).astype(np.float32)[:, None]
    if covars is not None:
        C = np.hstack([np.ones((len(a), 1), np.float32), covars.astype(np.float32)])
        # least-squares residuals; rcond=None silences the future-warning and uses
        # the stable default cutoff
        R = R - C @ np.linalg.lstsq(C, R, rcond=None)[0]
        a = a - C @ np.linalg.lstsq(C, a, rcond=None)[0]
    R = R - R.mean(0); a = a - a.mean()
    rs = np.sqrt((R ** 2).sum(0)); as_ = np.sqrt((a ** 2).sum())
    denom = rs * as_
    out = np.zeros(R.shape[1], np.float32)
    ok = denom > 0
    out[ok] = ((R * a).sum(0)[ok]) / denom[ok]
    return out


def main():
    t0 = time.time()
    with h5py.File(MASTER, "r") as h:
        v = h["var"]; ik = v.attrs.get("_index", "_index")
        ik = ik.decode() if isinstance(ik, bytes) else ik
        vit = v[ik]; vraw = vit["values"][:] if isinstance(vit, h5py.Group) else vit[:]
        var = np.array([(x.decode() if isinstance(x, bytes) else str(x)).split(" (")[0]
                        for x in vraw])
        counts = _csr(h["layers"]["counts"])
        og = h["obs"]
        obs = pd.DataFrame({c: np.asarray(_read_col(og[c])).astype(str) for c in
                            ["condition", "site", "timepoint_true",
                             "predicted_leiden_sub_joint"]})
        # the ONLY MAGIC in this file -- a per-cell imputed APLNR, not a matrix
        magic_anchor = (np.asarray(_read_col(og[f"{ANCHOR}_magic"]), dtype=float)
                        if f"{ANCHOR}_magic" in og else None)

    # STAN BARCODE FILTER -- see stan_filter.py. Applied so this script uses the same
    # cell set as FIGURE2_CLEAN and the KD figure panels (29,419 cells, not 45,435).
    _keep = stan_mask(read_obs_names(MASTER), ROOT)
    counts = counts[_keep]
    obs = obs.loc[_keep].reset_index(drop=True)
    if magic_anchor is not None:
        magic_anchor = magic_anchor[_keep]

    obs["cond"]  = obs["condition"].map({"Ctrl": "Ctrl", "DOX": "KD"})
    obs["state"] = obs["predicted_leiden_sub_joint"]
    print(f"loaded {counts.shape} in {time.time()-t0:.0f}s"
          f"{'' if magic_anchor is not None else f'   [no {ANCHOR}_magic obs column -- raw only]'}")

    gidx = {g: i for i, g in enumerate(var)}
    if ANCHOR not in gidx:
        raise SystemExit(f"{ANCHOR} not in var")
    lib   = np.asarray(counts.sum(1)).ravel().astype(float); lib[lib == 0] = 1
    ngene = np.asarray((counts > 0).sum(1)).ravel().astype(float)

    base = ((obs.timepoint_true == TIMEPOINT) & obs.state.isin(MESO)).to_numpy()

    # ---- KD arm: is the APLNR axis PRESERVED under perturbation? ----
    # This is a control, not a discovery. If the genes that track APLNR in controls
    # still track it in the KD, then the knockdown lowered the receptor without
    # dismantling the programme it sits in -- which is what makes "the programme is
    # reduced" a coherent claim rather than a broken-cells artefact.
    # Previously this comparison existed ONLY as a figure (cell 37 of maturation_kd),
    # was pooled across repeats, and wrote no table. Now it is per repeat and tabulated.
    kd_cols = {}
    for site in ["Chicago", "Michigan"]:
        m = base & (obs.site == site).to_numpy() & (obs.cond == "KD").to_numpy()
        if m.sum() < 50:
            print(f"  [KD {site}: only {int(m.sum())} cells -- skipped]")
            continue
        pct = np.asarray((counts[m] > 0).mean(0)).ravel()
        keep = pct >= MIN_PCT; keep[gidx[ANCHOR]] = True
        sub = counts[m][:, keep]
        ln = sub.multiply(1.0 / lib[m][:, None]).tocsr(); ln.data = np.log1p(ln.data * 1e4)
        dense = np.asarray(ln.todense(), dtype=np.float32)
        kept = var[keep]; ai = int(np.where(kept == ANCHOR)[0][0])
        print(f"  KD {site}: {int(m.sum())} cells, {ANCHOR} in "
              f"{100*pct[gidx[ANCHOR]]:.1f}%, {keep.sum():,} genes")
        kd_cols[f"rho_KD_{site}"] = pd.Series(corr_vs_anchor(dense, dense[:, ai].copy()),
                                              index=kept)
        del dense, sub, ln

    base = base & (obs.cond == "Ctrl").to_numpy()

    res, pcts = {}, {}
    for site in ["Chicago", "Michigan"]:
        m = base & (obs.site == site).to_numpy()
        n = int(m.sum())
        pct = np.asarray((counts[m] > 0).mean(0)).ravel()
        keep = pct >= MIN_PCT
        keep[gidx[ANCHOR]] = True
        print(f"  {site}: {n} control mesoderm cells, {ANCHOR} detected in "
              f"{100*pct[gidx[ANCHOR]]:.1f}%, {keep.sum():,} genes pass {MIN_PCT:.0%}")

        sub = counts[m][:, keep]
        ln = sub.multiply(1.0 / lib[m][:, None]).tocsr(); ln.data = np.log1p(ln.data * 1e4)
        dense = np.asarray(ln.todense(), dtype=np.float32)
        kept_genes = var[keep]
        ai = int(np.where(kept_genes == ANCHOR)[0][0])
        anchor_raw = dense[:, ai].copy()

        cov = np.column_stack([np.log1p(lib[m]), np.log1p(ngene[m])]).astype(np.float32)
        r_raw = corr_vs_anchor(dense, anchor_raw)
        r_adj = corr_vs_anchor(dense, anchor_raw, covars=cov)

        d = pd.DataFrame({"gene": kept_genes,
                          f"rho_raw_{site}": r_raw,
                          f"rho_adj_{site}": r_adj}).set_index("gene")

        if magic_anchor is not None:
            # smoothed anchor vs UNsmoothed partners -- the asymmetry the notebook
            # actually introduces. Partners stay raw on purpose; that is the point.
            d[f"rho_mgc_{site}"] = corr_vs_anchor(dense, magic_anchor[m])
        del dense, sub, ln
        res[site] = d
        pcts[site] = pd.Series(pct[keep], index=kept_genes)

    out = res["Chicago"].join(res["Michigan"], how="inner")
    out["pct_Chicago"]  = pcts["Chicago"].reindex(out.index)
    out["pct_Michigan"] = pcts["Michigan"].reindex(out.index)

    c, mi = out["rho_raw_Chicago"], out["rho_raw_Michigan"]
    out["rho_raw_mean"] = (c + mi) / 2
    out["rho_raw_min"]  = np.minimum(c.abs(), mi.abs())
    out["up_both"]   = (c >=  MIN_RHO) & (mi >=  MIN_RHO)
    out["down_both"] = (c <= -MIN_RHO) & (mi <= -MIN_RHO)
    out["adj_agrees"] = (np.sign(out["rho_adj_Chicago"]) == np.sign(c)) & \
                        (np.sign(out["rho_adj_Michigan"]) == np.sign(mi))

    for k, v in kd_cols.items():
        out[k] = v.reindex(out.index)
    kdc = [k for k in kd_cols if k in out]
    if kdc:
        out["rho_KD_mean"] = out[kdc].mean(axis=1)
        out["ctrl_KD_sign_agree"] = np.sign(out["rho_raw_mean"]) == np.sign(out["rho_KD_mean"])

    # ---- backward-compatible aliases ----
    # maturation_kd.ipynb cells 34, 38 and 39 read `Chicago`, `Michigan` and `mean` from
    # the previous (producer-less) version of this file. Emitting them as aliases of the
    # raw columns keeps those cells working unchanged. Do not remove them without
    # grepping the notebook first.
    out["Chicago"] = out["rho_raw_Chicago"]
    out["Michigan"] = out["rho_raw_Michigan"]
    out["mean"] = out["rho_raw_mean"]

    out = out.sort_values("rho_raw_mean", ascending=False)
    p = TAB / "aplnr_percell_correlation_byrepeat.csv"
    out.reset_index().to_csv(p, index=False)

    # ---- diagnostics: the numbers that decide whether MAGIC mattered ----
    print(f"\nwrote {p.name}  ({len(out):,} genes)")
    print(f"  up_both   {int(out.up_both.sum()):,}"
          f"   down_both {int(out.down_both.sum()):,}"
          f"   (|rho| >= {MIN_RHO} same sign in both repeats, on RAW)")
    print(f"  of up_both, depth-adjusted keeps the sign in both repeats: "
          f"{100*out.loc[out.up_both,'adj_agrees'].mean():.0f}%")
    for s in ["Chicago", "Michigan"]:
        if f"rho_mgc_{s}" in out:
            rr, rm = out[f"rho_raw_{s}"], out[f"rho_mgc_{s}"]
            print(f"  {s}: mean |rho| raw {rr.abs().mean():.3f} vs magic-anchor "
                  f"{rm.abs().mean():.3f}  (inflation {rm.abs().mean()/max(rr.abs().mean(),1e-9):.2f}x)"
                  f"   spearman(raw, magic) = {pd.Series(rr).corr(pd.Series(rm), method='spearman'):.3f}")
    if "rho_KD_mean" in out:
        v = out.dropna(subset=["rho_KD_mean"])
        strong = v[v.up_both]
        print(f"\n  === APLNR axis in Ctrl vs KD (the preservation control) ===")
        print(f"  genes with both arms measurable: {len(v):,}")
        print(f"  pearson(rho_Ctrl, rho_KD)  = {v['rho_raw_mean'].corr(v['rho_KD_mean']):.3f}")
        print(f"  spearman                   = "
              f"{v['rho_raw_mean'].corr(v['rho_KD_mean'], method='spearman'):.3f}")
        print(f"  sign agreement, all genes  = {100*v['ctrl_KD_sign_agree'].mean():.0f}%")
        print(f"  sign agreement, up_both set= {100*strong['ctrl_KD_sign_agree'].mean():.0f}%"
              f"   (n={len(strong):,})")
        print(f"  => high values mean the KD lowered the receptor WITHOUT dismantling")
        print(f"     the programme it sits in. That is the point of this control.")
    print(f"\ntop 15 by raw mean rho (excluding {ANCHOR}):")
    show = out.drop(index=ANCHOR, errors="ignore").head(15)
    cols = [c for c in ["rho_raw_Chicago", "rho_raw_Michigan", "rho_adj_Michigan",
                        "rho_mgc_Michigan", "pct_Michigan"] if c in show]
    print(show[cols].round(3).to_string())
    print(f"\ndone in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
