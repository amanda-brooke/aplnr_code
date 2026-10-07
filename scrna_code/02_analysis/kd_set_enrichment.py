"""
kd_set_enrichment.py -- ARMS 3 and 3.5.

THE PROBLEM THIS SOLVES
The genome-wide gene list from the pseudobulk cannot be used. Permutation-tested, the
`reproducible` flag carries a ~29% false-discovery rate at the threshold where the
biology survives, and no threshold is both clean and useful (see kd_pseudobulk.py).
So "here are the top genes that went down" is not a claim this dataset can support.

What it CAN support is set-level tests. Asking "does this whole set of genes shift
relative to matched background?" is one hypothesis, not 17,000, so it does not pay the
multiple-testing cost that killed the gene list -- and it is robust to any individual
gene being noisy, which at n = 2 they all are.

  ARM 3   : does the APLNR-CORRELATED set fall in the KD?
            (co-expression -> dependence; the question the correlation analysis raises)
  ARM 3.5 : which BIOLOGICAL CATEGORIES fall in the KD, correlated with APLNR or not?
            (the honest replacement for "what else went down")

METHOD
For each set, compare its mean pseudobulk log2FC against N_PERM random gene sets drawn
to MATCH the set's baseMean decile distribution. Expression matching is not optional:
highly expressed genes have systematically different log2FC variance, and the KD's
compositional shift (KD cells are 25-38% less complex) hits abundance classes unequally.
An unmatched null would manufacture significance for any set of well-expressed genes.

Reported per set: observed mean log2FC, null mean and sd, z, empirical p, and BH across
sets within a tissue. Also n_down_0.3 vs its own null, which is less sensitive to a
single large outlier dragging the mean.

*** THE CATEGORY LISTS BELOW ARE PRE-SPECIFIED AND MUST STAY THAT WAY. ***
They are textbook membership, written down before looking at this dataset's log2FCs.
Do not add a gene because it looked good, and do not drop one because it did not --
either move makes the null invalid and turns a real test back into a fishing trip. If
a category genuinely needs revising, revise it, re-run everything, and say so.

The APLNR-correlated sets ARE data-derived, which is fine and is the whole point of
arm 3: they were defined by correlation in CONTROL cells only, and are then tested
against a KD-vs-Ctrl contrast. Control-derived, KD-tested -- no circularity.

Run: python kd_set_enrichment.py   (needs only KD_pseudobulk_results.csv; ~20 s)
"""
import numpy as np, pandas as pd
from pathlib import Path

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
from config import ROOT
TAB  = ROOT / "06_tables" / "FIGURE_KD_clean"

N_PERM     = 5000
N_BINS     = 20
DOWN_THR   = -0.3
TISSUES    = ["NasM", "all_mesoderm", "aPS-LC", "pPS-LC", "pluri_AE"]
PRIMARY    = "NasM"          # where APLNR is expressed in 75-83% of cells
SEED       = 0

# ----------------------------------------------------------- pre-specified categories
CATEGORIES = {
    "somite / fate specification": [
        "MESP1", "MESP2", "TBX6", "MSGN1", "DLL1", "DLL3", "LFNG", "HES7",
        "RIPPLY1", "RIPPLY2", "TCF15", "MEOX1", "ALDH1A2", "FOXC1", "FOXC2", "PAX3"],
    "EMT": [
        "SNAI1", "SNAI2", "TWIST1", "ZEB1", "ZEB2", "CDH2", "VIM", "FN1",
        "SERPINE2", "MMP2", "MMP14", "SPARC", "TGFB1", "LOXL2"],
    "migration / actomyosin": [
        "RAC1", "CDC42", "RHOA", "MYH9", "MYH10", "MYL9", "PFN1", "CFL1",
        "ARPC2", "ARPC3", "ACTR2", "ACTR3", "WASF1", "WASF2", "LIMK1", "ROCK1"],
    "filopodia": [
        "FSCN1", "MYO10", "ENAH", "VASP", "BAIAP2", "FMNL2", "FMNL3", "DAAM1",
        "EPS8", "RADIL", "CDC42", "IRSP53", "PLXNA2", "NRP1", "NRP2"],
    "WNT pathway": [
        "WNT3", "WNT3A", "WNT5A", "WNT5B", "WNT8A", "WNT11", "AXIN2", "SP5",
        "LEF1", "TCF7", "NKD1", "DKK1", "RSPO3", "FZD7", "LRP5", "LRP6"],
    "NODAL / TGFb": [
        "NODAL", "LEFTY1", "LEFTY2", "GDF3", "TDGF1", "SMAD2", "SMAD3",
        "FOXH1", "ACVR1B", "ACVR2B", "TGFBR1", "CER1"],
    "BMP": [
        "BMP2", "BMP4", "BMP7", "BMPR1A", "BMPR2", "ID1", "ID2", "ID3",
        "MSX1", "MSX2", "SMAD1", "SMAD5", "BAMBI", "BMPER"],
    "pluripotency (negative control)": [
        "POU5F1", "NANOG", "SOX2", "LIN28A", "DPPA4", "TDGF1", "ZFP42",
        "SALL4", "UTF1", "L1TD1"],
    "ribosome (artefact control)": [
        "RPS3", "RPS6", "RPL7", "RPL10", "RPL13A", "RPS19", "RPL23A", "RPS27A",
        "RPL5", "RPS4X", "RPL11", "RPS8", "RPL3", "RPS2"],
}


def matched_null(df, genes, n_perm, rng):
    """Mean log2FC and n-below-threshold for `genes`, vs baseMean-matched random sets.

    Vectorised over permutations: for each expression bin we draw all n_perm samples in
    one shot via argpartition on a random matrix (a without-replacement draw per row),
    then index the value array positionally. The obvious loop -- rng.choice inside a
    per-permutation loop with pandas .loc -- is ~500x slower and does the same thing.
    """
    sub = df.loc[df.index.intersection(genes)]
    if len(sub) < 4:
        return None
    vals = df["log2FC_shrunk"].to_numpy(float)
    bins = df["bin"].to_numpy()
    tgt = sub["bin"].value_counts()
    pool = {b: np.where(bins == b)[0] for b in tgt.index}
    if any(len(pool[b]) < k for b, k in tgt.items()):
        return None

    tot = np.zeros(n_perm); downs = np.zeros(n_perm); n_tot = 0
    for b, k in tgt.items():
        p = pool[b]
        if k == len(p):                      # whole bin taken; no randomness available
            pick = np.tile(p, (n_perm, 1))
        else:
            r = rng.random((n_perm, len(p)))
            pick = p[np.argpartition(r, k, axis=1)[:, :k]]
        v = vals[pick]                       # (n_perm, k)
        tot += v.sum(1); downs += (v < DOWN_THR).sum(1); n_tot += k
    means = tot / n_tot
    obs_m = float(sub["log2FC_shrunk"].mean())
    obs_d = int((sub["log2FC_shrunk"] < DOWN_THR).sum())
    return {
        "n_genes": len(sub), "n_missing": len(set(genes)) - len(sub),
        "mean_log2FC": obs_m,
        "null_mean": float(means.mean()), "null_sd": float(means.std()),
        "z": float((obs_m - means.mean()) / means.std()) if means.std() > 0 else np.nan,
        # two-sided, +1 correction so p is never exactly 0 with finite permutations
        "p_emp": float((np.abs(means - means.mean()) >= abs(obs_m - means.mean())).sum() + 1)
                 / (n_perm + 1),
        "n_down": obs_d, "null_n_down": float(downs.mean()),
        "z_ndown": float((obs_d - downs.mean()) / downs.std()) if downs.std() > 0 else np.nan,
        "genes_down": ", ".join(sub.index[sub["log2FC_shrunk"] < DOWN_THR]
                                .sort_values().tolist()),
        # full null draws, for plotting the distribution. main() strips these before
        # writing the CSV -- they are arrays, not scalars.
        "_null_means": means,
    }


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); q[o] = p[o] * n / (np.arange(n) + 1)
    q[o] = np.minimum.accumulate(q[o][::-1])[::-1]
    return np.clip(q, 0, 1)


def main():
    rng = np.random.default_rng(SEED)
    res = pd.read_csv(TAB / "KD_pseudobulk_results.csv")
    gcol = "gene" if "gene" in res.columns else res.columns[0]

    # ---- data-derived APLNR sets (defined in CONTROL cells only) ----
    sets = dict(CATEGORIES)
    p40 = TAB / "genes_rising_with_APLNR.csv"
    if p40.exists():
        g = pd.read_csv(p40)["gene"].tolist()
        sets["APLNR-correlated (top 40)"] = [x for x in g if x != "APLNR"]
    pbr = TAB / "aplnr_percell_correlation_byrepeat.csv"
    if pbr.exists():
        c = pd.read_csv(pbr)
        rho = "rho_raw_mean" if "rho_raw_mean" in c.columns else "mean"
        if "up_both" in c.columns:
            up = c[c.up_both & ~c.gene.astype(str).str.startswith("ENSG")]
            up = up.sort_values(rho, ascending=False).gene.tolist()
            sets["APLNR-correlated (up_both, top 200)"] = [x for x in up if x != "APLNR"][:200]

    rows = []
    for tissue in TISSUES:
        d = res[res.tissue == tissue].set_index(gcol)
        d = d[["log2FC_shrunk", "baseMean"]].dropna()
        if d.empty:
            continue
        d = d[~d.index.duplicated()]
        d["bin"] = pd.qcut(np.log10(d["baseMean"] + 1), N_BINS, labels=False,
                           duplicates="drop")
        for name, genes in sets.items():
            r = matched_null(d, genes, N_PERM, rng)
            if r:
                r.pop("_null_means", None)      # arrays don't belong in a CSV
                rows.append({"tissue": tissue, "set": name, **r})

    out = pd.DataFrame(rows)
    out["padj"] = np.nan
    for t, g in out.groupby("tissue"):
        out.loc[g.index, "padj"] = bh(g["p_emp"].to_numpy())
    out = out.sort_values(["tissue", "p_emp"])
    out.to_csv(TAB / "KD_set_enrichment.csv", index=False)

    show = ["set", "n_genes", "mean_log2FC", "null_mean", "z", "p_emp", "padj",
            "n_down", "null_n_down"]
    print(f"=== SET ENRICHMENT vs baseMean-matched background "
          f"({N_PERM:,} permutations) ===\n")
    print(f"--- {PRIMARY} (72h) --- APLNR expressed in 75-83% of these cells\n")
    print(out[out.tissue == PRIMARY][show].round(3).to_string(index=False))
    print("\n\n--- pluri_AE (72h) --- specificity control, like-for-like design\n")
    pa = out[out.tissue == "pluri_AE"]
    print(pa[show].round(3).to_string(index=False) if len(pa) else "  (not fitted)")
    print(f"\n\nwrote {(TAB/'KD_set_enrichment.csv').name}   "
          f"({len(out)} set x tissue tests)")
    print("\nRead the two control rows first. 'ribosome' should be UP or flat -- if it is\n"
          "significantly down, the compositional correction has failed. 'pluripotency'\n"
          "should be flat in mesoderm. If either misbehaves, no other row is trustworthy.")


if __name__ == "__main__":
    main()
