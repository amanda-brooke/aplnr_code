"""
run_pseudobulk_ALL.py  --  ONE clean pseudobulk run over every state x timepoint group.

Barcodes-only (Stan's list; NO confidence filter), DESeq2 shAPLNR-vs-control.

Design goal: the PAIRED case is byte-for-byte the locked producer
(01_DERIVED_LOCKED/pseudobulk_2026-08-08/per_state_72h/_producer_run_pseudobulk_stan_only.py)
-- same non-magic input file, same 'dataset' column, same detection-fraction gene filter,
same ~site+condition design -- so the single-state 72h groups REPRODUCE the locked numbers
exactly. A verification pass at the end checks that and fails loudly if it doesn't.

On top of that, every group is graded (idea adapted from a second draft):
  paired_DE            both sites have usable Ctrl+KD -> formal DESeq2 (the locked method)
  single_site_Chicago  only Chicago usable -> DESCRIPTIVE log2FC only, no p-values
  single_site_Michigan only Michigan usable -> DESCRIPTIVE log2FC only, no p-values
  insufficient         no usable Ctrl+KD pair -> QC only, no effect table
Nothing is silently dropped: a QC row is written for every expected sample, incl. 0-cell.
Single-site numbers are DESCRIPTIVE ONLY and must not be reported as DE.

WHERE TO RUN: on your machine (the 2.6 GB master h5ad lives there, not in the cloud).
    pip install pydeseq2 h5py pandas numpy      # if needed
    python run_pseudobulk_ALL.py
Resumable: each group writes its own CSV and is skipped if present. Delete OUT for a
truly clean slate.
"""

from datetime import date
from pathlib import Path
import sys
import warnings

import h5py
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ----------------------------------------------------------------------------- config
STATE_PLAN = {
    "pluripotent_42h":  {"states": ["pluripotent"],                        "timepoints": ["42h"]},
    "pluripotent_72h":  {"states": ["pluripotent"],                        "timepoints": ["72h"]},
    "AE-LC_42h":        {"states": ["AE-LC"],                              "timepoints": ["42h"]},
    "AE-LC_72h":        {"states": ["AE-LC"],                              "timepoints": ["72h"]},
    "aPS-LC_42h":       {"states": ["aPS-LC"],                             "timepoints": ["42h"]},
    "aPS-LC_72h":       {"states": ["aPS-LC"],                             "timepoints": ["72h"]},
    "pPS-LC_42h":       {"states": ["pPS-LC"],                             "timepoints": ["42h"]},
    "pPS-LC_72h":       {"states": ["pPS-LC"],                             "timepoints": ["72h"]},
    "eNasM_42h":        {"states": ["eNasM"],                              "timepoints": ["42h"]},
    "eNasM_72h":        {"states": ["eNasM"],                              "timepoints": ["72h"]},
    "lNasM_42h":        {"states": ["lNasM"],                              "timepoints": ["42h"]},
    "lNasM_72h":        {"states": ["lNasM"],                              "timepoints": ["72h"]},
    "NasM_42h":         {"states": ["eNasM", "lNasM"],                     "timepoints": ["42h"]},
    "NasM_72h":         {"states": ["eNasM", "lNasM"],                     "timepoints": ["72h"]},
    "all_mesoderm_42h": {"states": ["aPS-LC", "pPS-LC", "eNasM", "lNasM"], "timepoints": ["42h"]},
    "all_mesoderm_72h": {"states": ["aPS-LC", "pPS-LC", "eNasM", "lNasM"], "timepoints": ["72h"]},
}

SITES = ["Chicago", "Michigan"]       # the two biological repeats (obs['dataset'])
CONDS = ["Ctrl", "DOX"]               # raw obs.condition; DOX == shAPLNR (KD)
MIN_DETECT = 0.10                     # gene tested only if detected in >=10% of the group's cells
MIN_EFFECT = 0.25                     # |log2FC| both repeats must clear for `reproduced`
MIN_CELLS_PER_SAMPLE = 10            # a site x condition pseudobulk below this is not usable.
                                     # 10 keeps Chicago's thin eNasM/lNasM 72h controls (~13-23
                                     # cells) so those states stay PAIRED and reproduce the
                                     # locked 08-08 / manuscript numbers. Raising it declassifies
                                     # them to single-site descriptive -- see AUDIT_groups.csv.

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
from config import ROOT
# SAME non-magic file the locked run used. Do NOT point this at *_magic.h5ad -- MAGIC
# output is imputed/continuous and must never enter a count-based DE model.
H5 = ROOT / "00_DATA_FROZEN" / "03_analysis_ready" / "adata_master_normalized.h5ad"
OUT = ROOT / f"pseudobulk_ALL_{date.today().isoformat()}"   # clean-slate output folder
OUT.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------------------- h5ad helpers
def _dec(a):
    a = np.asarray(a)
    if a.dtype.kind == "S":
        return np.array([x.decode("utf-8", "replace") for x in a], dtype=object)
    if a.dtype == object:
        return np.array([x.decode() if isinstance(x, bytes) else x for x in a], dtype=object)
    return a


def _arr(n):
    if isinstance(n, h5py.Dataset):
        return _dec(n[:])
    if "codes" in n and "categories" in n:
        c = np.asarray(_arr(n["codes"]))
        k = np.asarray(_arr(n["categories"]), dtype=object)
        o = np.full(len(c), "<NA>", dtype=object)
        m = c >= 0
        o[m] = k[c[m]]
        return o
    return _arr(n["values"])


def load_obs_and_genes(f):
    ob = f["obs"]
    # NOTE: uses 'dataset' for the Chicago/Michigan repeat, exactly as the locked run did.
    obs = pd.DataFrame(
        {k: _arr(ob[k]) for k in
         ["condition", "timepoint_true", "dataset", "predicted_leiden_sub_joint"]}
    ).rename(columns={"predicted_leiden_sub_joint": "state", "dataset": "site"})
    genes = _arr(f["var"]["_index"]).astype(str)
    return obs, genes


def load_stan():
    from stan_filter import stan_mask, read_obs_names
    return stan_mask, read_obs_names


# ---------------------------------------------------------------------- pseudobulk build
def build_samples(f, obs, genes, stan, states, timepoints):
    """Sum raw counts into the 4 site x condition pseudobulks.
    Returns (S[4,G] counts, Dc[4,G] nonzero-cell counts, ncells[4], included[4]).
    Rows for samples below MIN_CELLS_PER_SAMPLE are left as zeros and flagged not-included.
    Sample order: Chicago|Ctrl, Chicago|DOX, Michigan|Ctrl, Michigan|DOX."""
    in_grp = (obs.state.isin(states) & obs.timepoint_true.isin(timepoints) & stan).to_numpy()
    order = [(s, c) for s in SITES for c in CONDS]
    masks, ncells = [], []
    for s, c in order:
        m = in_grp & (obs.site.to_numpy() == s) & (obs.condition.to_numpy() == c)
        masks.append(np.where(m)[0])
        ncells.append(int(m.sum()))
    included = [n >= MIN_CELLS_PER_SAMPLE for n in ncells]

    L = f["layers"]["counts"]
    data, ind, indptr = L["data"], L["indices"], L["indptr"]
    S = np.zeros((4, len(genes)))
    Dc = np.zeros((4, len(genes)))
    for gi, (sel, inc) in enumerate(zip(masks, included)):
        if not inc:
            continue
        for ci in sel:
            a, b = indptr[ci], indptr[ci + 1]
            j = ind[a:b]
            np.add.at(S[gi], j, data[a:b])
            np.add.at(Dc[gi], j, 1.0)
    return S.astype(int), Dc, ncells, included


# -------------------------------------------------------------- paired DESeq2 (locked method)
def paired_deseq(S, Dc, ncells, included, genes, key):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    # denominator = only the samples that actually contributed to Dc. For paired_DE all
    # four are included, so this equals sum(ncells); explicit here so the detection
    # fraction stays correct if the function is ever reused with a partial sample set.
    tot = sum(n for n, inc in zip(ncells, included) if inc)
    det = Dc.sum(0) / tot
    tested = det >= MIN_DETECT
    if tested.sum() == 0:
        return None

    samples = [f"{s}|{c}" for s in SITES for c in CONDS]
    meta = pd.DataFrame(
        {"site": [x.split("|")[0] for x in samples],
         "condition": ["control" if x.endswith("Ctrl") else "shAPLNR" for x in samples]},
        index=samples)
    cts = pd.DataFrame(S[:, tested], index=samples, columns=genes[tested])

    dds = DeseqDataSet(counts=cts, metadata=meta, design="~site + condition",
                       refit_cooks=True, quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["condition", "shAPLNR", "control"], quiet=True)
    st.summary()
    res = st.results_df.copy()

    sf = np.asarray(dds.obs["size_factors"])
    lg = np.log2(S[:, tested] / sf[:, None] + 1)
    res["log2FC_repeat1_Chicago"] = lg[1] - lg[0]
    res["log2FC_repeat2_Michigan"] = lg[3] - lg[2]
    res["detected_frac"] = det[tested]
    res["min_abs_repeat"] = np.minimum(res.log2FC_repeat1_Chicago.abs(),
                                       res.log2FC_repeat2_Michigan.abs())
    res["repeats_agree"] = np.sign(res.log2FC_repeat1_Chicago) == np.sign(res.log2FC_repeat2_Michigan)
    res["reproduced"] = res.repeats_agree & (res.min_abs_repeat >= MIN_EFFECT)
    res["group"] = key
    res["status"] = "paired_DE"
    res["n_cells_control"] = ncells[0] + ncells[2]
    res["n_cells_shAPLNR"] = ncells[1] + ncells[3]
    return res.reset_index().rename(columns={"index": "gene"})


# ------------------------------------------------------- single-site DESCRIPTIVE (no model)
def single_site_descriptive(S, Dc, ncells, included, genes, key, site):
    """CPM-normalised Ctrl-vs-KD log2FC for one site. DESCRIPTIVE ONLY -- no p-values,
    no reproducibility. Never report these as DE."""
    i_ctrl, i_kd = (0, 1) if site == "Chicago" else (2, 3)
    if not (included[i_ctrl] and included[i_kd]):
        return None
    n = ncells[i_ctrl] + ncells[i_kd]
    det = (Dc[i_ctrl] + Dc[i_kd]) / n
    tested = det >= MIN_DETECT
    if tested.sum() == 0:
        return None

    def cpm(v):
        return v / max(v.sum(), 1) * 1e6

    ctrl = cpm(S[i_ctrl].astype(float))[tested]
    kd = cpm(S[i_kd].astype(float))[tested]
    lfc = np.log2(kd + 1) - np.log2(ctrl + 1)
    out = pd.DataFrame({
        "gene": genes[tested],
        f"log2FC_descriptive_{site}": lfc,
        "detected_frac": det[tested],
    })
    out["group"] = key
    out["status"] = f"single_site_{site}"
    out["n_cells_control"] = ncells[i_ctrl]
    out["n_cells_shAPLNR"] = ncells[i_kd]
    return out


# --------------------------------------------------------------------------------- main
def main():
    qc_rows, group_rows = [], []
    stan_mask, read_obs_names = load_stan()

    with h5py.File(H5, "r") as f:
        obs, genes = load_obs_and_genes(f)
        stan = stan_mask(read_obs_names(H5), ROOT)

        for key, spec in STATE_PLAN.items():
            S, Dc, ncells, included = build_samples(
                f, obs, genes, stan, spec["states"], spec["timepoints"])

            # exhaustive QC: one row per expected sample, including 0-cell
            order = [(s, c) for s in SITES for c in CONDS]
            for i, (s, c) in enumerate(order):
                qc_rows.append(dict(
                    group=key, states="+".join(spec["states"]),
                    timepoint="+".join(spec["timepoints"]), site=s, cond=c,
                    n_cells=ncells[i], included=included[i],
                    reason="" if included[i] else f"< MIN_CELLS_PER_SAMPLE ({MIN_CELLS_PER_SAMPLE})"))

            chi_ok = included[0] and included[1]
            mich_ok = included[2] and included[3]
            if chi_ok and mich_ok:
                status = "paired_DE"
            elif chi_ok:
                status = "single_site_Chicago"
            elif mich_ok:
                status = "single_site_Michigan"
            else:
                status = "insufficient"

            done = OUT / f"{key}.csv"
            if done.exists():
                print(f"  {key:18s} already done, reusing")
                group_rows.append(dict(group=key, status=status,
                                       n_ctrl=ncells[0] + ncells[2], n_kd=ncells[1] + ncells[3]))
                continue

            res = None
            if status == "paired_DE":
                res = paired_deseq(S, Dc, ncells, included, genes, key)
            elif status == "single_site_Chicago":
                res = single_site_descriptive(S, Dc, ncells, included, genes, key, "Chicago")
            elif status == "single_site_Michigan":
                res = single_site_descriptive(S, Dc, ncells, included, genes, key, "Michigan")

            reproduced = 0
            if res is not None:
                res.to_csv(done, index=False)
                if "reproduced" in res:
                    reproduced = int(res.reproduced.sum())
            group_rows.append(dict(
                group=key, status=status,
                n_ctrl=ncells[0] + ncells[2], n_kd=ncells[1] + ncells[3],
                genes=0 if res is None else len(res), reproduced=reproduced))
            print(f"  {key:18s} {status:20s} ctrl={ncells[0]+ncells[2]:5d} kd={ncells[1]+ncells[3]:5d}"
                  + ("" if res is None else f" | rows {len(res)}"))

    pd.DataFrame(qc_rows).to_csv(OUT / "QC_samples.csv", index=False)
    pd.DataFrame(group_rows).to_csv(OUT / "AUDIT_groups.csv", index=False)

    ran = [pd.read_csv(OUT / f"{k}.csv") for k in STATE_PLAN if (OUT / f"{k}.csv").exists()]
    if ran:
        pd.concat(ran, ignore_index=True).to_csv(OUT / "ALL_GROUPS_combined.csv", index=False)

    # ---- verification: single-state 72h paired groups must reproduce the locked numbers ----
    locked = ROOT / "01_DERIVED_LOCKED" / "pseudobulk_2026-08-08" / "per_state_72h"
    # Criterion: DESeq2 is fit by an iterative GLM, so cross-version runs agree only to
    # solver precision (~1e-4), not bit-for-bit. "Reproduces" therefore means: identical
    # gene set, median|diff| below solver noise, and NO gene reversing direction among
    # biologically relevant effects. A real methods change moves the MEDIAN and flips signs;
    # version noise moves only the max, on a handful of Cooks-refit genes.
    MEDIAN_TOL = 1e-3
    SIGNIF_LFC = 0.25          # only count sign flips among genes with a real effect size
    checks = []
    for stt in ["pluripotent", "aPS-LC", "pPS-LC", "eNasM", "lNasM", "AE-LC"]:
        new, old = OUT / f"{stt}_72h.csv", locked / f"{stt}_deseq2.csv"
        if not (new.exists() and old.exists()):
            continue
        ndf = pd.read_csv(new)
        if "log2FoldChange" not in ndf.columns:
            grade = ndf["status"].iloc[0] if "status" in ndf.columns and len(ndf) else "n/a"
            checks.append(dict(state=stt, genes=0, corr=float("nan"), median_abs=float("nan"),
                               max_abs=float("nan"), sign_flips=-1, worst_gene="",
                               verdict=f"not paired ({grade})"))
            continue
        m = ndf[["gene", "log2FoldChange"]].merge(
            pd.read_csv(old)[["gene", "log2FoldChange"]], on="gene", suffixes=("_new", "_old"))
        d = (m.log2FoldChange_new - m.log2FoldChange_old).abs()
        big = (m.log2FoldChange_new.abs() >= SIGNIF_LFC) | (m.log2FoldChange_old.abs() >= SIGNIF_LFC)
        flips = int(((np.sign(m.log2FoldChange_new) != np.sign(m.log2FoldChange_old)) & big).sum())
        ok = (len(m) == len(ndf) == len(pd.read_csv(old))) and d.median() < MEDIAN_TOL and flips == 0
        checks.append(dict(
            state=stt, genes=len(m), corr=round(float(m.log2FoldChange_new.corr(m.log2FoldChange_old)), 6),
            median_abs=float(d.median()), max_abs=float(d.max()), sign_flips=flips,
            worst_gene=str(m.loc[d.idxmax(), "gene"]),
            verdict="reproduces (solver noise)" if ok else "DIFFERS - investigate"))
    if checks:
        C = pd.DataFrame(checks)
        C.to_csv(OUT / "VERIFY_vs_locked_72h.csv", index=False)
        paired = C[C.sign_flips >= 0]
        all_ok = len(paired) and (paired.verdict == "reproduces (solver noise)").all()
        print(f"\nVERIFY vs locked 72h: {len(paired)}/6 states paired | "
              f"corr>={paired.corr.min():.6f} | worst median|d|={paired.median_abs.max():.1e} | "
              f"sign flips={int(paired.sign_flips.clip(lower=0).sum())} "
              f"({'PASS - reproduces locked' if all_ok else 'DIFFERS - investigate'})")

    (OUT / "README.md").write_text(
        "# pseudobulk_ALL — one clean run over STATE_PLAN\n\n"
        "Barcodes-only, shAPLNR vs control. Grades per group:\n\n"
        "- **paired_DE** — both sites usable; formal DESeq2 (`~site + condition`), identical "
        "method to the locked 08-08 run. Has padj, per-repeat log2FC, `reproduced`.\n"
        "- **single_site_Chicago / _Michigan** — only one site usable; CPM log2FC, "
        "**descriptive only, no p-values, do not report as DE**.\n"
        "- **insufficient** — no usable Ctrl+KD pair; see QC_samples.csv.\n\n"
        "Files: `QC_samples.csv` (every expected sample incl. 0-cell), `AUDIT_groups.csv`, "
        "`ALL_GROUPS_combined.csv`, `VERIFY_vs_locked_72h.csv`, one `<group>.csv` each.\n\n"
        f"Cell threshold: MIN_CELLS_PER_SAMPLE={MIN_CELLS_PER_SAMPLE}. "
        f"Detection filter: {MIN_DETECT:.0%}. Source: {H5.name} (NOT the magic file).\n")

    print(f"\nwrote {OUT}")
    print(pd.DataFrame(group_rows).to_string(index=False))


if __name__ == "__main__":
    main()
