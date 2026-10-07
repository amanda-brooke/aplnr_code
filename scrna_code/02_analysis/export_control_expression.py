"""
export_control_expression.py  --  CONTROL-only pseudobulk EXPRESSION across tissue x time.

Companion to run_pseudobulk_ALL.py. That script gives shAPLNR-vs-control log2 fold-changes
(a contrast). This one gives expression LEVELS in control cells only, so you can show normal
developmental progression: each gene across states (tissues) and timepoints, no KD involved.

Same source object and same Stan barcode filter as the DE pipeline, so cell definitions match.
Values are CPM from summed raw counts (library-size normalised), reported as log2(CPM+1).
Pseudobulks pool the two biological repeats (Chicago+Michigan); per-repeat matrices are also
written so you can check reproducibility. Nothing here is a statistical test.

Run on your machine:  python export_control_expression.py
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
STATES = ["pluripotent", "aPS-LC", "pPS-LC", "eNasM", "lNasM", "AE-LC"]
TPS = ["42h", "72h"]
SITES = ["Chicago", "Michigan"]
COND = "Ctrl"               # control only
MIN_DETECT = 0.10           # keep a gene if detected in >=10% of cells in at least one group
MIN_CELLS = 10              # a state x tp control pseudobulk below this is written but FLAGGED

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
from config import ROOT
H5 = ROOT / "00_DATA_FROZEN" / "03_analysis_ready" / "adata_master_normalized.h5ad"
OUT = ROOT / f"control_expression_{date.today().isoformat()}"
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


def main():
    from stan_filter import stan_mask, read_obs_names

    with h5py.File(H5, "r") as f:
        ob = f["obs"]
        obs = pd.DataFrame(
            {k: _arr(ob[k]) for k in
             ["condition", "timepoint_true", "dataset", "predicted_leiden_sub_joint"]}
        ).rename(columns={"predicted_leiden_sub_joint": "state", "dataset": "site"})
        genes = _arr(f["var"]["_index"]).astype(str)
        stan = stan_mask(read_obs_names(H5), ROOT)

        L = f["layers"]["counts"]
        data, ind, indptr = L["data"], L["indices"], L["indptr"]

        # columns in state-major, time order so progression reads left->right
        cols = [f"{st}_{tp}" for st in STATES for tp in TPS]
        cpm_pool = pd.DataFrame(index=genes, dtype=float)
        det_pool = pd.DataFrame(index=genes, dtype=float)
        cpm_site = {s: pd.DataFrame(index=genes, dtype=float) for s in SITES}
        qc = []

        base = (obs.condition.to_numpy() == COND) & stan
        for st in STATES:
            for tp in TPS:
                col = f"{st}_{tp}"
                grp = base & (obs.state.to_numpy() == st) & (obs.timepoint_true.to_numpy() == tp)

                # pooled (both repeats)
                sel = np.where(grp)[0]
                S = np.zeros(len(genes)); D = np.zeros(len(genes))
                for ci in sel:
                    a, b = indptr[ci], indptr[ci + 1]
                    j = ind[a:b]
                    np.add.at(S, j, data[a:b]); np.add.at(D, j, 1.0)
                n = len(sel)
                if n > 0:
                    cpm_pool[col] = np.log2(S / max(S.sum(), 1) * 1e6 + 1)
                    det_pool[col] = D / n
                qc.append(dict(state=st, timepoint=tp, site="pooled", n_cells=n,
                               flag="" if n >= MIN_CELLS else f"< {MIN_CELLS} cells"))

                # per repeat
                for s in SITES:
                    ssel = np.where(grp & (obs.site.to_numpy() == s))[0]
                    Ss = np.zeros(len(genes))
                    for ci in ssel:
                        a, b = indptr[ci], indptr[ci + 1]
                        np.add.at(Ss, ind[a:b], data[a:b])
                    if len(ssel) > 0:
                        cpm_site[s][col] = np.log2(Ss / max(Ss.sum(), 1) * 1e6 + 1)
                    qc.append(dict(state=st, timepoint=tp, site=s, n_cells=len(ssel),
                                   flag="" if len(ssel) >= MIN_CELLS else f"< {MIN_CELLS} cells"))

    # keep genes detected in >=1 group; order columns
    keep = (det_pool.reindex(columns=[c for c in cols if c in det_pool]) >= MIN_DETECT).any(axis=1)
    order = [c for c in cols if c in cpm_pool.columns]
    cpm_pool = cpm_pool.loc[keep, order]
    det_pool = det_pool.loc[keep, [c for c in cols if c in det_pool.columns]]

    cpm_pool.to_csv(OUT / "control_expression_log2cpm_pooled.csv")
    det_pool.to_csv(OUT / "control_expression_detected_frac.csv")
    for s in SITES:
        m = cpm_site[s]
        m.loc[m.index.intersection(cpm_pool.index),
              [c for c in cols if c in m.columns]].to_csv(OUT / f"control_expression_log2cpm_{s}.csv")
    pd.DataFrame(qc).to_csv(OUT / "QC_control_cells.csv", index=False)

    (OUT / "README.md").write_text(
        "# control_expression — normal progression across tissue x time\n\n"
        "Control cells only, Stan-filtered, same source object as the DE pipeline.\n"
        "Values = log2(CPM+1) from summed raw counts. NOT a statistical test.\n\n"
        "- `control_expression_log2cpm_pooled.csv` — genes x `state_timepoint` (both repeats "
        "pooled). **Main matrix for heatmaps / progression lines.**\n"
        "- `control_expression_detected_frac.csv` — fraction of cells expressing each gene, same layout.\n"
        "- `control_expression_log2cpm_{Chicago,Michigan}.csv` — per repeat, to check reproducibility.\n"
        "- `QC_control_cells.csv` — control cell counts per group; rows flagged when < "
        f"{MIN_CELLS} cells (eNasM/lNasM at 42h are near-empty — read those with care).\n\n"
        "Columns are ordered tissue-major then time (42h, 72h), so each tissue's progression "
        "reads left-to-right. Genes kept if detected in >=10% of cells in at least one group.\n")

    print(f"wrote {OUT}")
    print(f"  matrix: {cpm_pool.shape[0]} genes x {cpm_pool.shape[1]} groups")
    print(pd.DataFrame(qc).query("site=='pooled'")[["state", "timepoint", "n_cells", "flag"]].to_string(index=False))


if __name__ == "__main__":
    main()
