"""
master_io.py -- the single entry point for reading the master object.

WHY THIS EXISTS
    00_DATA_FROZEN/03_analysis_ready/adata_master_normalized.h5ad holds all 45,435
    cells. Stan's barcode filter is opt-in, so forgetting it fails SILENTLY -- which is
    what produced the 2026-08-02 pseudobulk on the wrong cell set.

    Every function here applies the filter BY DEFAULT. Getting unfiltered cells requires
    passing allow_unfiltered=True, which is greppable and shows up in review.

USAGE
    from master_io import open_master, master_cells

    h5, keep = open_master()              # 29,419 Stan cells
    obs = master_cells()                  # filtered obs DataFrame, with .keep_idx

    h5, keep = open_master(allow_unfiltered=True)   # must be explicit

CELL SET POLICY (decided 2026-08-07)
    Stan's barcode list ONLY. No label-transfer confidence filter -- Stan does not use
    one, and chat_scores_maturation_bidirectional.py never applied one either. Applying
    it in the pseudobulk but not the score was an inconsistency; this is the fix.

    If a specific analysis needs a confidence cut, apply it locally and SAY SO in that
    script's docstring. Do not add it here.
"""
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from stan_filter import stan_mask, read_obs_names

from config import ROOT
MASTER = ROOT / "00_DATA_FROZEN" / "03_analysis_ready" / "adata_master_normalized.h5ad"
MASTER_MAGIC = ROOT / "00_DATA_FROZEN" / "03_analysis_ready" / "adata_master_normalized_magic.h5ad"

EXPECTED_TOTAL = 45435
EXPECTED_STAN = 29419


def _dec(a):
    a = np.asarray(a)
    if a.dtype.kind == "S":
        return np.array([x.decode("utf-8", "replace") for x in a], dtype=object)
    if a.dtype == object:
        return np.array([x.decode() if isinstance(x, bytes) else x for x in a], dtype=object)
    return a


def _arr(node):
    if isinstance(node, h5py.Dataset):
        return _dec(node[:])
    if "codes" in node and "categories" in node:
        c = np.asarray(_arr(node["codes"]))
        k = np.asarray(_arr(node["categories"]), dtype=object)
        out = np.full(len(c), "<NA>", dtype=object)
        m = c >= 0
        out[m] = k[c[m]]
        return out
    return _arr(node["values"])


def stan_keep_idx(path=MASTER, verbose=True):
    """Integer positions of the cells that pass Stan's barcode filter."""
    names = read_obs_names(path)
    if len(names) != EXPECTED_TOTAL:
        raise RuntimeError(f"{path.name} has {len(names):,} cells, expected {EXPECTED_TOTAL:,}")
    keep = np.flatnonzero(np.asarray(stan_mask(names, ROOT)))
    if len(keep) != EXPECTED_STAN:
        raise RuntimeError(f"Stan filter kept {len(keep):,}, expected {EXPECTED_STAN:,}")
    if verbose:
        print(f"[master_io] {len(keep):,}/{len(names):,} cells (Stan barcodes, no confidence cut)")
    return keep


def open_master(magic=False, allow_unfiltered=False, verbose=True):
    """Return (open h5py.File, keep_idx). Caller closes the file.

    allow_unfiltered=True returns every cell and is deliberately awkward to type.
    """
    path = MASTER_MAGIC if magic else MASTER
    f = h5py.File(path, "r")
    if allow_unfiltered:
        n = len(read_obs_names(path))
        if verbose:
            print(f"[master_io] WARNING: unfiltered, all {n:,} cells")
        return f, np.arange(n)
    return f, stan_keep_idx(path, verbose=verbose)


def master_cells(columns=None, magic=False, allow_unfiltered=False, verbose=True):
    """Filtered obs as a DataFrame. `.attrs['keep_idx']` holds the source positions."""
    f, keep = open_master(magic=magic, allow_unfiltered=allow_unfiltered, verbose=verbose)
    try:
        ob = f["obs"]
        cols = columns or [k for k in ob.keys() if k != "_index"]
        d = {}
        for k in cols:
            if k not in ob:
                continue
            v = np.asarray(_arr(ob[k]))
            d[k] = v[keep]
        idx = np.asarray(_arr(ob["_index"]))[keep]
        df = pd.DataFrame(d, index=idx)
    finally:
        f.close()
    df.attrs["keep_idx"] = keep
    return df


if __name__ == "__main__":
    df = master_cells(columns=["predicted_leiden_sub_joint", "timepoint_true",
                               "condition", "site"])
    print(df.shape)
    print(pd.crosstab(df.timepoint_true, df.condition))
