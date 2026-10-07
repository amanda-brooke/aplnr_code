"""
dataconfig.py — load a dataset with its analysis config (exclusions applied ONCE).

The exclusions / reference window / units live in analysis.yml, not in the notebooks.
Notebooks call load_dataset("shaplnr_mps") and get a filtered, coordinate-ready df.
For a quick exploration you can override any field: load_dataset(..., overrides={...}).
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import yaml

_HERE = Path(__file__).resolve().parent
DEFAULT_CFG = _HERE / "analysis.yml"


def load_configs(cfg_path=DEFAULT_CFG) -> dict:
    with open(cfg_path) as f:
        return yaml.safe_load(f)


def get_config(name, cfg_path=DEFAULT_CFG, overrides: dict | None = None) -> dict:
    cfg = dict(load_configs(cfg_path)[name])
    if overrides:
        cfg.update(overrides)
    return cfg


def load_dataset(name, data_dir="data", cfg_path=DEFAULT_CFG, overrides: dict | None = None):
    """Return (df, cfg): df filtered per config, with radius / dist_from_edge / dnorm columns."""
    cfg = get_config(name, cfg_path, overrides)
    df = pd.read_pickle(Path(data_dir) / cfg["radial_file"])

    # --- exclusions (the one place they live) ---
    for off in (cfg.get("exclude_offsets") or []):
        df = df[~df["File_ID"].str.contains(off)]
    sf = cfg.get("size_filter")
    if sf:
        for col, vals in sf.items():
            df = df[df[col].isin(vals)]

    df = df.copy()
    df["radius"] = df["diameter_um"] / 2.0
    df["dist_from_edge"] = df["radius"] - df["radial_distance"]          # microns
    df["dnorm"] = 1.0 - df["radial_distance"] / df["radius"]             # fraction of radius
    return df, cfg


def coord_and_edges(cfg, n_um=30, n_norm=25):
    """Pick the coordinate column + bin edges implied by cfg['units']."""
    if cfg["units"] == "um":
        return "dist_from_edge", np.linspace(0, 300, n_um + 1)
    elif cfg["units"] == "fraction_radius":
        return "dnorm", np.linspace(0, 1.0, n_norm + 1)
    raise ValueError(f"unknown units: {cfg['units']}")


# --------------------------------------------------------------------------- master
def _read_any(path):
    """Read a table whether it's .parquet or .pkl."""
    path = Path(path)
    return pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_pickle(path)


def _find(path_noext, required=False):
    """Return the existing .parquet or .pkl for a base path. Raise a clear error if required & missing."""
    for ext in (".parquet", ".pkl"):
        p = Path(str(path_noext) + ext)
        if p.exists():
            return p
    if required:
        import os
        base = Path(path_noext).parent
        have = sorted(f.name for f in base.glob("*")) if base.exists() else "(directory does not exist)"
        raise FileNotFoundError(
            f"Master file '{Path(path_noext).name}.(parquet|pkl)' not found.\n"
            f"  looked in : {base.resolve()}\n"
            f"  cwd       : {os.getcwd()}\n"
            f"  present   : {have}\n"
            f"  -> run `python pipeline_v2/build_master.py` (or check you're in the repo root).")
    return None


def load_master(experiment=None, processed_dir="data/processed",
                include_excluded=False, **filters):
    """
    Load the master measurements table (built by build_master.py).

    experiment       : one experiment ("parent"/"sox2"/"shaplnr") or list, or None for all.
    include_excluded : default False -> applies the GLOBAL exclusion set (df[~df.excluded]),
                       so every analysis inherits the same exclusions automatically.
    **filters        : simple equality filters on any column, e.g.
                       load_master("parent", cell_type="mesoderm", geometry="circle").

    Returns one tidy DataFrame. Pair with load_tracks()/load_boundaries() for the
    per-track summaries and the mu+k*sigma boundary values.
    """
    mdir = Path(processed_dir) / "master"
    exps = [experiment] if isinstance(experiment, str) else experiment
    if exps:
        files = [_find(mdir / e) for e in exps]
    else:
        files = sorted(mdir.glob("*.parquet")) + sorted(mdir.glob("*.pkl"))
    files = [f for f in files if f is not None]
    if not files:
        raise FileNotFoundError(f"no master table in {mdir} — run build_master.py first")
    df = pd.concat((_read_any(f) for f in files), ignore_index=True)
    if not include_excluded and "excluded" in df.columns:
        df = df[~df["excluded"]].copy()
    for col, val in filters.items():
        df = df[df[col].isin(val) if isinstance(val, (list, tuple, set)) else df[col] == val]
    return df.reset_index(drop=True)


def load_tracks(processed_dir="data/processed", include_excluded=False, **filters):
    """Per-track summary table (one row per cell): duration, straightness, max_dist_edge, ever_beyond_1_5..."""
    df = _read_any(_find(Path(processed_dir) / "master_tracks", required=True))
    if not include_excluded and "excluded" in df.columns:
        df = df[~df["excluded"]].copy()
    for col, val in filters.items():
        df = df[df[col].isin(val) if isinstance(val, (list, tuple, set)) else df[col] == val]
    return df.reset_index(drop=True)


def load_boundaries(processed_dir="data/processed"):
    """The mu / sigma / mu+k*sigma boundary values per dataset & k (the numbers behind beyond_*)."""
    return _read_any(_find(Path(processed_dir) / "boundaries", required=True))
