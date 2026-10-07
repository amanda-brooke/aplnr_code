"""
prep — the "load data once" layer.

Goal: every figure gets its data from the SAME two tables, built once here, instead of
each plot re-loading and re-filtering on its own. That removes the risk that two panels
quietly filter differently.

Two levels, two tables:
  - edges     : one row per STEP (per-frame movement). Edge-level metrics live here
                (SPEED, radial_velocity, cos_theta). Base track filter already applied.
  - per_track : one row per CELL (whole trajectory), derived from edges. Track-level
                metrics live here (mean speed, net displacement, confinement ratio,
                penetration reach, directional bias).

Both tables carry File_ID (gastruloid) and Group (treatment/control) as columns, so
"gastruloid" and "condition" are just things you group/filter on — not separate tables.

Time/radial binning is NOT precomputed (bin size changes per plot); use `time_summary`
so the binning RULE (file-weighted median + IQR + a min-gastruloids-per-bin cut) is
identical everywhere even though the bin size is passed in per call.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

from tracking.filtering import apply_filters, filters

TRACK = "CLEAN_TRACK_ID"


def load_edges(pkl_path, group_rename=None, exclude_files=None, exclude_contains=None,
               exclude_reason="", base_filter=False, verbose=True):
    """Load a radial-edges pickle, standardize Group, optionally apply the base track filter,
    add dist_from_edge. Returns the edge-level table every edge plot should use.

    base_filter : if True, apply filtering.apply_filters (drop tracks < 20 frames + radial
                  jumps > 200um). DEFAULT FALSE, because your original figure notebooks had
                  this commented out — so leaving it off reproduces your published numbers.
                  Turn it on if you want short/glitchy tracks dropped (changes numbers
                  consistently across every panel).
    exclude_files    : explicit File_IDs to drop.
    exclude_contains : substrings; drop any File_ID containing one (e.g. ['offset_29'] to
                       drop a whole QC-failed batch across BOTH conditions — symmetric, and
                       you can't miss a gastruloid by hand).
    The exclusion is applied HERE, once, so every panel inherits it — never silently per-plot.
    Pass exclude_reason so the justification is recorded (printed + df.attrs)."""
    df = pd.read_pickle(pkl_path)
    if group_rename:
        df["Group"] = df["Group"].replace(group_rename)
    drop = set(exclude_files or [])
    if exclude_contains:
        for fid in df["File_ID"].unique():
            if any(sub in fid for sub in exclude_contains):
                drop.add(fid)
    if drop:
        before = df["File_ID"].nunique()
        df = df[~df["File_ID"].isin(drop)].copy()
        df.attrs["excluded"] = {"files": sorted(drop), "reason": exclude_reason}
        if verbose:
            print(f"Excluded {before - df['File_ID'].nunique()} gastruloid(s): {sorted(drop)}")
            print(f"  reason: {exclude_reason or '(none given — document this!)'}")
    if base_filter:
        df = apply_filters(df, verbose=verbose)          # off by default (matches your originals)
    if {"diameter_um", "radial_distance"} <= set(df.columns):
        df["dist_from_edge"] = df["diameter_um"] / 2.0 - df["radial_distance"]
    return df


def summarize_tracks(edges, track_col=TRACK):
    """
    Collapse edges -> one row per (File_ID, Group, track) with track-level metrics.
    Every track-level panel should read from THIS table so they all share the same cells.

    Columns produced:
      n_steps           number of edges (steps) in the track
      duration_h        time span of the track (h)
      t_start, t_end    first/last EDGE_TIME
      mean_speed        mean per-step SPEED (um/min)
      mean_cos_theta    directional bias; sign per the cos_theta column convention
      net_radial_disp   last - first radial_distance (negative = net inward)
      max_dist_from_edge furthest the cell got from the edge (penetration reach, um)
      confinement_ratio net displacement / path length (1 = straight, ->0 = wandering)
    """
    tc = track_col if track_col in edges.columns else "TRACK_ID"
    df = edges.sort_values(["File_ID", tc, "EDGE_TIME"])
    g = df.groupby(["File_ID", "Group", tc], sort=False)

    out = g.agg(
        n_steps=("EDGE_TIME", "size"),
        t_start=("EDGE_TIME", "min"),
        t_end=("EDGE_TIME", "max"),
        mean_speed=("SPEED", "mean"),
        mean_cos_theta=("cos_theta", "mean"),
        path_len=("DISPLACEMENT", "sum"),
        r_first=("radial_distance", "first"),
        r_last=("radial_distance", "last"),
        x0=("source_x", "first"), y0=("source_y", "first"),
        x1=("target_x", "last"), y1=("target_y", "last"),
    )
    if "dist_from_edge" in df.columns:
        out["max_dist_from_edge"] = g["dist_from_edge"].max()
    out["duration_h"] = out["t_end"] - out["t_start"]
    out["net_radial_disp"] = out["r_last"] - out["r_first"]
    net = np.hypot(out["x1"] - out["x0"], out["y1"] - out["y0"])
    out["confinement_ratio"] = np.where(out["path_len"] > 0, net / out["path_len"], np.nan)
    out = out.drop(columns=["r_first", "r_last", "x0", "y0", "x1", "y1"]).reset_index()
    return out


def time_summary(df, value, bin_hours, *, min_files=3, time_col="EDGE_TIME",
                 file_col="File_ID", group_col="Group"):
    """
    File-weighted time course with ONE shared rule:
      1) bin time into `bin_hours` windows,
      2) per (group, bin) take each gastruloid's median of `value`,
      3) across gastruloids report median + IQR (q25, q75),
      4) keep only bins with >= `min_files` contributing gastruloids.

    Returns tidy DataFrame: Group, bin_mid, median, q25, q75, n_files.
    Pass the same `min_files` everywhere (default 3) and the rule can't drift between plots.
    """
    d = df[np.isfinite(df[value])].copy()
    d["_bin"] = (d[time_col] // bin_hours) * bin_hours + bin_hours / 2.0   # bin center
    per_file = (d.groupby([group_col, "_bin", file_col])[value]
                  .median().reset_index())
    rows = []
    for (grp, b), sub in per_file.groupby([group_col, "_bin"]):
        if sub[file_col].nunique() < min_files:
            continue
        v = sub[value].values
        rows.append({group_col: grp, "bin_mid": b,
                     "median": np.median(v),
                     "q25": np.percentile(v, 25), "q75": np.percentile(v, 75),
                     "n_files": sub[file_col].nunique()})
    return pd.DataFrame(rows).sort_values([group_col, "bin_mid"]).reset_index(drop=True)


def penetration_by_window(per_track, windows, thresh_um=150.0):
    """Per-gastruloid % of tracks reaching >= thresh_um from edge, within each time window.
    Uses max_dist_from_edge from per_track. Returns Group, window, File_ID, percentage."""
    rows = []
    for (g, fid), sub in per_track.groupby(["Group", "File_ID"]):
        for (t0, t1) in windows:
            inwin = sub[(sub["t_start"] < t1) & (sub["t_end"] > t0)]
            if len(inwin) == 0:
                continue
            pct = 100.0 * (inwin["max_dist_from_edge"] >= thresh_um).mean()
            rows.append({"Group": g, "window": f"{t0}-{t1}h", "File_ID": fid, "percentage": pct})
    return pd.DataFrame(rows)
