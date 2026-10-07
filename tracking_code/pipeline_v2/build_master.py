"""
build_master.py — assemble ONE tidy master measurements table from every source pickle.

Contract: MASTER_SCHEMA.md.  One row per cell-per-frame, every subset as a column,
exclusions as GLOBAL flags (never deletions), coordinates computed once.  Output is
partitioned Parquet you load with `load_master()` (in dataconfig.py).

    python build_master.py                    # all active datasets in analysis.yml
    python build_master.py parent sox2        # just some (by config name)
    python build_master.py --data-dir /path/to/pickles --out data/processed

Writes:
    <out>/master/<experiment>.parquet         per-frame master (partitioned by experiment)
    <out>/master_tracks.parquet               per-track summaries
    <out>/boundaries.parquet                  mu/sigma/mu+k*sigma per dataset & k
    <out>/build_report.txt                    reconciliation log (counts, unmapped cols)

Nothing is dropped: the build STOPS if a source column is neither mapped nor listed
in COLUMNS_INTENTIONALLY_DROPPED, so a stray measurement can't silently vanish.
"""
from __future__ import annotations
import sys, os, argparse, datetime, textwrap
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dataconfig as dc
import ps_escape as pe

def _save(df, path_noext):
    """Write parquet if an engine is available, else fall back to pickle. Returns the path written."""
    try:
        p = path_noext + ".parquet"; df.to_parquet(p, index=False); return p
    except Exception:
        p = path_noext + ".pkl"; df.to_pickle(p); return p


KS = (1.0, 1.5, 2.0)                      # boundary multipliers to precompute
COLUMNS_INTENTIONALLY_DROPPED = {         # known source cols we deliberately don't carry
    "MANUAL_EDGE_COLOR",
    "EDGE_X_LOCATION", "EDGE_Y_LOCATION", "EDGE_Z_LOCATION",
    "diameter_bbox_um", "diameter_pctl_um", "diameter_firstframe_um", "diameter_boundary_um",
}
# extra source cols carried because notebook panels reference them
EXTRA_KEEP = ["diameter_source", "target_x", "target_y"]
# JOIN KEYS back to the spot-level table (per-spot mean fluorescence intensity lives there).
# Keep these so intensity can be merged in later without rebuilding.
SPOT_KEYS = ["SPOT_SOURCE_ID", "SPOT_TARGET_ID"]
# raw per-frame measurements carried straight through (present-if-there)
RAW_MEAS = ["SPEED", "DISPLACEMENT", "radial_velocity", "radial_velocity_theta",
            "cos_theta", "radial_persistence", "angle_to_center",
            "DIRECTIONAL_CHANGE_RATE", "LINK_COST", "u", "v", "LABEL", "TRACK_ID"]
COORD_KEEP = ["EDGE_TIME", "source_x", "source_y", "center_x", "center_y",
              "radial_distance", "diameter_um"]


# --------------------------------------------------------------------------- labels
def _apply_labels(df, cfg):
    """Attach experiment / cell_type / geometry / condition / z_stack from cfg rules."""
    fid = df["File_ID"].astype(str)
    df["experiment"] = cfg["experiment"]
    df["cell_type"]  = cfg["cell_type"]
    df["z_stack"]    = bool(cfg.get("z_stack", False))

    if "geometry" in cfg:
        df["geometry"] = cfg["geometry"]
    elif "geometry_from_substring" in cfg:
        rule = cfg["geometry_from_substring"]; default = rule.get("default")
        geo = pd.Series(default, index=df.index)
        for sub, val in rule.items():
            if sub == "default": continue
            geo = geo.mask(fid.str.contains(sub, na=False), val)
        df["geometry"] = geo
    else:
        df["geometry"] = "unknown"

    if "condition" in cfg:
        df["condition"] = cfg["condition"]
    elif cfg.get("condition_from"):
        df["condition"] = df[cfg["condition_from"]].astype(str)
    elif "condition_from_substring" in cfg:
        rule = cfg["condition_from_substring"]; default = rule.get("default")
        con = pd.Series(default, index=df.index)
        for sub, val in rule.items():
            if sub == "default": continue
            con = con.mask(fid.str.contains(sub, na=False), val)
        df["condition"] = con
    else:
        df["condition"] = "NA"

    # z-plane + parent link (best-effort; refine per real z File_IDs)
    if df["z_stack"].iloc[0]:
        df["z_plane"] = _z_plane_from_fileid(fid)
        df["parent_File_ID"] = fid.str.replace(r"[_-]?(z\d+|top|middle|bottom|apical|mid|basal).*$",
                                               "", regex=True, case=False)
    else:
        df["z_plane"] = pd.NA
        df["parent_File_ID"] = pd.NA
    return df


def _z_plane_from_fileid(fid):
    low = fid.str.lower()
    out = pd.Series(pd.NA, index=fid.index, dtype="object")
    out = out.mask(low.str.contains("top|apical", na=False), "Apical")
    out = out.mask(low.str.contains("mid", na=False), "Mid")
    out = out.mask(low.str.contains("bottom|basal", na=False), "Basal")
    # merged upper stacks (2-plane datasets: shaplnr "toplayers", parent "upper"):
    # a slab of everything above the basal plane, NOT a true apical plane -> label "Upper".
    # Applied last so it overrides the Apical match that "toplayers" would otherwise get,
    # while a bare "top" (the parental 3-plane top plane) stays "Apical".
    out = out.mask(low.str.contains("upper|toplayers", na=False), "Upper")
    return out


# ----------------------------------------------------------------- coordinates / physics
def _add_coords(df):
    df["radius"] = df["diameter_um"] / 2.0
    df["dist_from_edge"] = df["radius"] - df["radial_distance"]
    df["dnorm"] = 1.0 - df["radial_distance"] / df["radius"]
    # tangential velocity: component of (u,v) perpendicular to the radial direction
    if {"u", "v"}.issubset(df.columns):
        dx = df["source_x"] - df["center_x"]; dy = df["source_y"] - df["center_y"]
        r = np.sqrt(dx**2 + dy**2).replace(0, np.nan)
        df["tangential_velocity"] = (-dy * df["u"] + dx * df["v"]) / r
    else:
        df["tangential_velocity"] = np.nan
    return df


# ------------------------------------------------------------------------- exclusions
def _flag_exclusions(df, cfg):
    """GLOBAL flags only. excluded=True + exclude_reason; NOTHING dropped."""
    excluded = pd.Series(False, index=df.index)
    reason = pd.Series("", index=df.index, dtype="object")
    fid = df["File_ID"].astype(str)
    for off in (cfg.get("exclude_offsets") or []):
        hit = fid.str.contains(off, na=False)
        reason = reason.mask(hit & (reason == ""), off)
        excluded = excluded | hit
    sf = cfg.get("size_filter")
    if sf:
        for col, vals in sf.items():
            hit = ~df[col].isin(vals)
            reason = reason.mask(hit & (reason == ""), f"size_filter_{col}")
            excluded = excluded | hit
    df["excluded"] = excluded
    df["exclude_reason"] = reason.replace("", pd.NA)
    return df


# ------------------------------------------------------------------------- boundaries
def _fit_boundaries(df_kept, cfg, name):
    """Fit mu/sigma once (on non-excluded rows) in the dataset's preferred coord; boundary for each k.
    If cfg['reference_group'] is set, fit ONLY on that condition (e.g. Control) — the reference arm."""
    coord, edges = dc.coord_and_edges(cfg)
    ref = df_kept
    rg = cfg.get("reference_group")
    if rg and "condition" in df_kept.columns:
        ref = df_kept[df_kept["condition"] == rg]
        print(f"    [boundary fit on reference_group='{rg}': {ref['File_ID'].nunique()} gastruloids]")
    rows = []
    try:
        mu, sg, _ = pe.reference_boundary(ref, coord, edges,
                                          tuple(cfg["reference_window"]), 1.0)
        for k in KS:
            rows.append(dict(dataset=name, experiment=cfg["experiment"], coord=coord,
                             k=k, mu=mu, sigma=sg, boundary=mu + k * sg))
    except Exception as e:
        print(f"    [boundary fit failed for {name}: {type(e).__name__}: {e}]")
    return coord, rows


# ---------------------------------------------------------------------- per-track table
def _track_summaries(df, coord):
    """One row per (File_ID, CLEAN_TRACK_ID): duration, path, straightness, max dist, net radial..."""
    df = df.sort_values("EDGE_TIME")
    g = df.groupby(["experiment", "File_ID", "CLEAN_TRACK_ID"], sort=False)
    first = g.first(); last = g.last()
    path_len = g["DISPLACEMENT"].sum() if "DISPLACEMENT" in df.columns else g.size() * np.nan
    net = np.sqrt((last["source_x"] - first["source_x"])**2 +
                  (last["source_y"] - first["source_y"])**2)
    out = pd.DataFrame({
        "condition":  first["condition"].values,
        "cell_type":  first["cell_type"].values,
        "geometry":   first["geometry"].values,
        "excluded":   g["excluded"].max().values,          # excluded if any row is
        "duration_h": (last["EDGE_TIME"] - first["EDGE_TIME"]).values,
        "path_len":   path_len.values,
        "net_disp":   net.values,
        "mean_speed": g["SPEED"].mean().values if "SPEED" in df.columns else np.nan,
        "max_dist_edge": g["dist_from_edge"].max().values,
        "max_coord":  g[coord].max().values,
        "med_coord":  g[coord].median().values,
        "net_radial_disp": (last["dist_from_edge"] - first["dist_from_edge"]).values,
    }, index=first.index).reset_index()
    out["straightness"] = out["net_disp"] / out["path_len"].replace(0, np.nan)
    return out


# ------------------------------------------------------------------------------- build
def build(names=None, data_dir="data", out="data/processed"):
    cfgs = dc.load_configs()
    names = names or [n for n in cfgs if isinstance(cfgs[n], dict) and "experiment" in cfgs[n]]
    frame_parts, track_parts, bound_parts, report = [], [], [], []
    fitted = {}   # dataset -> {mu, sigma, coord}; lets other datasets inherit a boundary

    for name in names:
        cfg = cfgs.get(name)
        if not cfg or "experiment" not in cfg:
            report.append(f"{name:26s} SKIP (no experiment label — reserved stub?)"); continue
        path = os.path.join(data_dir, cfg["radial_file"])
        if not os.path.exists(path):
            report.append(f"{name:26s} SKIP (missing file {cfg['radial_file']})"); continue

        raw = pd.read_pickle(path)
        n_rows0, n_files0 = len(raw), raw["File_ID"].nunique()

        # --- column-map audit: every source col mapped, or explicitly dropped, else STOP ---
        known = (set(COORD_KEEP) | set(RAW_MEAS) | set(SPOT_KEYS) | set(EXTRA_KEEP)
                 | {"File_ID", "Group", "CLEAN_TRACK_ID"}
                 | COLUMNS_INTENTIONALLY_DROPPED)
        unmapped = [c for c in raw.columns if c not in known]
        if unmapped:
            raise SystemExit(f"[{name}] UNMAPPED source columns (map them or add to "
                             f"COLUMNS_INTENTIONALLY_DROPPED): {unmapped}")

        df = raw.copy()
        df = _apply_labels(df, cfg)
        df = _add_coords(df)
        df = _flag_exclusions(df, cfg)
        df["replicate"] = pd.NA                        # stub until File_ID->replicate map supplied
        df["source_file"] = cfg["radial_file"]
        df["dataset"] = name
        df["preferred_coord"] = dc.coord_and_edges(cfg)[0]

        # boundaries: fit on kept rows, OR inherit another dataset's (boundary_from)
        coord = dc.coord_and_edges(cfg)[0]
        bf = cfg.get("boundary_from")
        if bf:
            if bf not in fitted:
                raise SystemExit(f"[{name}] boundary_from='{bf}' not built in this run — "
                                 f"include it (or run all datasets).")
            src = fitted[bf]
            if src["coord"] != coord:
                raise SystemExit(f"[{name}] boundary_from='{bf}' unit mismatch "
                                 f"({src['coord']} vs {coord})")
            brows = [dict(dataset=name, experiment=cfg["experiment"], coord=coord, k=k,
                          mu=src["mu"], sigma=src["sigma"], boundary=src["mu"] + k * src["sigma"],
                          boundary_from=bf) for k in KS]
            print(f"    [boundary inherited from '{bf}': mu+1.5sigma={src['mu']+1.5*src['sigma']:.3g}]")
        else:
            _, brows = _fit_boundaries(df[~df["excluded"]], cfg, name)
        if brows:
            fitted[name] = dict(mu=brows[0]["mu"], sigma=brows[0]["sigma"], coord=coord)
        bound_parts += brows
        bmap = {r["k"]: r["boundary"] for r in brows}
        for k in KS:
            col = f"beyond_{str(k).replace('.', '_')}"
            df[col] = (df[coord] >= bmap[k]) if k in bmap else pd.NA

        # keep a tidy column set
        keep = (["dataset", "experiment", "cell_type", "geometry", "condition", "z_stack",
                 "z_plane", "parent_File_ID", "File_ID", "CLEAN_TRACK_ID", "replicate"]
                + SPOT_KEYS
                + COORD_KEEP + ["radius", "dist_from_edge", "dnorm", "tangential_velocity",
                                "preferred_coord"]
                + [c for c in RAW_MEAS if c in df.columns]
                + [c for c in EXTRA_KEEP if c in df.columns]
                + [f"beyond_{str(k).replace('.', '_')}" for k in KS]
                + ["excluded", "exclude_reason", "source_file"])
        df = df[[c for c in keep if c in df.columns]]

        # per-track summaries (+ dispersed label = ever crossed 1.5 boundary)
        tsum = _track_summaries(df, coord)
        b15 = bmap.get(1.5, np.nan)
        tsum["ever_beyond_1_5"] = tsum["max_coord"] >= b15
        track_parts.append(tsum)
        frame_parts.append(df)

        # reconciliation for this dataset
        assert len(df) == n_rows0, f"{name}: row count changed {n_rows0}->{len(df)}"
        assert df["File_ID"].nunique() == n_files0, f"{name}: lost a File_ID"
        nax = int(df["excluded"].sum())
        report.append(f"{name:26s} rows={n_rows0:>9,} files={n_files0:>3} "
                      f"excluded_rows={nax:>8,} boundary1.5={bmap.get(1.5, float('nan')):.3g}")

    if not frame_parts:
        print("nothing built.\n" + "\n".join(report)); return

    master = pd.concat(frame_parts, ignore_index=True)
    tracks = pd.concat(track_parts, ignore_index=True)
    bounds = pd.DataFrame(bound_parts)

    os.makedirs(os.path.join(out, "master"), exist_ok=True)
    _fmt = []
    for exp, part in master.groupby("experiment"):
        _fmt.append(_save(part, os.path.join(out, "master", str(exp))))
    _save(tracks, os.path.join(out, "master_tracks"))
    _save(bounds, os.path.join(out, "boundaries"))
    if _fmt and _fmt[0].endswith(".pkl"):
        print("[note] pyarrow unavailable -> wrote PICKLE (.pkl). load_master reads it fine.")

    # --- global reconciliation + report ---
    total_src = sum(int(l.split("rows=")[1].split()[0].replace(",", ""))
                    for l in report if "rows=" in l)
    header = [f"MASTER BUILD  {datetime.datetime.now():%Y-%m-%d %H:%M}",
              f"master rows: {len(master):,}   tracks: {len(tracks):,}   "
              f"experiments: {sorted(master['experiment'].unique())}",
              f"row reconciliation: sum(sources)={total_src:,}  master={len(master):,}  "
              f"{'OK' if total_src == len(master) else 'MISMATCH!!'}",
              "-" * 72]
    txt = "\n".join(header + report)
    with open(os.path.join(out, "build_report.txt"), "w") as f:
        f.write(txt + "\n")
    print(txt)
    print(f"\nWROTE {out}/master/<experiment>.parquet, master_tracks.parquet, boundaries.parquet")
    return master, tracks, bounds


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("names", nargs="*", help="dataset config names (default: all labeled)")
    ap.add_argument("--data-dir", default="data")
    ap.add_argument("--out", default="data/processed")
    a = ap.parse_args()
    build(a.names or None, data_dir=a.data_dir, out=a.out)
