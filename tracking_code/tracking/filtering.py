"""
filtering — ONE place that defines how tracks are filtered, applied ONCE after loading.

The point: track filtering (dropping short tracks, removing implausible tracking jumps)
should be a single global step, not something each plot does ad hoc. Apply `apply_filters`
right after you load a dataset, then every downstream plot uses the same cleaned data by
construction — so you never again have to wonder "which plot filtered short tracks?"

Two knobs, both reported so the filtering is auditable:
  - min_track_frames        : drop tracks with fewer than N frames (short-lived tracks)
  - max_step_radial_jump_um : drop individual steps where radial position jumps more than
                              this between consecutive frames (TrackMate linking artifacts)

Usage
-----
    from tracking.filtering import filters, apply_filters

    filters.min_track_frames        = 20    # ~3.3 h at a 10-min interval
    filters.max_step_radial_jump_um = 200

    edges = apply_filters(edges)            # filtered once; all plots use `edges` after this

A QC plot that intentionally shows what was excluded can call with overrides, e.g.
`apply_filters(edges_raw, min_track_frames=0)` — deviations are then explicit, never silent.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np


@dataclass
class FilterConfig:
    min_track_frames: int = 20            # drop tracks shorter than this many frames
    max_step_radial_jump_um: float = 200  # drop steps with |Δ radial position| above this
    track_col: str = "CLEAN_TRACK_ID"     # falls back to TRACK_ID if absent
    time_col: str = "EDGE_TIME"
    radial_col: str = "radial_distance"


# the single shared instance notebooks import and tweak
filters = FilterConfig()


def _track_col(df, cfg):
    if cfg.track_col in df.columns:
        return cfg.track_col
    if "TRACK_ID" in df.columns:
        return "TRACK_ID"
    raise KeyError(f"No track id column ({cfg.track_col} or TRACK_ID) in dataframe")


def apply_filters(edges, cfg: FilterConfig = filters, *,
                  min_track_frames: int | None = None,
                  max_step_radial_jump_um: float | None = None,
                  verbose: bool = True):
    """
    Apply the global track filters to an edges dataframe and return the filtered copy.

    Per-call overrides (min_track_frames, max_step_radial_jump_um) let a specific plot
    deviate explicitly; otherwise the shared config values are used. Prints before/after
    counts when verbose. Attaches the applied criteria to df.attrs['filter_record'] so a
    figure can report exactly what was done.
    """
    mtf = cfg.min_track_frames if min_track_frames is None else min_track_frames
    mj  = cfg.max_step_radial_jump_um if max_step_radial_jump_um is None else max_step_radial_jump_um
    tcol = _track_col(edges, cfg)
    df = edges.copy()

    n_rows0 = len(df)
    n_trk0 = df.groupby(["File_ID", tcol]).ngroups if "File_ID" in df.columns else df[tcol].nunique()
    record = {"min_track_frames": mtf, "max_step_radial_jump_um": mj,
              "track_col": tcol, "rows_before": n_rows0, "tracks_before": n_trk0}

    # 1) radial-jump removal (per step within a track), if we have radial positions
    n_jump = 0
    if mj is not None and cfg.radial_col in df.columns:
        grp = ["File_ID", tcol] if "File_ID" in df.columns else [tcol]
        df = df.sort_values(grp + [cfg.time_col])
        dr = df.groupby(grp)[cfg.radial_col].diff()
        jump = dr.abs() > mj            # NaN (first frame of a track) stays False -> kept
        n_jump = int(jump.sum())
        df = df[~jump.fillna(False)]
    elif mj is not None:
        if verbose:
            print(f"  (skipped radial-jump filter: no '{cfg.radial_col}' column)")

    # 2) short-track removal (after jump removal, so frame counts reflect kept steps)
    n_short_tracks = 0
    if mtf:
        grp = ["File_ID", tcol] if "File_ID" in df.columns else [tcol]
        sizes = df.groupby(grp)[cfg.time_col].transform("size")
        keep = sizes >= mtf
        n_short_tracks = (df[~keep].groupby(grp).ngroups
                          if (~keep).any() else 0)
        df = df[keep]

    n_rows1 = len(df)
    n_trk1 = df.groupby(["File_ID", tcol]).ngroups if "File_ID" in df.columns else df[tcol].nunique()
    record.update({"rows_after": n_rows1, "tracks_after": n_trk1,
                   "steps_dropped_jump": n_jump, "short_tracks_dropped": n_short_tracks})
    df.attrs["filter_record"] = record

    if verbose:
        print("Track filtering applied:")
        print(f"  short-track cutoff   : < {mtf} frames  -> dropped {n_short_tracks} tracks")
        print(f"  radial-jump cutoff   : |Δr| > {mj} um   -> dropped {n_jump} steps")
        print(f"  tracks : {n_trk0:,} -> {n_trk1:,}")
        print(f"  rows   : {n_rows0:,} -> {n_rows1:,}")
    return df


def filter_caption(df) -> str:
    """Short human-readable description of the filter applied, for a figure legend/methods."""
    r = df.attrs.get("filter_record")
    if not r:
        return "no filter record"
    return (f"tracks ≥ {r['min_track_frames']} frames; "
            f"steps with |Δr| > {r['max_step_radial_jump_um']} µm removed "
            f"(n = {r['tracks_after']:,} tracks)")
