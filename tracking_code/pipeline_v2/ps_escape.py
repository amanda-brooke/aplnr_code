"""
ps_escape.py — PS-like core fitting and cell-escape metrics.

Single source of truth for the escape method (supersedes the various
`calculate_penetration_*` / `dynamic_gaussian_escape` copies scattered in the old
plotter and notebooks). Pure data functions — NO plotting, NO styling.

Method (see manuscript methods / PS_like_project_summary.md):
  1. Build a radial density profile per gastruloid per time window, counting UNIQUE
     CLEAN_TRACK_ID per bin (fraction of that gastruloid's cells).
  2. Fit a Gaussian *core* (composite: Gaussian + one-sided inward exp tail) to a
     reference window -> mu, sigma.  Escape boundary = mu + k*sigma  (k=1.5 default).
  3. Window crossing fraction = fraction of tracks that cross the (fixed) boundary per window.
  4. EXCESS CROSSING (adjusted escape) = observed - prediction from each gastruloid's own migrated/broadened
     core (Gaussian with escape-free edge-side sigma).  Isolates escape from
     region migration + broadening.

Works in either coordinate: distance-from-edge in um, or fraction-of-radius
(pass the appropriate column + bins).  Same function for KD, parent, Sox2, z-plane.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import norm
import warnings
from scipy.optimize import OptimizeWarning
warnings.filterwarnings('ignore', category=OptimizeWarning)

TRACK_ID = "CLEAN_TRACK_ID"


# ----------------------------- models ------------------------------------------
def gaussian(x, A, mu, s):
    return A * np.exp(-(x - mu) ** 2 / (2 * s ** 2))


def gauss_tail_inward(x, A, mu, s, B, tau):
    """Gaussian core + one-sided exponential tail toward the CENTER (larger coord)."""
    return gaussian(x, A, mu, s) + np.where(x > mu, B * np.exp(-(x - mu) / tau), 0.0)


# ----------------------------- profiles ----------------------------------------
def make_bins(lo, hi, n):
    edges = np.linspace(lo, hi, n + 1)
    return edges, 0.5 * (edges[:-1] + edges[1:])


def radial_profile(sub, coord, edges, id_col=TRACK_ID, normalize=True, reduce_tracks=True):
    """
    Fraction-of-cells profile over `coord`.
    reduce_tracks=True (default): collapse each track to ONE position per window (its median
    position) before binning, so a moving cell contributes to a single bin (one cell = one
    vote). reduce_tracks=False reverts to the old edge-based binning (a track spans several
    bins) -- kept only for comparison.
    """
    lo, hi = edges[0], edges[-1]
    ctr = 0.5 * (edges[:-1] + edges[1:])
    if reduce_tracks:
        pos = sub.groupby(id_col)[coord].median()          # one position per track this window
        pos = pos[(pos >= lo) & (pos < hi)]
        tot = int(pos.size)
        if tot == 0:
            return None
        b = pd.cut(pos, bins=edges, right=False, labels=False)
        p = np.zeros(len(ctr))
        for k, c in b.value_counts().items():
            if pd.notna(k):
                p[int(k)] = c / tot
    else:
        gv = sub[(sub[coord] >= lo) & (sub[coord] < hi)]
        tot = gv[id_col].nunique()
        if tot == 0:
            return None
        b = pd.cut(gv[coord], bins=edges, right=False, labels=False)
        p = np.zeros(len(ctr))
        for k, c in gv.groupby(b)[id_col].nunique().items():
            if pd.notna(k):
                p[int(k)] = c / tot
    if normalize and p.max() > 0:
        p = p / p.max()
    return p


def median_profile(df, coord, edges, id_col=TRACK_ID):
    profs = [radial_profile(g, coord, edges, id_col) for _, g in df.groupby("File_ID")]
    profs = [p for p in profs if p is not None]
    return np.median(profs, axis=0) if profs else None


# ----------------------------- fits --------------------------------------------
def fit_core(y, ctr):
    """Composite fit -> (mu, sigma) of the Gaussian core."""
    span = ctr[-1] - ctr[0]
    p0 = [y.max(), ctr[np.argmax(y)], span * 0.12, y.max() * 0.3, span * 0.2]
    lb = [0, ctr[0], span * 0.02, 0, span * 0.03]
    ub = [y.max() * 3, ctr[-1], span * 0.5, y.max() * 3, span * 2]
    p, _ = curve_fit(gauss_tail_inward, ctr, y, p0=p0, bounds=(lb, ub), maxfev=40000)
    return p[1], abs(p[2])


def sigma_edge(y, ctr):
    """Escape-free width: fit Gaussian to the edge-facing flank only (coord <= peak)."""
    ip = int(np.argmax(y))
    m = ctr <= ctr[ip]
    p, _ = curve_fit(gaussian, ctr[m], y[m], p0=[y.max(), ctr[ip], (ctr[-1] - ctr[0]) * 0.12], maxfev=20000)
    return ctr[ip], abs(p[2])


def reference_boundary(df, coord, edges, ref_window, k=1.5, reference_group=None, group_col="Group"):
    """Boundary = mu + k*sigma from the composite fit to the reference-window median profile.
    If reference_group is given (e.g. "Control"), the boundary is fit on THAT condition only,
    then held fixed for measuring escape in every group."""
    t0, t1 = ref_window
    ref = df[(df.EDGE_TIME >= t0) & (df.EDGE_TIME < t1)]
    if reference_group is not None and group_col in ref.columns:
        ref = ref[ref[group_col] == reference_group]
    ctr = 0.5 * (edges[:-1] + edges[1:])
    med = median_profile(ref, coord, edges)
    mu, sg = fit_core(med, ctr)
    return mu, sg, mu + k * sg


# ----------------------------- escape metrics ----------------------------------
def escape_by_window(df, coord, edges, boundary, time_bin_h=2, id_col=TRACK_ID):
    """
    Per gastruloid, per time bin:
      obs        = median-position occupancy (% of tracks whose MEDIAN position is past boundary; matches exp)
      cross_frac = any-edge crossing fraction (% of tracks that EVER cross the boundary this window; descriptive)
      exp    = % predicted by that gastruloid's own migrated+broadened core (no escape)
      excess = max(obs - exp, 0)  = EXCESS CROSSING (adjusted escape)
    Returns tidy DataFrame: File_ID, t_bin, n_tracks, obs, exp, excess.
    """
    ctr = 0.5 * (edges[:-1] + edges[1:])
    d = df.copy()
    d["t_bin"] = (d.EDGE_TIME // time_bin_h) * time_bin_h + time_bin_h / 2
    rows = []
    for (fid, tb), s in d.groupby(["File_ID", "t_bin"]):
        tot = s[id_col].nunique()
        if tot == 0:
            continue
        # one position per track (median) so obs matches the profile/exp basis (fixes the
        # obs-vs-exp definition mismatch): fraction of tracks whose MEDIAN position is past boundary
        pos = s.groupby(id_col)[coord].median()
        obs = 100 * float((pos >= boundary).mean())            # median-position occupancy (used in excess)
        cross_frac = 100 * s[s[coord] >= boundary][id_col].nunique() / tot  # any-edge crossing (descriptive)
        prof = radial_profile(s, coord, edges, id_col)
        exp = np.nan
        if prof is not None and prof.max() > 0:
            try:
                pk, se = sigma_edge(prof, ctr)
                exp = 100 * (1 - norm.cdf((boundary - pk) / se))
            except Exception:
                pass
        rows.append(dict(File_ID=fid, t_bin=tb, n_tracks=tot,
                         obs=obs, cross_frac=cross_frac, exp=exp,
                         excess=max(obs - exp, 0) if np.isfinite(exp) else np.nan))
    return pd.DataFrame(rows)


# def run_escape(df, coord, edges, ref_window, k=1.5, group_col=None, reference_group=None):
#     """Convenience: reference boundary + per-window escape (optionally carrying a group column).
#     reference_group (e.g. "Control") fits the boundary on that condition only, then measures
#     escape for every group against that fixed boundary."""
#     mu, sg, boundary = reference_boundary(df, coord, edges, ref_window, k,
#                                           reference_group=reference_group,
#                                           group_col=(group_col or "Group"))
#     res = escape_by_window(df, coord, edges, boundary)
#     if group_col and group_col in df.columns:
#         gmap = df.groupby("File_ID")[group_col].first()
#         res["group"] = res["File_ID"].map(gmap)
#     res.attrs.update(dict(mu=mu, sigma=sg, boundary=boundary, k=k, coord=coord))
#     return res, dict(mu=mu, sigma=sg, boundary=boundary)

def run_escape(
    df, coord, edges, ref_window,
    k=1.5,
    group_col=None,
    reference_group=None,
    time_bin_h=2,
):
    """Reference boundary + per-window crossing metrics."""
    
    mu, sg, boundary = reference_boundary(
        df, coord, edges, ref_window, k,
        reference_group=reference_group,
        group_col=(group_col or "Group"),
    )

    res = escape_by_window(
        df, coord, edges, boundary,
        time_bin_h=time_bin_h,
    )

    if group_col and group_col in df.columns:
        gmap = df.groupby("File_ID")[group_col].first()
        res["group"] = res["File_ID"].map(gmap)

    res.attrs.update(
        dict(
            mu=mu,
            sigma=sg,
            boundary=boundary,
            k=k,
            coord=coord,
            time_bin_h=time_bin_h,
        )
    )

    return res, dict(mu=mu, sigma=sg, boundary=boundary)


# ----------------------------- non-Gaussianity (inward-persistence) ------------
def non_gaussianity(prof, ctr):
    """
    Deviation of a density profile from Gaussian (Alex's 'inward persistence' idea).
    Returns: skew (signed, +=inward tail), kurt (excess), kl (KL from a Gaussian matched
    to the profile's own mean+variance), inward_excess (% mass beyond that Gaussian, inward side).
    """
    p = np.nan_to_num(prof, nan=0.0).astype(float)
    if p.sum() <= 0:
        return dict(skew=np.nan, kurt=np.nan, kl=np.nan, inward_excess=np.nan)
    p = p / p.sum()
    m = float(np.sum(p * ctr)); v = float(np.sum(p * (ctr - m) ** 2))
    if v <= 0:
        return dict(skew=np.nan, kurt=np.nan, kl=np.nan, inward_excess=np.nan)
    sd = np.sqrt(v)
    g = np.exp(-(ctr - m) ** 2 / (2 * v)); g = g / g.sum()
    eps = 1e-9
    return dict(
        skew=float(np.sum(p * ((ctr - m) / sd) ** 3)),
        kurt=float(np.sum(p * ((ctr - m) / sd) ** 4) - 3),
        kl=float(np.sum(p * np.log((p + eps) / (g + eps)))),
        inward_excess=float(100 * np.sum(np.clip(p - g, 0, None)[ctr > m])),
    )


def nongauss_by_window(df, coord, edges, time_bin_h=2, id_col=TRACK_ID, group_col=None):
    """Per gastruloid, per time bin: skew / kurt / kl / inward_excess of the profile."""
    ctr = 0.5 * (edges[:-1] + edges[1:])
    d = df.copy()
    d["t_bin"] = (d.EDGE_TIME // time_bin_h) * time_bin_h + time_bin_h / 2
    rows = []
    for (fid, tb), s in d.groupby(["File_ID", "t_bin"]):
        prof = radial_profile(s, coord, edges, id_col, normalize=False)
        if prof is None:
            continue
        r = non_gaussianity(prof, ctr)
        r.update(File_ID=fid, t_bin=tb, n_tracks=s[id_col].nunique())
        rows.append(r)
    res = pd.DataFrame(rows)
    if group_col and group_col in df.columns and len(res):
        res["group"] = res["File_ID"].map(df.groupby("File_ID")[group_col].first())
    return res
