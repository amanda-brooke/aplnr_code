"""
ps_profile.py — radial profiles + Gaussian core (primary) + FWHM (supplement)
VERBATIM (ast-extracted) from figure_1_clean_3_usethisone_2.ipynb & figure_3_clean_3.ipynb.
Pure data functions; no plotting.
"""
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy import stats as _stats


# ==========================================================================
# GAUSSIAN CORE — PRIMARY METHOD (mu, sigma of the PS-like core)
# ==========================================================================

def gaussian_curve(x, amplitude, mu, sigma):
    return amplitude * np.exp(
        -0.5 * ((x - mu) / sigma) ** 2
    )


def gaussian_core(x, amplitude, mu, sigma):
    return amplitude * np.exp(
        -0.5 * ((x - mu) / sigma) ** 2
    )


def robust_gaussian_core_fit(positions, bins):
    """
    Robustly fit the dominant Gaussian-like peak.

    The soft-L1 loss reduces the influence of the long inward tail.
    """
    positions = np.asarray(positions, dtype=float)

    counts, _ = np.histogram(positions, bins=bins)
    y = counts.astype(float) / counts.sum()
    x = (bins[:-1] + bins[1:]) / 2

    # Only used to identify an initial peak.
    y_smooth = gaussian_filter1d(y, sigma=1)
    peak_i = int(np.argmax(y_smooth))

    amplitude0 = float(y_smooth[peak_i])
    mu0 = float(x[peak_i])

    # Initial sigma estimated from the approximate half-maximum width.
    above_half = np.where(
        y_smooth >= 0.5 * amplitude0
    )[0]

    if len(above_half) >= 2:
        initial_fwhm = (
            x[above_half[-1]] -
            x[above_half[0]]
        )
        sigma0 = max(initial_fwhm / 2.355, 10)
    else:
        sigma0 = 25.0

    p0 = [amplitude0, mu0, sigma0]

    lower = [
        0.0,          # amplitude
        bins[0],      # mu
        5.0,          # sigma
    ]

    upper = [
        1.0,
        bins[-1],
        120.0,
    ]

    # Scale the robust loss relative to the profile height.
    robust_scale = max(0.01, amplitude0 * 0.10)

    fit = least_squares(
        lambda p: gaussian_core(x, *p) - y,
        x0=p0,
        bounds=(lower, upper),
        loss="soft_l1",
        f_scale=robust_scale,
        max_nfev=20000,
    )

    amplitude, mu, sigma = fit.x
    predicted = gaussian_core(
        x,
        amplitude,
        mu,
        sigma,
    )

    # Evaluate fit quality only around the fitted core.
    # We expect the far-right tail not to fit.
    core_mask = predicted >= 0.20 * amplitude

    if core_mask.sum() >= 3:
        ss_res = np.sum(
            (y[core_mask] - predicted[core_mask]) ** 2
        )
        ss_tot = np.sum(
            (y[core_mask] - y[core_mask].mean()) ** 2
        )
        core_r2 = (
            1 - ss_res / ss_tot
            if ss_tot > 0
            else np.nan
        )
    else:
        core_r2 = np.nan

    # ESCAPE IS COMPUTED ONLY IN ps_escape (canonical mu+1.5sigma, median-position profile).
    # The old mu+2sigma escape has been removed here so ps_profile cannot be used for escape.
    # This function returns the CORE FIT only (mu / sigma / variance / core_r2).
    return {
        "x": x,
        "observed": y,
        "predicted": predicted,
        "amplitude": amplitude,
        "mu": mu,
        "sigma": sigma,
        "variance": sigma ** 2,
        "core_r2": core_r2,
    }


def robust_gaussian_profile_fit(
    x,
    y,
    smooth_sigma=1.0,
    min_sigma_um=5.0,
    max_sigma_um=100.0,
):
    """
    Robust Gaussian fit to the dominant density peak.
    The soft-L1 loss reduces the effect of the inward tail.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]

    if len(x) < 5 or np.nanmax(y) <= 0:
        return None

    y_smooth = gaussian_filter1d(
        np.nan_to_num(y, nan=0.0),
        sigma=smooth_sigma,
    )

    peak_idx = int(np.argmax(y_smooth))
    amplitude0 = float(y_smooth[peak_idx])
    mu0 = float(x[peak_idx])

    above_half = np.where(
        y_smooth >= amplitude0 / 2
    )[0]

    if len(above_half) >= 2:
        fwhm0 = x[above_half[-1]] - x[above_half[0]]
        sigma0 = max(fwhm0 / 2.355, min_sigma_um)
    else:
        sigma0 = 20.0

    try:
        result = least_squares(
            lambda p: gaussian_curve(x, *p) - y,
            x0=[amplitude0, mu0, sigma0],
            bounds=(
                [0.0, x.min(), min_sigma_um],
                [max(1.0, amplitude0 * 5), x.max(), max_sigma_um],
            ),
            loss="soft_l1",
            f_scale=max(0.005, amplitude0 * 0.10),
            max_nfev=20000,
        )
    except Exception:
        return None

    amplitude, mu, sigma = result.x
    predicted = gaussian_curve(x, amplitude, mu, sigma)

    core_mask = predicted >= 0.20 * amplitude

    if core_mask.sum() >= 3:
        ss_res = np.sum(
            (y[core_mask] - predicted[core_mask]) ** 2
        )
        ss_tot = np.sum(
            (y[core_mask] - y[core_mask].mean()) ** 2
        )

        core_r2 = (
            1 - ss_res / ss_tot
            if ss_tot > 0
            else np.nan
        )
    else:
        core_r2 = np.nan

    return {
        "amplitude": amplitude,
        "mu": mu,
        "sigma": sigma,
        "left": mu - sigma,
        "right": mu + sigma,
        "predicted": predicted,
        "core_r2": core_r2,
    }


def fit_gaussians_to_profiles(
    x,
    profiles,
    smooth_sigma=1.0,
    minimum_r2=0.60,
):
    records = []

    for file_id, profile in profiles:
        fit = robust_gaussian_profile_fit(
            x,
            profile,
            smooth_sigma=smooth_sigma,
        )

        if fit is None:
            continue

        if minimum_r2 is not None:
            if (
                not np.isfinite(fit["core_r2"])
                or fit["core_r2"] < minimum_r2
            ):
                continue

        records.append(
            {
                "File_ID": file_id,
                **fit,
            }
        )

    return records




# ==========================================================================
# RADIAL DENSITY PROFILES (per gastruloid, per window)
# ==========================================================================

def radial_profiles_by_file_normalized(df, time_range, edge_bins_um, min_valid_frac=0.2):
    t0, t1 = time_range
    df_time = df[(df["EDGE_TIME"] >= t0) & (df["EDGE_TIME"] < t1)].copy()
    df_time = df_time[df_time["TRACK_ID"].notna()].copy()
    profiles = []
    for file_id, file_df in df_time.groupby("File_ID"):
        if ("center_x" not in file_df.columns) or ("center_y" not in file_df.columns):
            continue
        cx = file_df["center_x"].iloc[0]
        cy = file_df["center_y"].iloc[0]
        radius = file_radius_um(file_df, cx, cy)
        r_from_center = np.sqrt((file_df["source_x"] - cx) ** 2 + (file_df["source_y"] - cy) ** 2)
        dist_from_edge = radius - r_from_center
        valid = (dist_from_edge >= edge_bins_um[0]) & (dist_from_edge < edge_bins_um[-1])
        if valid.mean() < min_valid_frac:
            continue
        file_df = file_df.loc[valid].copy()
        file_df["dist_from_edge"] = dist_from_edge[valid].to_numpy()
        file_df["edge_bin"] = pd.cut(file_df["dist_from_edge"], bins=edge_bins_um, right=False, labels=False)
        total_cells = file_df["TRACK_ID"].nunique()
        if total_cells == 0:
            continue
        counts = file_df.groupby("edge_bin")["TRACK_ID"].nunique()
        # within this gastruloid's radius, empty bins are real zeros (0 cells -> 0 density);
        # bins beyond the radius stay NaN (outside the tissue, not measured)
        prof = np.where(np.asarray(edge_bins_um[:-1], float) < radius, 0.0, np.nan)
        for b, c in counts.items():
            if pd.isna(b):
                continue
            prof[int(b)] = c / total_cells
        profiles.append((file_id, prof))
    edge_centers = (edge_bins_um[:-1] + edge_bins_um[1:]) / 2.0
    return edge_centers, profiles


def build_profile_df(df, time_windows_lab, edge_bins_um):
    centers = (edge_bins_um[:-1] + edge_bins_um[1:]) / 2.0
    tcol = "CLEAN_TRACK_ID" if "CLEAN_TRACK_ID" in df.columns else "TRACK_ID"
    rows = []
    for (t0, t1, lab) in time_windows_lab:
        win = df[(df["EDGE_TIME"] >= t0) & (df["EDGE_TIME"] < t1) & df[tcol].notna()]
        for (fid, grp), g in win.groupby(["File_ID", "Group"]):
            cx, cy = g["center_x"].iloc[0], g["center_y"].iloc[0]
            radius = g["diameter_um"].iloc[0] / 2.0
            de = radius - np.sqrt((g["source_x"] - cx) ** 2 + (g["source_y"] - cy) ** 2)
            m = (de >= edge_bins_um[0]) & (de < edge_bins_um[-1])
            gg = g.loc[m].copy(); gg["_de"] = de[m].to_numpy()
            total = gg[tcol].nunique()
            if total == 0:
                continue
            gg["_bin"] = pd.cut(gg["_de"], bins=edge_bins_um, right=False, labels=False)
            cnt = gg.groupby("_bin")[tcol].nunique()
            for b in range(len(centers)):
                rows.append({"File_ID": fid, "Group": grp, "time_bin": lab,
                             "dist_bin_center": float(centers[b]),
                             "density": float(cnt.get(b, 0)) / total})
    return pd.DataFrame(rows)


def summarize_profiles_IQR(profiles):
    if len(profiles) == 0:
        return None, None, None
    stack = np.vstack([p for _, p in profiles])
    return (np.nanmedian(stack, axis=0),
            np.nanpercentile(stack, 25, axis=0),
            np.nanpercentile(stack, 75, axis=0))


def extract_ring_features(profile_df, sigma=1):
    rows = []

    for (file_id, group, tb), sub in profile_df.groupby(
        ["File_ID", "Group", "time_bin"], observed=False
    ):
        sub = sub.sort_values("dist_bin_center")
        x = sub["dist_bin_center"].values.astype(float)
        y = sub["density"].values.astype(float)

        if len(x) < 5 or np.nanmax(y) <= 0:
            continue

        # Smooth density profile
        ys = gaussian_filter1d(y, sigma=sigma)

        # Normalize only for thresholding
        ymax = np.nanmax(ys)
        yn = ys / ymax

        # Peak of densest annulus
        peak_idx = np.nanargmax(yn)
        d_peak = x[peak_idx]

        # Find FWHM around the main peak only
        half = 0.5

        left_idx = peak_idx
        while left_idx > 0 and yn[left_idx] >= half:
            left_idx -= 1

        right_idx = peak_idx
        while right_idx < len(yn) - 1 and yn[right_idx] >= half:
            right_idx += 1

        d_outer_halfmax = x[left_idx]   # closer to edge
        d_inner_halfmax = x[right_idx]  # farther inward

        ring_center = 0.5 * (d_outer_halfmax + d_inner_halfmax)
        ring_width = d_inner_halfmax - d_outer_halfmax

        rows.append({
            "File_ID": file_id,
            "Group": group,
            "time_bin": tb,
            "d_peak": d_peak,
            "d_outer_halfmax": d_outer_halfmax,
            "d_inner_halfmax": d_inner_halfmax,
            "d_ring_center": ring_center,
            "ring_width": ring_width,
            "peak_density": ymax,
        })

    return pd.DataFrame(rows)




# ==========================================================================
# FWHM / HALF-MAX — SECONDARY (supplement only; NOT the escape boundary)
# ==========================================================================

def fwhm_bounds(x, y):
    if y is None or np.all(np.isnan(y)):
        return (np.nan, np.nan)
    y0 = np.nan_to_num(y, nan=0.0)
    if y0.max() <= 0:
        return (np.nan, np.nan)
    above = y0 >= y0.max() / 2.0
    if above.sum() == 0:
        return (np.nan, np.nan)
    return x[above][0], x[above][-1]


def fwhm_center(h, c):
    h = np.asarray(h, float)
    if not np.isfinite(h).any() or h.max() <= 0: return np.nan
    half = h.max() / 2.0
    above = np.where(h >= half)[0]
    if len(above) == 0: return np.nan
    lo, hi = above[0], above[-1]
    left  = c[0]  if lo == 0       else np.interp(half, [h[lo-1], h[lo]], [c[lo-1], c[lo]])
    right = c[-1] if hi == len(h)-1 else np.interp(half, [h[hi+1], h[hi]], [c[hi+1], c[hi]])
    return 0.5 * (left + right)




# --- non-Gaussianity metrics live in ps_escape (they reuse radial_profile);
# re-exported here so they are findable as profile-shape measures too ---
from ps_escape import non_gaussianity, nongauss_by_window  # noqa: F401
