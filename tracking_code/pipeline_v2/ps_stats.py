"""
ps_stats.py — significance + summary helpers
VERBATIM (ast-extracted) from figure_1_clean_3_usethisone_2.ipynb & figure_3_clean_3.ipynb.
Pure data functions; no plotting.
"""
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy import stats as _stats


# ==========================================================================
# SIGNIFICANCE
# ==========================================================================

def bh_fdr(pvals):                      # Benjamini–Hochberg (no statsmodels needed)
    p = np.asarray(pvals, float); n = len(p); order = np.argsort(p)
    adj = np.empty(n); prev = 1.0
    for i in range(n - 1, -1, -1):
        prev = min(prev, p[order[i]] * n / (i + 1)); adj[order[i]] = prev
    return adj


def sig_bracket(ax, x1, x2, data, text="***", lw=0.8, fontsize=22, frac=0.06):
    """Place bracket just above the data and auto-expand the y-axis so it
    doesn't collide with the top spine (no manual offsets needed)."""
    y0, y1 = ax.get_ylim(); yr = y1 - y0
    top = np.nanmax(np.concatenate([np.asarray(d, float) for d in data]))
    y = top + frac * yr
    h = 0.03 * yr
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color="black", lw=lw)
    ax.text((x1 + x2) / 2, y + h, text, ha="center", va="bottom",
            fontsize=fontsize, fontweight="bold")
    ax.set_ylim(y0, max(y1, y + h + 0.12 * yr))




# ==========================================================================
# BINNING / VALID-WINDOW HELPERS
# ==========================================================================

def shared_valid_bins(
    per_file_a, per_file_b, min_a=8, min_b=3, time_bin_col="time_bin"
):
    """
    Returns bins that have at least min_a files in dataset A and min_b in dataset B.
    Also returns a counts table for inspection.
    """
    cA = bin_counts(per_file_a, time_bin_col).rename("n_A")
    cB = bin_counts(per_file_b, time_bin_col).rename("n_B")
    counts = pd.concat([cA, cB], axis=1).fillna(0).astype(int).reset_index()
    valid = counts.loc[
        (counts["n_A"] >= min_a) & (counts["n_B"] >= min_b), time_bin_col
    ]
    return valid.to_numpy(), counts


def shared_bins(df1, df2, time_bin_col="time_bin"):
    b1 = set(df1[time_bin_col].dropna().unique())
    b2 = set(df2[time_bin_col].dropna().unique())
    return np.array(sorted(b1.intersection(b2)))




def time_summary(df, value, bin_hours, *, min_files=3, time_col="EDGE_TIME",
                 file_col="File_ID", group_col="Group", track_col="CLEAN_TRACK_ID"):
    """
    Group-aware time course — the ONE shared rule (matches Fig 1's file_weighted helper):
      1) bin time,
      2) collapse each TRACK to its mean in (group,bin,file)  -> one cell = one vote,
      3) per (group,bin,file) take the gastruloid MEDIAN across its tracks,
      4) across gastruloids report median + IQR (q25,q75),
      5) keep bins with >= min_files gastruloids.
    Returns tidy: Group, bin_mid, median, q25, q75, n_files. (b) copy in pipeline_v2 so
    tracking/prep.py stays frozen.
    """
    d = df[np.isfinite(df[value])].copy()
    d["_bin"] = (d[time_col] // bin_hours) * bin_hours + bin_hours / 2.0        # bin center
    if track_col in d.columns:                                                  # 2. per-track mean
        d = d.groupby([group_col, "_bin", file_col, track_col])[value].mean().reset_index()
    per_file = d.groupby([group_col, "_bin", file_col])[value].median().reset_index()  # 3. gastruloid median
    rows = []
    for (grp, b), sub in per_file.groupby([group_col, "_bin"]):
        if sub[file_col].nunique() < min_files:
            continue
        v = sub[value].values
        rows.append({group_col: grp, "bin_mid": b, "median": np.median(v),
                     "q25": np.percentile(v, 25), "q75": np.percentile(v, 75),
                     "n_files": sub[file_col].nunique()})
    return pd.DataFrame(rows).sort_values([group_col, "bin_mid"]).reset_index(drop=True)
