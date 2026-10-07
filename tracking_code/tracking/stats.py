"""
stats — quick, methods-consistent summary statistics for tracking metrics.

Follows the thesis convention: aggregate each metric PER GASTRULOID (File_ID) first,
then compare groups with a Mann-Whitney U test and report median + IQR. Aggregating per
gastruloid (not per track/per row) keeps n = independent gastruloids, which is the correct
unit of replication.

Usage
-----
    from tracking.stats import summary_stats
    summary_stats(edges, value_cols=["SPEED", "radial_velocity"])          # all timepoints
    summary_stats(edges, value_cols=["SPEED", "radial_velocity"], time_range=(28, 48))
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import stats as _ss


def _per_gastruloid(df, value_col, group_col, file_col, agg):
    """Median (or mean) of value_col per File_ID, carrying its Group label."""
    f = (df.groupby([group_col, file_col])[value_col]
           .agg(agg).reset_index().rename(columns={value_col: "value"}))
    return f


def summary_stats(edges, value_cols=("SPEED", "radial_velocity"),
                  group_col="Group", file_col="File_ID",
                  time_col="EDGE_TIME", time_range=None,
                  per_file_agg="median", test="mannwhitney", verbose=True):
    """
    For each metric: aggregate per gastruloid, then summarize by group and (if exactly two
    groups) run a Mann-Whitney U test. Returns a tidy DataFrame; prints a readable report.
    """
    df = edges
    if time_range is not None:
        lo, hi = time_range
        df = df[(df[time_col] >= lo) & (df[time_col] < hi)]

    rows = []
    for vc in value_cols:
        if vc not in df.columns:
            if verbose: print(f"[skip] '{vc}' not in dataframe")
            continue
        per = _per_gastruloid(df.dropna(subset=[vc]), vc, group_col, file_col, per_file_agg)
        groups = list(per[group_col].unique())
        # per-group summary
        for g in groups:
            v = per.loc[per[group_col] == g, "value"].to_numpy()
            q1, med, q3 = np.percentile(v, [25, 50, 75]) if len(v) else (np.nan,)*3
            rows.append({"metric": vc, "group": g, "n_gastruloids": len(v),
                         "median": med, "IQR_low": q1, "IQR_high": q3,
                         "mean": np.mean(v) if len(v) else np.nan})
        # pairwise test if exactly two groups
        pval = np.nan
        if len(groups) == 2:
            a = per.loc[per[group_col] == groups[0], "value"].to_numpy()
            b = per.loc[per[group_col] == groups[1], "value"].to_numpy()
            if len(a) and len(b):
                if test == "mannwhitney":
                    pval = _ss.mannwhitneyu(a, b, alternative="two-sided").pvalue
                elif test == "ks":
                    pval = _ss.ks_2samp(a, b).pvalue
        for r in rows:
            if r["metric"] == vc and len(groups) == 2:
                r["p_value"] = pval

    out = pd.DataFrame(rows)
    if verbose and not out.empty:
        tw = f" (t in [{time_range[0]}, {time_range[1]}) h)" if time_range else ""
        print(f"Per-gastruloid {per_file_agg} summary{tw}:\n")
        for vc in out["metric"].unique():
            sub = out[out["metric"] == vc]
            print(f"  {vc}:")
            for _, r in sub.iterrows():
                print(f"    {r['group']:<12} n={int(r['n_gastruloids'])}  "
                      f"median={r['median']:.3g}  IQR=[{r['IQR_low']:.3g}, {r['IQR_high']:.3g}]")
            if "p_value" in sub.columns and sub["p_value"].notna().any():
                p = sub["p_value"].dropna().iloc[0]
                print(f"    Mann-Whitney U p = {p:.3g}{'  (*)' if p < 0.05 else ''}")
            print()
    return out
