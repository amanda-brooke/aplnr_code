"""
stats_export — single source of truth for the numbers a figure reports.

The problem this solves: right now every statistic in the manuscript (means,
SEMs, p-values, n's) is typed by hand into the .tex. Re-run the analysis and a
value changes, and you have to hunt through prose to update it. This module lets
a figure-producing cell *record* every number it reports into a small JSON file
next to the figure's panels. `build_macros.py` then turns those JSONs into a
LaTeX file of \\newcommand macros, so the manuscript references a macro instead
of a literal. Re-run -> regenerate -> recompile, and every number updates itself.

Typical use in a figure notebook
--------------------------------
    from tracking.stats_export import FigureStats

    fs = FigureStats("fig3", outroot="../manuscript/figures")   # writes figures/fig3/stats.json
    # ... after you compute the values you're going to report ...
    fs.add_pm("penet_control", control.mean(), control.std(ddof=1))  # "8.9 \\pm 12.5"
    fs.add_pm("penet_kd",      kd.mean(),      kd.std(ddof=1))
    fs.add_p ("penet_p",       p_value)                              # "1.0 \\times 10^{-6}"
    fs.add_num("penet_U",      U, sig=3)
    fs.save()

The values stored are already LaTeX-ready strings, so the build step is a dumb,
reliable mapping with no formatting logic of its own.
"""
from __future__ import annotations

import json
import math
import datetime
from pathlib import Path


# ---------------------------------------------------------------------------
# Formatters — produce LaTeX-ready strings (meant to sit inside math mode $...$)
# ---------------------------------------------------------------------------
def fmt_sig(x: float, sig: int = 3) -> str:
    """Round to `sig` significant figures and render without scientific notation."""
    x = float(x)
    if x == 0 or not math.isfinite(x):
        return "0"
    decimals = sig - int(math.floor(math.log10(abs(x)))) - 1
    decimals = max(decimals, 0)
    return f"{round(x, decimals):.{decimals}f}"


def pm(mean: float, err: float, sig: int = 3) -> str:
    """Mean +/- error, e.g. '0.536 \\pm 0.065'. Use inside math mode."""
    return f"{fmt_sig(mean, sig)} \\pm {fmt_sig(err, sig)}"


def pm_pct(mean: float, err: float, sig: int = 3) -> str:
    """Percentage mean +/- error, e.g. '8.9\\% \\pm 12.5\\%'. Use inside math mode."""
    return f"{fmt_sig(mean, sig)}\\% \\pm {fmt_sig(err, sig)}\\%"


def pval(p: float, floor: float = 1e-10) -> str:
    """
    LaTeX p-value:
      - below `floor`            -> '< 10^{-10}'
      - below 1e-3               -> 'a.b \\times 10^{-n}'
      - otherwise                -> plain (e.g. '0.004', '0.12')
    Use inside math mode.
    """
    p = float(p)
    if p < floor:
        exp = int(round(math.log10(floor)))
        return f"< 10^{{{exp}}}"
    if p < 1e-3:
        exp = int(math.floor(math.log10(p)))
        mant = p / 10 ** exp
        return f"{mant:.1f} \\times 10^{{{exp}}}"
    return f"{p:.3g}"


# ---------------------------------------------------------------------------
# Recorder
# ---------------------------------------------------------------------------
class FigureStats:
    """Collects LaTeX-ready stat strings for one figure and writes stats.json."""

    def __init__(self, figure_id: str, outroot="../manuscript/figures"):
        self.figure_id = figure_id
        self.outroot = Path(outroot)
        self.values: dict[str, str] = {}
        self.meta: dict[str, object] = {}

    # -- recording ---------------------------------------------------------
    def add(self, key: str, latex: str) -> str:
        """Store a pre-formatted LaTeX string. Returns it (handy for inline display)."""
        self.values[key] = str(latex)
        return self.values[key]

    def add_pm(self, key: str, mean: float, err: float, sig: int = 3) -> str:
        return self.add(key, pm(mean, err, sig))

    def add_pm_pct(self, key: str, mean: float, err: float, sig: int = 3) -> str:
        return self.add(key, pm_pct(mean, err, sig))

    def add_p(self, key: str, p: float) -> str:
        return self.add(key, pval(p))

    def add_num(self, key: str, x: float, sig: int = 3) -> str:
        return self.add(key, fmt_sig(x, sig))

    def set_meta(self, **kw) -> "FigureStats":
        self.meta.update(kw)
        return self

    # -- output ------------------------------------------------------------
    def figdir(self) -> Path:
        d = self.outroot / self.figure_id
        d.mkdir(parents=True, exist_ok=True)
        return d

    def save(self) -> Path:
        payload = {
            "figure_id": self.figure_id,
            "generated": datetime.datetime.now().isoformat(timespec="seconds"),
            "meta": self.meta,
            "values": self.values,
        }
        path = self.figdir() / "stats.json"
        path.write_text(json.dumps(payload, indent=2))
        return path
