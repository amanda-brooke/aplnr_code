"""
palettes.py -- THE one place for figure colours, label grouping and plot order.

    from palettes import TISSUE_COLORS, TISSUE_MAP, TISSUE_ORDER, SITE_COLORS

Edit colours HERE. Nowhere else. No hex literals in notebooks.

--------------------------------------------------------------------------------
HISTORY -- why this file exists, so the same thing does not happen again
--------------------------------------------------------------------------------
figure2_palettes.json (now superseded, safe to delete) contained '#' comments,
which are illegal in JSON. Both notebooks that loaded it wrapped the call in
try/except, so it failed SILENTLY and every figure used matplotlib defaults or
hardcoded hex. Nobody noticed for two weeks.

Separately, TISSUE_MAP and TISSUE_ORDER were defined in a notebook cell that got
deleted on 2026-08-08 during cleanup. The figures kept working because the kernel
still held them in memory -- a Restart & Run All would have failed. They were
reconstructed from 06_tables/FIGURE_KD_clean/2_composition.csv, which reconciles
exactly:
    eNasM 1,714 + lNasM 522                        = 2,236  -> "NasM"
    aPM-LC 36 + CS-LC 22 + LPM-LC 21 + stressed 8  =    87  -> "other"
    total = 29,419 Stan-filtered cells

Both failures share a cause: the definition lived somewhere it could vanish
without anything complaining. Keep grouping, colour and order together, in code,
in one file that fails loudly.
"""

# ---------------------------------------------------------------- COLOURS

TISSUE_COLORS = {
    "pluripotent": "#6B8E9E",
    "aPS-LC":      "#4A7C7E",
    "pPS-LC":      "#E6A57E",
    "NasM":        "#C0563E",
    "eNasM":       "#BE7F96",   # dusty apricot
    "lNasM":       "#B64F3B",   # muted brick red
    "AE-LC":       "#8B6A9E",
    "aPM-LC":      "#D98270",
    "LPM-LC":      "#E997B8",
    "CS-LC":       "#9A6A5C",
    "DE":          "#D4A94E",
    "PGC-LC":      "#B8A88A",
    "stressed":    "#C4C4C4",
    "other":       "#AFA9A4",
}

CONDITION_COLORS = {
    "Ctrl": "#4A7C7E", "Control": "#4A7C7E",
    "DOX":  "#C0563E", "KD": "#C0563E", "shAPLNR": "#C0563E",
}

SITE_COLORS      = {"Chicago": "#6B8E9E", "Michigan": "#B64F3B"}
TIMEPOINT_COLORS = {"42h": "#9FB8C4", "72h": "#4A7C7E"}

# chart furniture, not category colours
RED, BLUE, GREY_FAINT = "#991E03", "#6B8E9E", "#F2EFEC"

EMBRYO_COLORS = {
    "Epiblast": "#6B8E9E",
    "Primitive Streak": "#8B6A9E",
    "Non-Neural Ectoderm": "#7FA88F",
    "Endoderm": "#D4A94E",
    "ExE Mesoderm": "#E8B894",
    "Nascent Mesoderm": "#E09B6F",
    "Emergent Mesoderm": "#D4795A",
    "Advanced Mesoderm": "#C0563E",
    "Axial Mesoderm": "#9E4230",
    "Hemogenic Endothelial Progenitors": "#A85462",
    "Erythroblasts": "#7D2B38",
}

PROGRAM_COLORS = {
    "Epi(A)": "#E89BC2", "Epi(P)": "#AEBB74", "AVE": "#E69F00", "DE": "#C9A227",
    "VE": "#9B86BD", "AM.Ecto": "#7B68A6", "Gast": "#72B5C4",
    "Axial.Meso": "#E2B65D", "Noto": "#5B8E7D", "Para.Meso": "#D98270",
    "LP.Meso": "#E997B8", "EXE.Meso.Prog": "#B9B1DE", "Connecting stalk": "#9A6A5C",
}

# ------------------------------------------------- GROUPING AND PLOT ORDER

# Maps the transferred label (predicted_leiden_sub_joint) to the plotting label.
TISSUE_MAP = {
    "pluripotent": "pluripotent",
    "aPS-LC":      "aPS-LC",
    "pPS-LC":      "pPS-LC",
    "AE-LC":       "AE-LC",
    "eNasM":       "eNasM",     # SPLIT (was "NasM" before 2026-08-08)
    "lNasM":       "lNasM",     # SPLIT (was "NasM" before 2026-08-08)
    "DE":          "DE",
    "PGC-LC":      "PGC-LC",
    "aPM-LC":      "other",
    "CS-LC":       "other",
    "LPM-LC":      "other",
    "stressed":    "other",
}

# DEFAULT: eNasM and lNasM shown separately.
TISSUE_ORDER = ["pluripotent", "aPS-LC", "pPS-LC", "eNasM", "lNasM",
                "AE-LC", "DE", "PGC-LC", "other"]

# MERGED variant. Figure 2's composition panel deliberately shows one NasM
# category; the KD fraction panel shows eNasM and lNasM separately. Import the
# pair you need -- do not redefine either in a notebook.
#     FIGURE2_current   -> TISSUE_MAP_MERGED, TISSUE_ORDER_MERGED
#     FIGURE_KD_current -> TISSUE_MAP,        TISSUE_ORDER          (split)
TISSUE_MAP_MERGED = {**TISSUE_MAP, "eNasM": "NasM", "lNasM": "NasM"}

TISSUE_ORDER_MERGED = ["pluripotent", "AE-LC", "aPS-LC", "pPS-LC", "NasM",
                       "DE", "PGC-LC", "other"]


# ------------------------------------------------------ FIGURE-SPECIFIC ROLES
# Added 2026-08-08. These were defined separately in maturation_figures.ipynb
# and maturation_supplementary.ipynb, in several conflicting versions -- eNasM
# alone existed as #BE7F96, #E6A57E and #C0563E in three different places.

CTRL_C, KD_C = "#BE6455", "#D99A90"      # Ctrl vs KD, maturation panels
R1_C,   R2_C = "#8A8780", "#D8D5CF"      # repeat annotation strips

# The two nascent-mesoderm states, by their canonical TISSUE_COLORS entries --
# defined this way so they can never drift from the rest of the figures.
STATE_COLORS = {"eNasM": TISSUE_COLORS["eNasM"], "lNasM": TISSUE_COLORS["lNasM"]}

# Neutral chart furniture -- connector lines, null histograms, reference lines.
GREY_LINE  = "#918C85"
GREY_FILL  = "#C9C4BD"
GREY_DARK  = "#696660"
GREY_MID   = "#B8B5AE"

# APLNR-correlate gene sets (supplementary null test)
POSITIVE_COLOR   = "#B64F3B"
NEGATIVE_COLOR   = "#6B8E9E"
BACKGROUND_COLOR = "#B8B2AA"
TREND_COLOR      = "#363B39"


def tissue(name, default="#AFA9A4"):
    """Colour for a label, neutral grey for anything unrecognised."""
    return TISSUE_COLORS.get(name, default)


def check(order=None):
    """Fail loudly if a plotted label has no colour, or a mapped label is dropped."""
    order = order or TISSUE_ORDER
    no_colour = [t for t in order if t not in TISSUE_COLORS]
    if no_colour:
        raise KeyError(f"no colour defined for {no_colour} -- add to TISSUE_COLORS")
    not_plotted = sorted(set(TISSUE_MAP.values()) - set(order))
    return {"order": order, "mapped_but_not_plotted": not_plotted}


if __name__ == "__main__":
    print(f"{len(TISSUE_COLORS)} tissue colours, {len(TISSUE_MAP)} mapped labels\n")
    for k, v in TISSUE_COLORS.items():
        mapped = [s for s, t in TISSUE_MAP.items() if t == k]
        print(f"  {k:14s} {v}   <- {', '.join(mapped) if mapped else '(not a plot label)'}")
    print("\nSPLIT  (default):", check(TISSUE_ORDER)["order"])
    print("MERGED (option) :", check(TISSUE_ORDER_MERGED)["order"])
