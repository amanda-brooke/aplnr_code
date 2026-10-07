"""
config.py -- the one place that says where the data lives.

All scripts and notebooks import ROOT from here instead of hardcoding a path.

ROOT is the project data folder, i.e. the folder that contains 00_DATA_FROZEN/
(raw/annotated h5ad objects, Stan's barcode list, trained models), 01_references/
(IH gastruloid timeseries, Tyser embryo) and where outputs are written (06_tables/,
07_figures/, pseudobulk_*/).

Set it with an environment variable:
    export SC_WORKFLOW_ROOT=/path/to/sc_workflow        (macOS/Linux)
    set SC_WORKFLOW_ROOT=V:\\path\\to\\sc_workflow        (Windows)

If unset, the folder that contains this repository is used (the original layout,
where the code folder sits inside the data folder).
"""
import os
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]

ROOT = Path(os.environ.get("SC_WORKFLOW_ROOT", _REPO.parent)).expanduser().resolve()
# Parent of ROOT: holds a few inputs that pre-date the sc_workflow folder
# (the Cell Ranger outputs and the IH reference scVI model used in preprocessing).
PROJECT_DIR = ROOT.parent

if not (ROOT / "00_DATA_FROZEN").is_dir():
    raise FileNotFoundError(
        f"No 00_DATA_FROZEN/ under {ROOT}. Set SC_WORKFLOW_ROOT to the data folder.")

# Raw Cell Ranger outputs (only needed to rerun preprocessing from scratch).
CELLRANGER_DIR = Path(os.environ.get(
    "CELLRANGER_DIR", "/mnt/e/scRNA Michigan and Chicago/cellranger_results"))
