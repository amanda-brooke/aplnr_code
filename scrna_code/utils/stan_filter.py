"""
Helper functions for scRNAseq data preprocessing.

Imported in:
- 02_analysis/kd_pseudobulk.py
- 02_analysis/kd_aplnr_correlation.py
- 02_analysis/run_pseudobulk_ALL.py
- 02_analysis/export_control_expression.py
- 04_supplement/qc_supplement.py
- 04_supplement/ih_validation.py
- 04_supplement/supp_figures.py

- sample_barcodes.csv includes cells passed qc_pipeline

USAGE
    from stan_filter import stan_mask
    keep = stan_mask(obs_names, ROOT)      # boolean array, same length as obs_names
    M, obs = M[keep], obs.loc[keep].reset_index(drop=True)
"""
import re
import numpy as np
import pandas as pd
from anndata import AnnData
from pathlib import Path

BARCODE_RE = re.compile(r"[ACGT]{14,}")
CSV_NAME = "sample_barcodes.csv"

# Canonical home is the frozen data folder. The project-root location is the
# pre-freeze path and is kept as a fallback so older checkouts still work.
CSV_CANDIDATES = (
    Path("00_DATA_FROZEN") / "01_source" / CSV_NAME,
    Path(CSV_NAME),
)


def parse_barcode(x):
    """First run of >=14 ACGT chars, or None. Matches FIGURE2_CLEAN's bc16()."""
    h = BARCODE_RE.findall(str(x))
    return h[0] if h else None


def stan_barcodes(root):
    p = next((Path(root) / c for c in CSV_CANDIDATES if (Path(root) / c).exists()), None)
    if p is None:
        tried = "\n  ".join(str(Path(root) / c) for c in CSV_CANDIDATES)
        raise FileNotFoundError(
            f"sample_barcodes.csv not found. Tried:\n  {tried}\n"
            f"Stan's barcode list is the single cell filter for this "
            f"project; without it the KD pipeline silently analyses 45,435 cells while "
            f"the figure panels use 29,419.")
    s = pd.read_csv(p, index_col=0)
    bcs = {parse_barcode(x) for x in s.index}
    bcs.discard(None)
    return bcs


def stan_mask(obs_names, root, verbose=True):
    """Boolean mask over obs_names for cells present in Stan's barcode list."""
    bcs = stan_barcodes(root)
    parsed = np.array([parse_barcode(x) for x in obs_names], dtype=object)
    keep = np.array([b is not None and b in bcs for b in parsed])
    n_unparsed = int(sum(b is None for b in parsed))
    if verbose:
        print(f"Stan filter: {keep.sum():,}/{len(keep):,} cells kept "
              f"({100*keep.mean():.1f}%)"
              + (f"   [{n_unparsed:,} obs_names had no parseable barcode]"
                 if n_unparsed else ""))
    if keep.sum() == 0:
        raise RuntimeError("Stan filter kept ZERO cells -- barcode formats do not match. "
                           "Check that obs_names still contain the raw 16bp barcode.")
    if keep.mean() < 0.4:
        print(f"  WARNING: retention {100*keep.mean():.1f}% is lower than the expected "
              f"~64.7%. Verify the barcode list is the current one.")
    return keep


def read_obs_names(h5_path):
    """obs_names from an .h5ad without loading the matrix."""
    import h5py
    with h5py.File(h5_path, "r") as h:
        og = h["obs"]
        ik = og.attrs.get("_index", "_index")
        ik = ik.decode() if isinstance(ik, bytes) else ik
        it = og[ik]
        raw = it["values"][:] if isinstance(it, h5py.Group) else it[:]
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in raw])


# Library size gates, raw counts, inclusive on both sides
LIBRARY_SIZE_GATE = {
    "Chicago_Ctrl_48h":  (15_000, 35_000),
    "Chicago_DOX_48h":   (10_000, 35_000),
    "Chicago_Ctrl_72h":  ( 4_500, 20_000),
    "Chicago_DOX_72h":   ( 5_000, 20_000),
    "Michigan_Ctrl_48h": ( 7_000, 25_000),
    "Michigan_DOX_48h":  ( 8_000, 25_000),
    "Michigan_Ctrl_72h": ( 5_000, 30_000),
    "Michigan_DOX_72h":  ( 4_000, 20_000),
}

def library_size_gate(sample_id: np.ndarray | list, total_counts: np.ndarray | list):
    """Library size filter (based on LIBRARY_SIZE_GATE)

    Parameters
    ----------
    sample_id : np.ndarray | list
        Sample labels (e.g. "Michigan_Ctrl_72h")
    total_counts : np.ndarray | list
        UMIs per cell

    Returns
    -------
    A boolean array to keep cells within the library size bounds
    """
    sample_id = np.asarray(sample_id, dtype=object)
    total_counts = np.asarray(total_counts, dtype=float)
    unknown = sorted(set(sample_id) - set(LIBRARY_SIZE_GATE))
    if unknown:
        raise KeyError(f"No library-size gate defined for {unknown}",
                       f"Valid sample ids: {list(LIBRARY_SIZE_GATE.keys())}")
    lo = np.array([LIBRARY_SIZE_GATE[s][0] for s in sample_id], float)
    hi = np.array([LIBRARY_SIZE_GATE[s][1] for s in sample_id], float)
    return (total_counts >= lo) & (total_counts <= hi)

def qc_pipeline(data_path: str | Path) -> AnnData:
    from pathlib import Path
    import scanpy as sc
    """FUll QC pipeline starting from the cellranger outputs.
    Author: Zhiyuan Yu
    ./sample_barcodes.csv includes cells thare are included from this pipeline

    Parameters
    ----------
    data_path: str | Path
        path to folder with cellranger output
        liang_amanda_data_raw
        ├── chicago
        │   ├── chicago-sample1-raw
        │   │   ├── barcodes.tsv.gz
        │   │   ├── features.tsv.gz
        │   │   └── matrix.mtx.gz
        │   ├── chicago-sample2-raw
        │   │   ├── ...
        │   ├── chicago-sample3-raw
        │   │   ├── ...
        │   └── chicago-sample4-raw
        │       └── ...
        └── michigan
            ├── michigan-sample1-raw
            │   ├── ...
            ├── michigan-sample2-raw
            │   ├── ...
            ├── michigan-sample3-raw
            │   ├── ...
            └── michigan-sample4-raw
                └── ...
    """

    """
    =============== data import ================
    - Load cellranger mtx output
    - Index each cell with label_barcode
        e.g. AAACCCTGTAATCGCT-1_Chicago_Ctrl_48h
    ============================================
    """
    data_path = Path(data_path)
    condition_map = {
        "chicago-sample1-raw": "Chicago_Ctrl_48h",
        "chicago-sample2-raw": "Chicago_DOX_48h",
        "chicago-sample3-raw": "Chicago_Ctrl_72h",
        "chicago-sample4-raw": "Chicago_DOX_72h",
        "michigan-sample1-raw": "Michigan_Ctrl_48h",
        "michigan-sample2-raw": "Michigan_DOX_48h",
        "michigan-sample3-raw": "Michigan_Ctrl_72h",
        "michigan-sample4-raw": "Michigan_DOX_72h",
    }
    data_list = {}
    for ss, ll in condition_map.items():
        sample_source = ss.split("-")[0]
        data_list[ll] = sc.read_10x_mtx(data_path / sample_source / ss)
        print(data_list[ll].shape)
        data_list[ll].obs.index = data_list[ll].obs.index.str.cat(
            [ll] * data_list[ll].shape[0], sep="_"
        )

    """
    ============ Library filter =============
    - Keep cells with total UMIs in [lo, hi]
    =========================================
    """
    def filter_adata_lib(adata: AnnData, thresholds: tuple):
        lib_size = np.array(adata.X.sum(axis=1)).squeeze()
        include_idx = np.logical_and(lib_size >= thresholds[0], lib_size <= thresholds[1])
        return adata[include_idx, :].copy()

    for ll, (lo, hi) in LIBRARY_SIZE_GATE.items():
        print(f"Before lib filt {ll}: {data_list[ll].shape}")
        data_list[ll] = filter_adata_lib(data_list[ll], (lo, hi))
        print(f"After lib filt {ll}: {data_list[ll].shape}")



    """
    ============ Filter lowly expressed genes =============
    - anndata to pandas dataframe as scprep expects that
    - Genes expressed in fewer than 10 cells are dropped
    =======================================================
    """
    import scprep
    def adata_to_df(adata):
        var_names = adata.var_names.str.cat([f" ({g})" for g in adata.var["gene_ids"]])
        df = adata.to_df()
        df.columns = var_names
        return df
    for ll, dat in data_list.items():
        dat = adata_to_df(dat)
        data_list[ll] = scprep.filter.filter_rare_genes(dat, min_cells=10)

    """
    ============== Filter cells with high MT ===============
    - Compute UMI fraction of mitochondrial genes per cell
    - Drop cells with MT fraction > 20%
    ========================================================
    """
    def get_total_mt_fraction(data_frame: pd.DataFrame, mt_list: list[str]):
        mt_list = [g for g in mt_list if g in data_frame.columns]
        mt_expr = data_frame.loc[:, mt_list].sum(axis=1)
        total_expr = data_frame.sum(axis=1)
        return mt_expr / total_expr

    # This explicitly matches mito genes start with MT- (e.g. MT-CO1)
    # Genes such as MT1A (Metallothionein 1A) wont be matched
    # Since gene filtering was done in the previous step,
    # here the intersection of MT genes across samples is used
    mt_genes = list(
        set.intersection(
            *[
                set(s)
                for s in [
                    [g for g in dat.columns if g.startswith("MT-")]
                    for dat in data_list.values()
                ]
            ]
        )
    )
    for ll, dat in data_list.items():
        mt_expr = get_total_mt_fraction(dat, mt_genes)
        data_list[ll] = dat[mt_expr <= 0.2]

    # Generate a single anndata
    adata_list = []
    for ll, dat in data_list.items():
        adata = AnnData(dat)
        adata.obs["sample_labels"] = ll
        adata_list.append(adata)
    return sc.concat(adata_list, join="outer")

def verify_gate_matches_list(sample_id, total_counts, obs_names, root):
    """Prove the gate and the supplied barcode list agree. Returns a dict."""
    gate = library_size_gate(sample_id, total_counts)
    listed = stan_mask(obs_names, root, verbose=False)
    return {"gate": int(gate.sum()), "barcode_list": int(listed.sum()),
            "both": int((gate & listed).sum()),
            "gate_only": int((gate & ~listed).sum()),
            "list_only": int((~gate & listed).sum())}
