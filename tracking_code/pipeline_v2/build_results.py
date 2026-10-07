"""
build_results.py — compute escape for every dataset in analysis.yml and write ONE tidy table.

Run once (offline; the big pickles are slow). Every figure then just loads
data/processed/escape_results.(parquet|pkl) and filters by `dataset` — no recomputation,
consistent numbers across all figures.

    python pipeline_v2/build_results.py                 # all datasets in analysis.yml
    python pipeline_v2/build_results.py shaplnr_mps      # just some
"""
from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
import dataconfig as dc
import ps_escape as pe

OUT = "data/processed/escape_results"   # extension chosen at write time


def build(names=None, data_dir="data", out=OUT):
    cfgs = dc.load_configs()
    names = names or list(cfgs)
    frames = []
    for name in names:
        try:
            df, cfg = dc.load_dataset(name, data_dir=data_dir)
            coord, edges = dc.coord_and_edges(cfg)
            res, info = pe.run_escape(df, coord, edges,
                                      tuple(cfg["reference_window"]),
                                      cfg.get("sigma_k", 1.5),
                                      group_col=cfg.get("group_col"))
            res["dataset"]  = name
            res["units"]    = cfg["units"]
            res["mu"]       = info["mu"]
            res["sigma"]    = info["sigma"]
            res["boundary"] = info["boundary"]
            frames.append(res)
            print(f"  {name:28s} files={res['File_ID'].nunique():>2}  boundary={info['boundary']:.3g}  rows={len(res)}")
        except Exception as e:
            print(f"  {name:28s} SKIP ({type(e).__name__}: {e})")
    if not frames:
        print("nothing built"); return None
    allres = pd.concat(frames, ignore_index=True)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    try:
        path = out + ".parquet"; allres.to_parquet(path)
    except Exception:
        path = out + ".pkl"; allres.to_pickle(path)
    print(f"WROTE {path}  shape={allres.shape}  datasets={sorted(allres['dataset'].unique())}")
    return allres


if __name__ == "__main__":
    build(sys.argv[1:] or None)
