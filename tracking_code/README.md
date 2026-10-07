# Gastruloid T⁺ cell migration — data pipeline

Reproducibility code for the analysis **pipeline** behind the paper: raw TrackMate
exports are processed into a single consolidated "master" table (radial
coordinates, migration metrics, PS-like boundary). Figure notebooks are provided
separately; this repository is the upstream data pipeline they read from.

## Pipeline order

1. **Ingest raw tracks → radial-annotated pickles**

   ```
   python run.py --dataset <name> --out ./data
   ```

   Datasets are defined in `datasets.yml`. For each dataset this loads the
   TrackMate spots/edges/tracks, cleans them (splits divisions), and computes the
   radial coordinate system: the colony **center** (boundary-circle fit),
   `radial_distance`, `dist_from_edge`, `angle_to_center`, and per-step
   velocities. Output: `data/<dataset>_radial_edges.pkl`.

2. **Build the consolidated master table**

   ```
   python pipeline_v2/build_master.py
   ```

   Reads every `data/*_radial_edges.pkl` according to `pipeline_v2/analysis.yml`
   (per-dataset parameters and exclusions — the single source of truth), applies
   the PS-like boundary fit (Gaussian core + inward tail, μ + 1.5σ), and writes a
   partitioned master to `data/processed/master/` plus per-track summaries and the
   boundary table. Schema is documented in `pipeline_v2/MASTER_SCHEMA.md`.

3. **(Optional) Precompute escape results**

   ```
   python pipeline_v2/build_results.py
   ```

   Derives the escape/dispersal result tables from the master.

Figures are generated from the master table in separate notebooks (not included
in this pipeline repository).

## Layout

| Path | Role |
|------|------|
| `run.py`, `datasets.yml` | raw-track ingestion entry point + dataset manifest |
| `tracking/` | upstream pipeline: loading (`data_loader.py`), division splitting (`split_tracks.py`), radial/center computation (`radial_analyzer.py`), per-track migration summaries (`prep.py`), filtering/stats helpers |
| `pipeline_v2/build_master.py` | assembles the consolidated master table |
| `pipeline_v2/analysis.yml` | per-dataset parameters + exclusions (single source) |
| `pipeline_v2/ps_escape.py` | PS-like core fitting + escape metrics (canonical μ + 1.5σ) |
| `pipeline_v2/ps_profile.py`, `ps_stats.py` | radial-profile fitting + statistics helpers |
| `pipeline_v2/dataconfig.py` | master/table loaders (`load_master`, `load_tracks`, `load_boundaries`) |
| `pipeline_v2/build_results.py` | optional escape-result precomputation |
| `pipeline_v2/MASTER_SCHEMA.md` | master table schema/contract |
| `pipeline_v2/test_escape_regression.py` | regression test (escape numbers) |

## Environment

```
pip install -r requirements.txt
```

## Data availability

Raw TrackMate exports and the built master table are not included in this
repository — see the Data Availability statement in the paper for the deposit.
