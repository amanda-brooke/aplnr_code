from __future__ import annotations

import pickle
from pathlib import Path
import pandas as pd

from tracking.data_loader import BaseDataLoader, discover_offsets, discover_group_folders
from tracking.radial_analyzer import RadialAnalyzer
from tracking.split_tracks import split_radial_tracks_on_division


def save_pickle(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump(df, f)
    print(f"Saved: {path}")


class TrackingPipeline:
    """
    Loads pre-processed tracking CSVs (spots/edges/tracks), combines across offsets,
    optionally cleans (splits divisions), optionally computes radial metrics, and saves
    dataset-prefixed pickles.
    """

    def __init__(self, base_data_path: Path, output_path: Path, dataset_id: str, layout: str):
        self.base_data_path = base_data_path
        self.output_path = output_path
        self.dataset_id = dataset_id

        layout = layout.lower()
        if layout not in ("flat", "tc"):
            raise ValueError("layout must be 'flat' or 'tc'")
        self.layout = layout

    def _load_flat(self, data_type: str, group_label: str) -> pd.DataFrame:
        """flat layout: offset_XX/<data_type>/*.csv"""
        dfs = []
        for offset_folder in discover_offsets(self.base_data_path):
            folder = offset_folder / data_type
            if not folder.exists():
                continue

            offset_str = offset_folder.name.split("_")[-1]
            time_offset = float(offset_str)

            loader = BaseDataLoader(
                folder=folder,
                group_label=group_label,
                data_type=data_type,
                time_offset=time_offset,
                layout="flat",
                dataset_id=self.dataset_id,
            )
            dfs.append(loader.load_all())

        return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()

    def _load_tc(self, data_type: str) -> pd.DataFrame:
        """tc layout: offset_XX/(control|treatment)/<data_type>/*.csv"""
        dfs = []
        for offset_folder in discover_offsets(self.base_data_path):
            groups = discover_group_folders(offset_folder)  # { "Control": {"spots": Path, ...}, ... }
            offset_str = offset_folder.name.split("_")[-1]
            time_offset = float(offset_str)

            for group_label, dt_dict in groups.items():
                folder = dt_dict.get(data_type, None)
                if folder is None or not folder.exists():
                    continue

                loader = BaseDataLoader(
                    folder=folder,
                    group_label=group_label,  # "Control" or "Treatment"
                    data_type=data_type,
                    time_offset=time_offset,
                    layout="tc",
                    dataset_id=self.dataset_id,
                )
                dfs.append(loader.load_all())

        return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()

    def load_combined(self, group_label_for_flat: str = "WT") -> dict[str, pd.DataFrame]:
        """
        Returns dict with keys spots/edges/tracks combined across offsets.
        For flat layout, everything gets Group=group_label_for_flat.
        For tc layout, Group comes from folder names (Control/Treatment).
        """
        all_data = {}
        for data_type in ["spots", "edges", "tracks"]:
            if self.layout == "flat":
                df = self._load_flat(data_type=data_type, group_label=group_label_for_flat)
            else:
                df = self._load_tc(data_type=data_type)
            all_data[data_type] = df
        return all_data

    def clean_edges(self, edges: pd.DataFrame) -> pd.DataFrame:
        """Split branched trajectories at divisions and add CLEAN_TRACK_ID."""
        if edges.empty:
            return edges
        return split_radial_tracks_on_division(edges)

    def radial_on_combined(self, spots: pd.DataFrame, edges: pd.DataFrame) -> pd.DataFrame:
        """Compute radial metrics grouped by File_ID."""
        if spots.empty or edges.empty:
            return pd.DataFrame()
        analyzer = RadialAnalyzer(spot_table=spots, edge_table=edges)
        processed_edges = analyzer.process_all_files()
        analyzer.calculate_radial_persistence(processed_edges)
        return processed_edges

    def run(self, group_label_for_flat: str = "WT", do_clean: bool = True, do_radial: bool = False) -> dict[str, pd.DataFrame]:
        """
        Run full pipeline for ONE dataset.
        Saves raw combined pickles; saves cleaned edges; optionally saves radial edges.
        """
        combined = self.load_combined(group_label_for_flat=group_label_for_flat)

        # Save raw combined
        save_pickle(combined["spots"],  self.output_path / f"{self.dataset_id}_spots_raw.pkl")
        save_pickle(combined["edges"],  self.output_path / f"{self.dataset_id}_edges_raw.pkl")
        save_pickle(combined["tracks"], self.output_path / f"{self.dataset_id}_tracks_raw.pkl")

        # Clean
        edges_for_analysis = combined["edges"]
        if do_clean:
            edges_for_analysis = self.clean_edges(edges_for_analysis)
            save_pickle(edges_for_analysis, self.output_path / f"{self.dataset_id}_edges_cleaned.pkl")
        else:
            save_pickle(edges_for_analysis, self.output_path / f"{self.dataset_id}_edges_uncleaned.pkl")

        # Optional radial (only meaningful for circles)
        if do_radial:
            radial_edges = self.radial_on_combined(combined["spots"], edges_for_analysis)
            save_pickle(radial_edges, self.output_path / f"{self.dataset_id}_radial_edges.pkl")
            combined["radial_edges"] = radial_edges

        combined["edges_for_analysis"] = edges_for_analysis
        return combined
