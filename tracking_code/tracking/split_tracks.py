import pandas as pd
import networkx as nx

def split_radial_tracks_on_division(radial_edges: pd.DataFrame) -> pd.DataFrame:
    required = {"File_ID", "SPOT_SOURCE_ID", "SPOT_TARGET_ID"}
    missing = required - set(radial_edges.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")

    clean_edges_list = []
    track_counter = -1

    for _, file_edges in radial_edges.groupby("File_ID"):
        G = nx.DiGraph()
        G.add_edges_from(zip(file_edges["SPOT_SOURCE_ID"], file_edges["SPOT_TARGET_ID"]))

        all_sources = set(file_edges["SPOT_SOURCE_ID"])
        all_targets = set(file_edges["SPOT_TARGET_ID"])
        root_nodes = list(all_sources - all_targets)

        file_track_records = []

        for root in root_nodes:
            track_counter += 1
            stack = [(root, track_counter)]

            while stack:
                current_node, current_track_id = stack.pop()
                children = list(G.successors(current_node))

                if len(children) <= 1:
                    for child in children:
                        match = file_edges[
                            (file_edges["SPOT_SOURCE_ID"] == current_node) &
                            (file_edges["SPOT_TARGET_ID"] == child)
                        ].copy()
                        match["CLEAN_TRACK_ID"] = current_track_id
                        file_track_records.append(match)
                        stack.append((child, current_track_id))
                else:
                    for child in children:
                        track_counter += 1
                        new_id = track_counter
                        match = file_edges[
                            (file_edges["SPOT_SOURCE_ID"] == current_node) &
                            (file_edges["SPOT_TARGET_ID"] == child)
                        ].copy()
                        match["CLEAN_TRACK_ID"] = new_id
                        file_track_records.append(match)
                        stack.append((child, new_id))

        if file_track_records:
            clean_edges_list.append(pd.concat(file_track_records, ignore_index=True))

    return pd.concat(clean_edges_list, ignore_index=True) if clean_edges_list else pd.DataFrame()

if __name__ == "__main__":
    raise SystemExit("This module provides functions only. Run from a notebook/script.")
