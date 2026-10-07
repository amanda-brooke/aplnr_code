import os
import pandas as pd
from pathlib import Path


class BaseDataLoader:
    def __init__(
        self,
        folder: Path,
        group_label: str,
        data_type: str,
        time_offset: float = 0.0,
        layout: str = "flat",
        dataset_id: str | None = None
    ):
        self.folder = folder
        self.group_label = group_label
        self.data_type = data_type
        self.time_offset = time_offset
        self.dataset_id = dataset_id

        layout = layout.lower()
        if layout not in ("flat", "tc"):
            raise ValueError("layout must be 'flat' or 'tc'")
        self.layout = layout

    def load_and_prepare_data(self, file_path: Path) -> pd.DataFrame:
        print(f"Loading data from {file_path}")
        data_table = pd.read_csv(file_path, header=0, skiprows=[1, 2])
        data_table = data_table.apply(pd.to_numeric, errors='coerce', axis=0)
        data_table.fillna(0, inplace=True)

        # Get the file stem (e.g., "c1_spots" or "c1_edges")
        stem = file_path.stem
        # If the stem ends with _<data_type>, remove that part to get the common base id.
        suffix = f"_{self.data_type}"
        if stem.endswith(suffix):
            base_id = stem[:-len(suffix)]
        else:
            base_id = stem

        # Determine replicate id based on folder layout
        if self.layout == "flat":
            replicate_id = self.folder.parent.name  # offset_27
        else:  # "tc"
            replicate_id = f"{self.folder.parent.parent.name}_{self.folder.parent.name}"  # offset_27_control

        prefix = f"{self.dataset_id}_" if self.dataset_id else ""
        data_table["File_ID"] = f"{prefix}{base_id}_{replicate_id}"
        data_table["Group"] = self.group_label


        if self.data_type == 'edges':
            data_table['EDGE_TIME'] = (data_table['EDGE_TIME'] / 3600) + self.time_offset
            data_table['SPEED'] = data_table['SPEED'] * 60
        elif self.data_type == 'tracks':
            data_table['TRACK_START_HOURS'] = (data_table['TRACK_START'] / 3600) + self.time_offset
            data_table['TRACK_STOP_HOURS'] = (data_table['TRACK_STOP'] / 3600) + self.time_offset
            data_table['TRACK_DURATION_HOURS'] = data_table['TRACK_DURATION'] / 3600
        elif self.data_type == 'spots':
            data_table['SPOT_TIME'] = (data_table['POSITION_T'] / 3600) + self.time_offset

        return data_table


    def load_all(self) -> pd.DataFrame:
        data_tables = []
        for file in os.listdir(self.folder):
            if file.endswith('.csv'):
                file_path = self.folder / file
                dt = self.load_and_prepare_data(file_path)
                data_tables.append(dt)
        if data_tables:
            self.combined_data = pd.concat(data_tables, ignore_index=True)
        else:
            self.combined_data = pd.DataFrame()
        return self.combined_data


def discover_group_folders(offset_folder: Path) -> dict:
    """
    Given a root offset folder, discover group subfolders.
    Returns a dictionary with keys like 'Treatment', 'Control', etc., mapping to
    a dictionary of data types and their respective folder paths.
    For example:
    {
        'Treatment': {'spots': Path(...), 'edges': Path(...), 'tracks': Path(...)},
        'Control': {'spots': Path(...), 'edges': Path(...), 'tracks': Path(...)}
    }
    """
    groups = {}
    # Assume groups are direct subfolders (e.g., treatment, control)
    for group_dir in offset_folder.iterdir():
        if group_dir.is_dir():
            group_name = group_dir.name.capitalize()  # e.g., Treatment, Control
            groups[group_name] = {}
            # For each data type, check if a subfolder exists
            for data_type in ['spots', 'edges', 'tracks']:
                dt_folder = group_dir / data_type
                if dt_folder.exists():
                    groups[group_name][data_type] = dt_folder
    return groups



def discover_offsets(base_data_path: Path) -> list[Path]:
    return sorted(base_data_path.glob("offset_*"))
