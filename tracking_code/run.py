"""
Single entry point for the TrackMate analysis pipeline.

Reads a dataset definition from a manifest (datasets.yml) and runs the full
pipeline for that dataset: load -> combine across offsets -> clean (split
divisions) -> optional radial metrics -> save dataset-prefixed pickles.

Run from the repository root, e.g.:

    python run.py --dataset sox2_live --out ./pickles
    python run.py --manifest datasets.yml --dataset shaplnr_mps --out ./pickles
"""
import argparse
from pathlib import Path

import yaml

from tracking.pipeline import TrackingPipeline


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest", default="datasets.yml",
                        help="Path to the dataset manifest YAML (default: datasets.yml)")
    parser.add_argument("--dataset", required=True,
                        help="Dataset key in the manifest")
    parser.add_argument("--out", required=True,
                        help="Output directory for pickles")
    args = parser.parse_args()

    with open(args.manifest, "r") as f:
        cfg = yaml.safe_load(f)

    datasets = cfg["datasets"]
    if args.dataset not in datasets:
        raise ValueError(
            f"Dataset '{args.dataset}' not found in {args.manifest}. "
            f"Available: {sorted(datasets)}"
        )

    ds = datasets[args.dataset]

    output_dir = Path(args.out)
    output_dir.mkdir(parents=True, exist_ok=True)

    pipeline = TrackingPipeline(
        base_data_path=Path(ds["data_root"]),
        output_path=output_dir,
        dataset_id=args.dataset,
        layout=ds["layout"],
    )

    pipeline.run(
        group_label_for_flat=ds.get("group_label", "WT"),  # only used for flat layout
        do_clean=ds.get("clean", True),                    # optional; defaults to True
        do_radial=ds.get("radial", False),
    )


if __name__ == "__main__":
    main()
