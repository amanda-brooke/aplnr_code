#!/usr/bin/env python3
"""Cross-platform runner for the MPS/APS bulk RNA-seq pipeline.

Run from the project folder, for example:
  python scripts/run_pipeline.py --config config.json status
  python scripts/run_pipeline.py --config config.json verify
  python scripts/run_pipeline.py --config config.json qc
  python scripts/run_pipeline.py --config config.json quant
  python scripts/run_pipeline.py --config config.json deseq2
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

REQUIRED_MANIFEST_COLUMNS = {
    "sample_id", "celltype", "replicate",
    "condition", "day", "include", "r1", "r2"
}


def load_config(config_path: str | Path) -> tuple[dict[str, Any], Path]:
    config_path = Path(config_path).expanduser().resolve()
    with config_path.open(encoding="utf-8") as handle:
        cfg = json.load(handle)
    root_setting = Path(cfg.get("project_root", "."))
    project_root = root_setting if root_setting.is_absolute() else (config_path.parent / root_setting)
    return cfg, project_root.resolve()


def resolve_path(project_root: Path, value: str) -> Path:
    path = Path(os.path.expandvars(os.path.expanduser(value)))
    return path.resolve() if path.is_absolute() else (project_root / path).resolve()


def command_value(cfg: dict[str, Any], name: str) -> str:
    return str(cfg.get("commands", {}).get(name, name))


def command_exists(command: str) -> bool:
    p = Path(command)
    return p.exists() if (p.is_absolute() or p.parent != Path('.')) else shutil.which(command) is not None


def run_command(cmd: list[str], log_path: Path | None = None) -> None:
    printable = subprocess.list2cmdline(cmd) if os.name == "nt" else " ".join(repr(x) if " " in x else x for x in cmd)
    print(f"\nRUNNING: {printable}")
    if log_path is None:
        subprocess.run(cmd, check=True)
        return
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8", errors="replace") as log:
        log.write(f"COMMAND: {printable}\n\n")
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd, output=f"See {log_path}")


def read_manifest(cfg: dict[str, Any], root: Path, include_only: bool = True) -> list[dict[str, str]]:
    manifest_path = resolve_path(root, cfg["manifest"])
    with manifest_path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"Manifest is empty: {manifest_path}")
    missing = REQUIRED_MANIFEST_COLUMNS - set(rows[0])
    if missing:
        raise ValueError(f"Manifest is missing columns: {sorted(missing)}")
    if include_only:
        rows = [r for r in rows if r["include"].strip().lower() in {"true", "1", "yes", "y"}]
    return rows


def validate_manifest(rows: list[dict[str, str]]) -> None:
    sample_ids = [r["sample_id"].strip() for r in rows]
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("sample_id values are not unique")
    fastqs = [r[k].strip() for r in rows for k in ("r1", "r2")]
    if len(fastqs) != len(set(fastqs)):
        raise ValueError("A FASTQ filename appears more than once in the manifest")
    allowed = {
        "celltype": {"MPS", "APS"},
        "condition": {"CTRL", "DOX"},
        "day": {"D1", "D2"},
    }
    for key, values in allowed.items():
        observed = {r[key].strip() for r in rows}
        if not observed <= values:
            raise ValueError(f"Unexpected {key} values: {sorted(observed - values)}")
    # Confirm every replicate is a complete 2 x 2 x 2 block (celltype x condition x day).
    expected = {(ct, cond, day) for ct in allowed["celltype"] for cond in allowed["condition"] for day in allowed["day"]}
    for rep in sorted({r["replicate"] for r in rows}):
        rep_rows = [r for r in rows if r["replicate"] == rep]
        observed = {(r["celltype"], r["condition"], r["day"]) for r in rep_rows}
        if observed != expected:
            raise ValueError(f"Replicate {rep} is not a complete MPS/APS x CTRL/DOX x D1/D2 block")


def expected_paths(cfg: dict[str, Any], root: Path, rows: list[dict[str, str]]) -> list[tuple[str, str, Path]]:
    raw = resolve_path(root, cfg["raw_fastq_dir"])
    return [(r["sample_id"], read, raw / r[read]) for r in rows for read in ("r1", "r2")]


def build_tx2gene(fasta: Path, output: Path) -> None:
    print(f"Building transcript-to-gene map from {fasta}")
    opener = gzip.open if fasta.suffix == ".gz" else open
    output.parent.mkdir(parents=True, exist_ok=True)
    with opener(fasta, "rt", encoding="utf-8") as inp, gzip.open(output, "wt", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["target_id", "gene_id", "gene_name"])
        n = 0
        for line in inp:
            if not line.startswith(">"):
                continue
            target = line[1:].strip().split()[0]
            parts = target.split("|")
            if len(parts) < 6 or not parts[0].startswith("ENST") or not parts[1].startswith("ENSG"):
                raise ValueError(f"Unexpected GENCODE FASTA header: {target}")
            writer.writerow([target, parts[1], parts[5]])
            n += 1
    print(f"Wrote {n:,} transcript mappings to {output}")


def prepare(cfg: dict[str, Any], root: Path) -> None:
    for key in ["kallisto_output_dir", "fastqc_output_dir", "multiqc_output_dir", "deseq2_output_dir", "sensitivity_output_dir", "expression_output_dir", "plot_output_dir", "log_dir"]:
        resolve_path(root, cfg[key]).mkdir(parents=True, exist_ok=True)
    rows = read_manifest(cfg, root)
    validate_manifest(rows)
    tx2gene = resolve_path(root, cfg["tx2gene"])
    fasta = resolve_path(root, cfg["transcript_fasta"])
    if not tx2gene.exists():
        if not fasta.exists():
            raise FileNotFoundError(f"Neither tx2gene nor transcript FASTA exists. Missing: {tx2gene} and {fasta}")
        build_tx2gene(fasta, tx2gene)
    print(f"Manifest validated: {len(rows)} samples, {2 * len(rows)} FASTQ files expected.")
    print(f"Project root: {root}")
    print("Preparation complete.")


def status(cfg: dict[str, Any], root: Path) -> None:
    rows = read_manifest(cfg, root)
    validate_manifest(rows)
    entries = expected_paths(cfg, root, rows)
    present = [(s, r, p) for s, r, p in entries if p.exists()]
    missing = [(s, r, p) for s, r, p in entries if not p.exists()]
    total_bytes = sum(p.stat().st_size for _, _, p in present)
    print(f"FASTQs present: {len(present)}/{len(entries)} ({total_bytes / 1024**3:.1f} GiB downloaded)")
    if missing:
        print("\nMissing files:")
        for sample, read, path in missing:
            print(f"  {sample:20s} {read.upper()}: {path.name}")
    else:
        print("All expected FASTQ files are present.")


def parse_checksum_file(path: Path) -> dict[str, str]:
    checksums: dict[str, str] = {}
    with path.open(encoding="utf-8-sig") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            match = re.match(r"^([0-9a-fA-F]{64})\s+[* ]?(.+)$", line)
            if not match:
                continue
            checksums[Path(match.group(2).strip()).name] = match.group(1).lower()
    if not checksums:
        raise ValueError(f"No SHA256 entries could be parsed from {path}")
    return checksums


def sha256_file(path: Path, chunk_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify(cfg: dict[str, Any], root: Path) -> None:
    rows = read_manifest(cfg, root)
    validate_manifest(rows)
    missing = [p for _, _, p in expected_paths(cfg, root, rows) if not p.exists()]
    if missing:
        raise FileNotFoundError(f"Cannot verify checksums: {len(missing)} expected FASTQ files are still missing. Run 'status' first.")
    checksum_path = resolve_path(root, cfg["checksum_file"])
    if not checksum_path.exists():
        raise FileNotFoundError(f"Checksum file not found: {checksum_path}")
    expected = parse_checksum_file(checksum_path)
    out_path = resolve_path(root, cfg["log_dir"]) / "sha256_verification.tsv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    failures = []
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["filename", "expected_sha256", "observed_sha256", "status"])
        for i, (_, _, path) in enumerate(expected_paths(cfg, root, rows), start=1):
            if path.name not in expected:
                failures.append((path.name, "not listed in checksum file"))
                writer.writerow([path.name, "", "", "NOT_LISTED"])
                continue
            print(f"[{i:02d}/{2*len(rows):02d}] SHA256 {path.name}")
            observed = sha256_file(path)
            ok = observed == expected[path.name]
            writer.writerow([path.name, expected[path.name], observed, "OK" if ok else "MISMATCH"])
            if not ok:
                failures.append((path.name, "checksum mismatch"))
    if failures:
        raise RuntimeError(f"Checksum verification failed for {len(failures)} file(s). See {out_path}")
    print(f"All {2*len(rows)} FASTQ files passed SHA256 verification. Report: {out_path}")


def qc(cfg: dict[str, Any], root: Path, force: bool = False, available_only: bool = False) -> None:
    rows = read_manifest(cfg, root)
    all_entries = expected_paths(cfg, root, rows)
    missing = [p for _, _, p in all_entries if not p.exists()]
    if missing and not available_only:
        raise FileNotFoundError(
            f"FASTQC cannot run: {len(missing)} FASTQs are missing. "
            "Use --available-only to process completed downloads now."
        )
    entries = [(s, r, p) for s, r, p in all_entries if p.exists()]
    if not entries:
        print("No completed FASTQ files are available for FastQC yet.")
        return
    if missing:
        print(f"Partial QC mode: processing {len(entries)}/{len(all_entries)} available FASTQs; {len(missing)} are still missing.")

    fastqc_cmd = command_value(cfg, "fastqc")
    if not command_exists(fastqc_cmd):
        raise FileNotFoundError(f"FastQC command not found: {fastqc_cmd}. Set commands.fastqc in config.json.")
    out_dir = resolve_path(root, cfg["fastqc_output_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    def report_path(fastq: Path) -> Path:
        name = fastq.name
        for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
            if name.endswith(suffix):
                name = name[:-len(suffix)]
                break
        return out_dir / f"{name}_fastqc.zip"

    to_run = []
    for _, _, fastq in entries:
        report = report_path(fastq)
        if force or not report.exists():
            to_run.append(fastq)
    if to_run:
        cmd = [fastqc_cmd, "--threads", str(cfg.get("threads", 8)), "--outdir", str(out_dir)] + [str(p) for p in to_run]
        run_command(cmd, resolve_path(root, cfg["log_dir"]) / "fastqc.log")
        print(f"FastQC completed for {len(to_run)} newly processed file(s).")
    else:
        print("FastQC outputs already exist for every available FASTQ; skipping. Use --force to rerun.")
    run_multiqc(cfg, root)

def stranded_args(value: str) -> list[str]:
    value = value.strip().lower()
    if value in {"confirm_with_core", "required_confirmation", "confirm"}:
        raise ValueError("Library strandedness has not been confirmed. The form says directional but does not specify FR versus RF. Set strandedness in config.json before quantification.")
    if value == "unstranded":
        return []
    if value in {"fr", "fr-stranded", "forward"}:
        return ["--fr-stranded"]
    if value in {"rf", "rf-stranded", "reverse"}:
        return ["--rf-stranded"]
    raise ValueError("strandedness must be unstranded, fr-stranded, or rf-stranded")


def quant(cfg: dict[str, Any], root: Path, force: bool = False, available_only: bool = False) -> None:
    rows = read_manifest(cfg, root)
    validate_manifest(rows)
    kallisto_cmd = command_value(cfg, "kallisto")
    if not command_exists(kallisto_cmd):
        raise FileNotFoundError(f"kallisto command not found: {kallisto_cmd}. Set commands.kallisto in config.json.")
    index = resolve_path(root, cfg["kallisto_index"])
    if not index.exists():
        raise FileNotFoundError(f"Kallisto index not found: {index}")
    raw = resolve_path(root, cfg["raw_fastq_dir"])
    out_root = resolve_path(root, cfg["kallisto_output_dir"])
    log_root = resolve_path(root, cfg["log_dir"]) / "kallisto"
    out_root.mkdir(parents=True, exist_ok=True)
    extra = [str(x) for x in cfg.get("kallisto_extra_args", [])]
    strand = stranded_args(str(cfg.get("strandedness", "unstranded")))

    processed = 0
    already_complete = 0
    unavailable = 0
    for i, row in enumerate(rows, start=1):
        sample = row["sample_id"]
        r1, r2 = raw / row["r1"], raw / row["r2"]
        if not r1.exists() or not r2.exists():
            if available_only:
                missing_reads = []
                if not r1.exists():
                    missing_reads.append("R1")
                if not r2.exists():
                    missing_reads.append("R2")
                print(f"[{i:02d}/{len(rows):02d}] SKIP {sample}: missing {'/'.join(missing_reads)}")
                unavailable += 1
                continue
            raise FileNotFoundError(
                f"Missing FASTQ pair for {sample}: {r1}, {r2}. "
                "Use --available-only to quantify complete pairs that have already downloaded."
            )
        out_dir = out_root / sample
        complete = (out_dir / "abundance.tsv").exists() and (out_dir / "run_info.json").exists()
        if complete and not force:
            print(f"[{i:02d}/{len(rows):02d}] SKIP {sample}: complete output already exists")
            already_complete += 1
            continue
        if out_dir.exists():
            # Remove either a forced rerun or an incomplete output left by an interrupted/failed run.
            shutil.rmtree(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        cmd = [kallisto_cmd, "quant", "-i", str(index), "-o", str(out_dir), "-t", str(cfg.get("threads", 8))]
        cmd += strand + extra + [str(r1), str(r2)]
        print(f"[{i:02d}/{len(rows):02d}] {sample}")
        run_command(cmd, log_root / f"{sample}.log")
        processed += 1
    summarize_kallisto(cfg, root)
    print(
        f"Kallisto step finished: {processed} processed now, "
        f"{already_complete} already complete, {unavailable} waiting for FASTQ pairs."
    )

def summarize_kallisto(cfg: dict[str, Any], root: Path) -> None:
    rows = read_manifest(cfg, root)
    out_root = resolve_path(root, cfg["kallisto_output_dir"])
    output = out_root / "kallisto_run_summary.tsv"
    fields = ["sample_id", "n_processed", "n_pseudoaligned", "p_pseudoaligned", "n_unique", "p_unique", "n_targets", "kallisto_version", "index_version"]
    records = []
    for row in rows:
        info_path = out_root / row["sample_id"] / "run_info.json"
        if not info_path.exists():
            continue
        text = info_path.read_text(encoding="utf-8", errors="replace")
        try:
            info = json.loads(text)
        except json.JSONDecodeError:
            # Some Windows kallisto builds write unescaped backslashes in the
            # informational "call" field, making run_info.json invalid JSON.
            # The call field is not needed for the summary, so replace only that
            # line and parse the remaining metrics normally.
            cleaned_lines = []
            for line in text.splitlines():
                if re.match(r'^\s*"call"\s*:', line):
                    indent = line[: len(line) - len(line.lstrip())]
                    comma = "," if line.rstrip().endswith(",") else ""
                    cleaned_lines.append(f'{indent}"call": ""{comma}')
                else:
                    cleaned_lines.append(line)
            info = json.loads("\n".join(cleaned_lines))
        rec = {k: info.get(k, "") for k in fields}
        rec["sample_id"] = row["sample_id"]
        records.append(rec)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(records)
    print(f"Kallisto summary: {output} ({len(records)}/{len(rows)} samples)")


def run_multiqc(cfg: dict[str, Any], root: Path) -> None:
    multiqc_cmd = command_value(cfg, "multiqc")
    if not command_exists(multiqc_cmd):
        print(f"WARNING: MultiQC command not found ({multiqc_cmd}); skipping combined QC report.")
        return
    out_dir = resolve_path(root, cfg["multiqc_output_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    scan_dirs = [resolve_path(root, cfg["fastqc_output_dir"]), resolve_path(root, cfg["kallisto_output_dir"]), resolve_path(root, cfg["log_dir"])]
    cmd = [multiqc_cmd, "--force", "--outdir", str(out_dir)] + [str(x) for x in scan_dirs]
    run_command(cmd, resolve_path(root, cfg["log_dir"]) / "multiqc.log")


def deseq2(cfg: dict[str, Any], root: Path) -> None:
    rows = read_manifest(cfg, root)
    out_root = resolve_path(root, cfg["kallisto_output_dir"])
    missing = [r["sample_id"] for r in rows if not (out_root / r["sample_id"] / "abundance.tsv").exists()]
    if missing:
        raise FileNotFoundError(f"DESeq2 cannot run; missing kallisto abundance.tsv for {len(missing)} samples: {missing[:5]}")
    rscript = command_value(cfg, "rscript")
    if not command_exists(rscript):
        raise FileNotFoundError(f"Rscript command not found: {rscript}. Set commands.rscript in config.json.")
    script = root / "scripts" / "deseq2_analysis.R"
    config_path = root / "config.json"
    run_command([rscript, str(script), str(config_path)], resolve_path(root, cfg["log_dir"]) / "deseq2.log")


def tools(cfg: dict[str, Any], root: Path) -> None:
    for name in ["kallisto", "fastqc", "multiqc", "rscript"]:
        cmd = command_value(cfg, name)
        found = command_exists(cmd)
        print(f"{name:9s}: {'FOUND' if found else 'MISSING'}  ({cmd})")
    print(f"kallisto index: {'FOUND' if resolve_path(root, cfg['kallisto_index']).exists() else 'MISSING'}  ({resolve_path(root, cfg['kallisto_index'])})")
    print(f"tx2gene:       {'FOUND' if resolve_path(root, cfg['tx2gene']).exists() else 'MISSING'}  ({resolve_path(root, cfg['tx2gene'])})")


def main() -> None:
    parser = argparse.ArgumentParser(description="MPS/APS bulk RNA-seq pipeline")
    parser.add_argument("--config", default="config.json")
    parser.add_argument("command", choices=["prepare", "tools", "status", "verify", "qc", "quant", "summarize", "multiqc", "deseq2", "all"])
    parser.add_argument("--force", action="store_true", help="Rerun completed FastQC or kallisto steps")
    parser.add_argument("--available-only", action="store_true", help="For QC or quant, process only FASTQs/pairs that are currently present and skip missing downloads")
    args = parser.parse_args()
    cfg, root = load_config(args.config)
    if args.command == "prepare": prepare(cfg, root)
    elif args.command == "tools": tools(cfg, root)
    elif args.command == "status": status(cfg, root)
    elif args.command == "verify": verify(cfg, root)
    elif args.command == "qc": qc(cfg, root, force=args.force, available_only=args.available_only)
    elif args.command == "quant": quant(cfg, root, force=args.force, available_only=args.available_only)
    elif args.command == "summarize": summarize_kallisto(cfg, root)
    elif args.command == "multiqc": run_multiqc(cfg, root)
    elif args.command == "deseq2": deseq2(cfg, root)
    elif args.command == "all":
        prepare(cfg, root)
        verify(cfg, root)
        qc(cfg, root, force=args.force, available_only=False)
        quant(cfg, root, force=args.force, available_only=False)
        run_multiqc(cfg, root)
        deseq2(cfg, root)


if __name__ == "__main__":
    main()