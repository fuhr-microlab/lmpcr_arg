#!/usr/bin/env python3
"""
End-to-end Y-adapter ARG nanopore analysis pipeline.

Code autoannotated using Composer v2.5.

Orchestrates Steps 1–3 by calling the repository Python scripts with the same
logic and parameter choices used in the publication analysis:
  1. Pre-demux concatemer split + barcode demultiplexing
  2. Per-barcode primer/Y-adapter split, multi-ARG filter, target mapping,
     and start/end concatenation
  3. Clustering + RACON polishing of combined FASTAs

Author: O'Brien, M.E., Fuhrmeister, E.R., Marchand, J.A.
Date: 2026
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

from lib.default_params import EXAMPLE_PATHS, CUTADAPT_DEMUX

REPO_ROOT = Path(__file__).resolve().parent

# Publication analysis defaults (match HPC run parameters)
DEFAULT_BARCODE_WHITELIST = [
    "barcode01", "barcode02", "barcode03",
    "barcode04", "barcode05", "barcode06",
]
DEFAULT_TARGETS = ["ctx-m", "KPC", "OXA", "qnrS"]
DEFAULT_MIN_CLUSTER_SIZE = 15


def run(cmd, cwd=None):
    """Run a subprocess; raise on failure with captured output."""
    print("\n>>>", " ".join(str(c) for c in cmd), flush=True)
    try:
        completed = subprocess.run(
            [str(c) for c in cmd],
            check=True,
            capture_output=True,
            text=True,
            cwd=cwd,
        )
    except subprocess.CalledProcessError as e:
        print("=== Command failed ===", file=sys.stderr)
        print("Command:", " ".join(str(c) for c in cmd), file=sys.stderr)
        if e.stdout:
            print("--- STDOUT ---", file=sys.stderr)
            print(e.stdout, file=sys.stderr)
        if e.stderr:
            print("--- STDERR ---", file=sys.stderr)
            print(e.stderr, file=sys.stderr)
        raise
    if completed.stdout:
        print(completed.stdout, end="" if completed.stdout.endswith("\n") else "\n")
    if completed.stderr:
        print(completed.stderr, end="" if completed.stderr.endswith("\n") else "\n", file=sys.stderr)
    return completed


def should_process_barcode(name, whitelist):
    """True if name matches whitelist prefix (empty whitelist → all)."""
    if not whitelist:
        return True
    return any(name.startswith(w) for w in whitelist)


def list_demux_inputs(demux_dir):
    """Return demux FASTA/FASTQ paths, excluding unsorted."""
    files = []
    for pattern in ("*.fasta", "*.fastq"):
        files.extend(sorted(demux_dir.glob(pattern)))
    return [p for p in files if not p.name.startswith("unsorted.")]


def step1_split_and_demux(args, paths):
    """Step 1: pre-demux split + cutadapt demultiplexing."""
    print("=" * 70)
    print("STEP 1: Split concatemers and demultiplex")
    print("=" * 70)

    splittr_dir = paths["splittr_dir"]
    demux_dir = paths["demux_dir"]
    blast_db_dir = paths["blast_db_dir"]
    splittr_dir.mkdir(parents=True, exist_ok=True)
    demux_dir.mkdir(parents=True, exist_ok=True)
    blast_db_dir.mkdir(parents=True, exist_ok=True)

    splitter_out = splittr_dir / "split_reads.fasta"
    splitter_script = REPO_ROOT / "adapter_splitter.py"

    cmd = [
        sys.executable, splitter_script,
        args.input,
        args.barcodes,
        splitter_out,
        "--preset", "pre_demux",
        "--blast-db-dir", blast_db_dir,
    ]
    run(cmd)

    demux_script = REPO_ROOT / "cutadapt_demux.py"
    cmd = [
        sys.executable, demux_script,
        "--input", splitter_out,
        "--barcodes", args.barcodes,
        "--output-dir", demux_dir,
        "--cutadapt-bin", args.cutadapt_bin,
        "--error-rate", CUTADAPT_DEMUX["error_rate"],
        "--min-overlap", CUTADAPT_DEMUX["min_overlap"],
        "--min-len", CUTADAPT_DEMUX["min_len"],
        "--max-len", CUTADAPT_DEMUX["max_len"],
    ]
    if CUTADAPT_DEMUX["no_trim"]:
        cmd.append("--no-trim")
    if CUTADAPT_DEMUX["internal"]:
        cmd.append("--internal")
    run(cmd)

    print("Step 1 complete.")
    return demux_dir


def step2_per_barcode(args, paths, demux_dir):
    """Step 2: post-demux split, ARG filter, map, concatenate."""
    print("=" * 70)
    print("STEP 2: Per-barcode split, ARG filter, map, concatenate")
    print("=" * 70)

    split_dir = paths["split_dir"]
    arg_filter_dir = paths["arg_filter_dir"]
    arg_mapping_dir = paths["arg_mapping_dir"]
    arg_concat_dir = paths["arg_concat_dir"]
    split_dir.mkdir(parents=True, exist_ok=True)
    arg_filter_dir.mkdir(parents=True, exist_ok=True)
    arg_mapping_dir.mkdir(parents=True, exist_ok=True)
    arg_concat_dir.mkdir(parents=True, exist_ok=True)

    splitter_script = REPO_ROOT / "adapter_splitter.py"
    filter_script = REPO_ROOT / "adapter_filter.py"
    map_script = REPO_ROOT / "arg_filter_and_plot.py"

    demux_files = list_demux_inputs(demux_dir)
    if not demux_files:
        raise FileNotFoundError(f"No demux FASTA/FASTQ files found in {demux_dir}")

    # 2A — primer / Y-adapter splitting
    print("\n--- Step 2A: Primer / Y-adapter splitting ---")
    for demux_file in demux_files:
        barcode_name = demux_file.stem
        if not should_process_barcode(barcode_name, args.barcodes_whitelist):
            print(f"Skipping {demux_file.name} (not in barcode whitelist)")
            continue
        out_file = split_dir / f"{barcode_name}_split.fasta"
        run([
            sys.executable, splitter_script,
            demux_file,
            args.primers_and_y,
            out_file,
            "--preset", "post_demux",
        ])

    # 2B — multi-ARG hit filter
    print("\n--- Step 2B: Multi-ARG hit filtering ---")
    for split_file in sorted(split_dir.glob("*_split.fasta")):
        barcode_name = split_file.name[: -len("_split.fasta")]
        if not should_process_barcode(barcode_name, args.barcodes_whitelist):
            continue
        out_file = arg_filter_dir / f"{barcode_name}_filtered.fasta"
        run([
            sys.executable, filter_script,
            split_file,
            args.arg_sequences,
            out_file,
            "--max-hits", "1",
            "--min-length", "0",
            "--min-primer-overlap", "0.2",
            "--identity", "80",
            "--verbose",
        ])

    # 2C — ARG mapping per target
    print("\n--- Step 2C: ARG mapping / position filtering ---")
    if not map_script.exists():
        raise FileNotFoundError(
            f"Required script not found: {map_script}\n"
            "Add arg_filter_and_plot.py to the repository root to run Step 2C."
        )

    for target in args.targets:
        for filtered_file in sorted(arg_filter_dir.glob("*_filtered.fasta")):
            barcode_name = filtered_file.name[: -len("_filtered.fasta")]
            if not should_process_barcode(barcode_name, args.barcodes_whitelist):
                continue
            if filtered_file.stat().st_size == 0:
                print(f"WARNING: empty file {filtered_file}, skipping")
                continue
            barcode_out = arg_mapping_dir / target / barcode_name
            barcode_out.mkdir(parents=True, exist_ok=True)
            run([
                sys.executable, map_script,
                filtered_file,
                args.arg_sequences,
                barcode_out,
                "--verbose",
                "--debug",
                "--target", target,
            ])

    # 2D — concatenate start + end
    print("\n--- Step 2D: Concatenate start + end fragments ---")
    for target in args.targets:
        target_map = arg_mapping_dir / target
        if not target_map.is_dir():
            continue
        target_concat = arg_concat_dir / target
        target_concat.mkdir(parents=True, exist_ok=True)
        for barcode_dir in sorted(p for p in target_map.iterdir() if p.is_dir()):
            barcode_name = barcode_dir.name
            if not should_process_barcode(barcode_name, args.barcodes_whitelist):
                continue
            start_file = barcode_dir / f"start_{target}.fasta"
            end_file = barcode_dir / f"end_{target}.fasta"
            combined = target_concat / f"{barcode_name}_combined.fasta"
            if start_file.is_file() and end_file.is_file():
                combined.write_text(start_file.read_text() + end_file.read_text())
                print(f"Created {combined}")
            else:
                print(f"Skipping {barcode_name}/{target} (missing start or end file)")

    print("Step 2 complete.")
    return arg_concat_dir


def step3_cluster_polish(args, arg_concat_dir):
    """Step 3: cluster + polish each combined FASTA."""
    print("=" * 70)
    print("STEP 3: Cluster and polish")
    print("=" * 70)

    polish_script = REPO_ROOT / "cluster_polish.py"
    if not polish_script.exists():
        raise FileNotFoundError(f"Required script not found: {polish_script}")

    # Publication clustering searched combined FASTAs under arg_concat/<target>/
    combined_files = sorted(arg_concat_dir.glob("*/*_combined.fasta"))
    if not combined_files:
        print(f"WARNING: no *_combined.fasta files under {arg_concat_dir}")
        return

    for fasta_file in combined_files:
        filename = fasta_file.stem  # e.g. barcode01_FWD_combined
        match = re.search(r"barcode(\d+)", filename, re.IGNORECASE)
        if not match:
            print(f"ERROR: could not extract barcode number from {filename}; skipping")
            continue
        barcode_num = int(match.group(1))
        cluster_dirname = f"clusters{barcode_num:03d}"
        output_dir = fasta_file.parent / cluster_dirname

        if args.skip_existing_clusters and output_dir.is_dir():
            print(f"Skipping {filename} ({cluster_dirname} already exists)")
            continue

        output_dir.mkdir(parents=True, exist_ok=True)
        run([
            sys.executable, polish_script,
            fasta_file,
            output_dir,
            "--min-cluster-size", str(args.min_cluster_size),
            "--polish",
        ])

    print("Step 3 complete.")


def build_paths(output_dir: Path):
    return {
        "output_dir": output_dir,
        "splittr_dir": output_dir / "splittr",
        "demux_dir": output_dir / "demux",
        "split_dir": output_dir / "split",
        "arg_filter_dir": output_dir / "arg_filter",
        "arg_mapping_dir": output_dir / "arg_mapping",
        "arg_concat_dir": output_dir / "arg_concat",
        "blast_db_dir": output_dir / "tmp_blast_db",
    }


def parse_args():
    p = argparse.ArgumentParser(
        description="Run the full Y-adapter ARG nanopore analysis pipeline.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--input", required=True,
        help="FASTQ from nanopore basecalling (see publication for SRA accessions)",
    )
    p.add_argument(
        "--barcodes", default=EXAMPLE_PATHS["barcodes"],
        help="Barcode FASTA (FWD headers)",
    )
    p.add_argument(
        "--primers-and-y", default=EXAMPLE_PATHS["primers_and_y"],
        help="Primer + Y-adapter FASTA",
    )
    p.add_argument(
        "--arg-sequences", default=EXAMPLE_PATHS["arg_sequences"],
        help="ARG sequence FASTA",
    )
    p.add_argument(
        "--output-dir", default=EXAMPLE_PATHS["output_dir"],
        help="Root output directory",
    )
    p.add_argument(
        "--barcodes-whitelist", nargs="*", default=DEFAULT_BARCODE_WHITELIST,
        help="Process barcodes whose demux names start with these prefixes "
             "(empty list = process all). Default: barcode01–barcode06.",
    )
    p.add_argument(
        "--targets", nargs="+", default=DEFAULT_TARGETS,
        help="ARG targets for mapping / concatenation",
    )
    p.add_argument(
        "--min-cluster-size", type=int, default=DEFAULT_MIN_CLUSTER_SIZE,
        help="Minimum cluster size for Step 3",
    )
    p.add_argument(
        "--skip-existing-clusters", action="store_true",
        help="Skip clustering when clustersNNN directory already exists",
    )
    p.add_argument(
        "--cutadapt-bin", default="cutadapt",
        help="cutadapt executable",
    )
    p.add_argument(
        "--start-at", choices=["1", "2", "3"], default="1",
        help="Start at pipeline step (1=demux, 2=per-barcode, 3=cluster)",
    )
    return p.parse_args()


def main():
    args = parse_args()
    output_dir = Path(args.output_dir).expanduser().resolve()
    paths = build_paths(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    args.input = Path(args.input).expanduser()
    args.barcodes = Path(args.barcodes).expanduser()
    args.primers_and_y = Path(args.primers_and_y).expanduser()
    args.arg_sequences = Path(args.arg_sequences).expanduser()

    def resolve_existing(path):
        if path.is_file():
            return path
        alt = REPO_ROOT / path
        if alt.is_file():
            return alt
        raise FileNotFoundError(f"Input file not found: {path}")

    args.input = resolve_existing(args.input)
    args.barcodes = resolve_existing(args.barcodes)
    args.primers_and_y = resolve_existing(args.primers_and_y)
    args.arg_sequences = resolve_existing(args.arg_sequences)

    start = int(args.start_at)
    demux_dir = paths["demux_dir"]
    arg_concat_dir = paths["arg_concat_dir"]

    if start <= 1:
        demux_dir = step1_split_and_demux(args, paths)
    if start <= 2:
        if not demux_dir.is_dir():
            raise FileNotFoundError(f"Demux directory not found: {demux_dir}")
        arg_concat_dir = step2_per_barcode(args, paths, demux_dir)
    if start <= 3:
        if not arg_concat_dir.is_dir():
            raise FileNotFoundError(f"Concat directory not found: {arg_concat_dir}")
        step3_cluster_polish(args, arg_concat_dir)

    print("\n" + "=" * 70)
    print("PIPELINE COMPLETED SUCCESSFULLY")
    print("=" * 70)
    print(f"Outputs under: {output_dir}")


if __name__ == "__main__":
    main()
