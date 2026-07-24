#!/usr/bin/env python3

"""
Cutadapt Demultiplexing

Code autoannotated using Composer v2.5.

This code uses cutadapt to split and separate reads by barcode.

Author: O'Brien, M.E., Fuhrmeister, E.R., Marchand, J.A.
Date: 2026

"""
import os
import sys
import argparse
import subprocess
from datetime import datetime

from lib.default_params import CUTADAPT_DEMUX


def parse_forward_barcodes(barcode_file):
    forward_barcodes = {}
    with open(barcode_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                name = line.strip()[1:]
                seq = next(f).strip()
                if 'FWD' in name:
                    forward_barcodes[name] = seq
    return forward_barcodes


def detect_sequence_format(path):
    try:
        with open(path, 'r') as f:
            first = f.readline().strip()
        if first.startswith('@'):
            return 'fastq'
        if first.startswith('>'):
            return 'fasta'
    except FileNotFoundError:
        pass
    return 'fastq'


def count_reads(path, fmt):
    try:
        if fmt == 'fastq':
            with open(path, 'r') as f:
                lines = sum(1 for _ in f)
            return lines // 4
        else:
            count = 0
            with open(path, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        count += 1
            return count
    except FileNotFoundError:
        return 0


def run_cutadapt_batch(
    cutadapt_bin,
    input_fastq,
    named_adapters,
    output_template,
    unassigned_output,
    error_rate,
    min_overlap,
    min_len,
    max_len,
    internal_match,
    no_trim,
    times,
    output_format,
):
    cmd = [cutadapt_bin]

    # Adapter placement
    for adapter_name, adapter_seq in named_adapters:
        if internal_match:
            cmd += ['-b', f'{adapter_name}={adapter_seq}']
        else:
            cmd += ['-g', f'{adapter_name}={adapter_seq}']

    # Trimming behavior
    if no_trim:
        cmd.append('--no-trim')

    # Allow repeated occurrences
    if times is not None and times > 0:
        cmd += ['--times', str(times)]

    # Core filters
    if min_len is not None:
        cmd += ['-m', str(min_len)]
    if max_len is not None:
        cmd += ['-M', str(max_len)]
    cmd += ['-e', str(error_rate)]
    if min_overlap is not None:
        cmd += ['-O', str(min_overlap)]

    # Set output format to match input (cutadapt uses --fasta/--fastq)
    if output_format == 'fasta':
        cmd += ['--fasta']
    elif output_format == 'fastq':
        cmd += ['--fastq']

    # Outputs
    cmd += ['-o', output_template, '--untrimmed-output', unassigned_output, input_fastq]

    try:
        completed = subprocess.run(cmd, check=True, capture_output=True, text=True)
        return completed.stdout, ' '.join(cmd)
    except subprocess.CalledProcessError as e:
        # Surface cutadapt's stderr to help debugging when invoked via pipeline
        sys.stderr.write("\n[cutadapt error] Command failed:\n")
        sys.stderr.write(' '.join(cmd) + "\n")
        if e.stdout:
            sys.stderr.write("--- cutadapt STDOUT ---\n" + e.stdout + "\n")
        if e.stderr:
            sys.stderr.write("--- cutadapt STDERR ---\n" + e.stderr + "\n")
        raise


def main():
    ap = argparse.ArgumentParser(description='Demultiplex reads by forward primers using cutadapt.')
    ap.add_argument('--input', required=True, help='Input FASTA/FASTQ file')
    ap.add_argument('--barcodes', required=True, help='FASTA with barcodes (FWD entries used)')
    ap.add_argument('--output-dir', required=True, help='Output directory')
    ap.add_argument('--cutadapt-bin', default=CUTADAPT_DEMUX['cutadapt_bin'],
                    help=f"cutadapt executable (default: {CUTADAPT_DEMUX['cutadapt_bin']})")
    ap.add_argument('--error-rate', type=float, default=CUTADAPT_DEMUX['error_rate'],
                    help=f"Error rate (default: {CUTADAPT_DEMUX['error_rate']})")
    ap.add_argument('--min-overlap', type=int, default=CUTADAPT_DEMUX['min_overlap'],
                    help=f"Minimum overlap (default: {CUTADAPT_DEMUX['min_overlap']})")
    ap.add_argument('--min-len', type=int, default=CUTADAPT_DEMUX['min_len'],
                    help=f"Minimum length filter (default: {CUTADAPT_DEMUX['min_len']})")
    ap.add_argument('--max-len', type=int, default=CUTADAPT_DEMUX['max_len'],
                    help=f"Maximum length filter (default: {CUTADAPT_DEMUX['max_len']})")
    ap.add_argument('--internal', action='store_true', default=CUTADAPT_DEMUX['internal'],
                    help='Match adapters internally (-b) instead of anchored at start (-g)')
    ap.add_argument('--no-trim', action='store_true', default=CUTADAPT_DEMUX['no_trim'],
                    help='Do not trim matched barcode')
    ap.add_argument('--times', type=int, default=CUTADAPT_DEMUX['times'],
                    help='Remove up to N occurrences (concatemers)')

    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    forward_barcodes = parse_forward_barcodes(args.barcodes)
    if not forward_barcodes:
        print("No FWD barcodes found in barcode FASTA headers.", file=sys.stderr)
        sys.exit(1)

    named_adapters = list(forward_barcodes.items())
    input_fmt = detect_sequence_format(args.input)
    out_ext = 'fastq' if input_fmt == 'fastq' else 'fasta'
    output_template = os.path.join(args.output_dir, f'{{name}}.{out_ext}')
    unsorted_path = os.path.join(args.output_dir, f'unsorted.{out_ext}')

    stdout, cmd_str = run_cutadapt_batch(
        cutadapt_bin=args.cutadapt_bin,
        input_fastq=args.input,
        named_adapters=named_adapters,
        output_template=output_template,
        unassigned_output=unsorted_path,
        error_rate=args.error_rate,
        min_overlap=args.min_overlap,
        min_len=args.min_len,
        max_len=args.max_len,
        internal_match=args.internal,
        no_trim=args.no_trim,
        times=args.times,
        output_format=input_fmt,
    )

    # Write concise summary with counts
    summary_lines = []
    summary_lines.append(f'Demultiplex summary - {datetime.now().isoformat()}')
    summary_lines.append('')
    summary_lines.append('Parameters:')
    summary_lines.append(f'  input_file: {args.input}')
    summary_lines.append(f'  barcode_file: {args.barcodes}')
    summary_lines.append(f'  output_dir: {args.output_dir}')
    summary_lines.append(f'  io_format: {input_fmt}')
    summary_lines.append(f'  error_rate: {args.error_rate}')
    summary_lines.append(f'  min_overlap: {args.min_overlap}')
    summary_lines.append(f'  min_len: {args.min_len}')
    summary_lines.append(f'  max_len: {args.max_len}')
    summary_lines.append(f'  internal_match: {args.internal}')
    summary_lines.append(f'  no_trim: {args.no_trim}')
    summary_lines.append(f'  times: {args.times}')
    summary_lines.append('')
    summary_lines.append('Per-output read counts:')

    total_assigned = 0
    for fwd_name in forward_barcodes.keys():
        out_path = os.path.join(args.output_dir, f'{fwd_name}.{out_ext}')
        count = count_reads(out_path, input_fmt)
        total_assigned += count
        summary_lines.append(f'  {fwd_name}.{out_ext}: {count}')

    unsorted_count = count_reads(unsorted_path, input_fmt)
    summary_lines.append(f'  unsorted.{out_ext}: {unsorted_count}')
    summary_lines.append('')
    summary_lines.append('Totals:')
    summary_lines.append(f'  Assigned reads: {total_assigned}')
    summary_lines.append(f'  Unassigned (unsorted) reads: {unsorted_count}')
    summary_lines.append(f'  All outputs total: {total_assigned + unsorted_count}')
    summary_lines.append('')
    summary_lines.append('Command:')
    summary_lines.append(f'  {cmd_str}')

    summary_path = os.path.join(args.output_dir, 'demux_summary.txt')
    with open(summary_path, 'w') as f:
        f.write('\n'.join(summary_lines) + '\n')

    print(f"Summary written to: {summary_path}")


if __name__ == '__main__':
    main()
