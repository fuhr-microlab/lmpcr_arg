# Y-Adapter ARG Nanopore Analysis

This repository accompanies the publication:

> Megan E. O’Brien, Bradie Ahern, Piper Brase, Sooyeol Kim, Caroline McCormack, Denise Garcia, Erika Keim, Erin Dahl, Matthew Feck, Kelly Kauber, Jorge A. Marchand, Rose S. Kantor, Kara L. Nelson, Amy J. Pickering, Breanna McArdle, Erica R. Fuhrmeister *Capturing genomic contexts of clinically relevant ARGs in wastewater by Ligation-Mediation PCR* (2026).

Scripts in this pipeline were autoannotated using Composer v2.5 for publication purposes. Annotations were manually reviewed for accuracy.

The workflow processes multiplexed nanopore amplicons targeting antibiotic resistance genes (ARGs) generated with Y-adapter / barcode chemistry. It is intended for analyses of the type described in the associated publication, not as a general-purpose toolkit for unrelated designs. For exact parameter sets used, please refer to the publication.

## Overview

This pipeline takes raw nanopore reads to demultiplexed, sample-specific clusters of ARGs. Libraries can contain concatemers, barcode-spanning chimeras, and multi-ARG chimeric inserts; the workflow addresses those by:

1. Splitting chimeric/concatemeric products (BLAST+), then demultiplexing by barcode (cutadapt)
2. Re-splitting each barcode bin on primers / Y-adapters (BLAST+), removing multi-ARG hits and mapping remaining reads to ARG targets (BLAST+), then concatenating start/end fragments per target
3. Clustering concatenated reads (VSEARCH) and optionally polishing consensus sequences (RACON + minimap2)

## Requirements

- Python 3.8+
- [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi) (`blastn`, `makeblastdb`)
- [cutadapt](https://cutadapt.readthedocs.io/)
- [VSEARCH](https://github.com/torognes/vsearch)
- [RACON](https://github.com/lbcb-sci/racon) (for `--polish`)
- [minimap2](https://github.com/lh3/minimap2) (required by RACON polishing to align reads before consensus)
- Biopython: `pip install biopython`

All external binaries assumed to be on `PATH`.

## Repository layout

```text
.
├── run_pipeline.py              # end-to-end orchestrator
├── adapter_splitter.py          # Steps 1A and 2A (use --preset)
├── cutadapt_demux.py            # Step 1B
├── adapter_filter.py            # Step 2B
├── arg_filter_and_plot.py       # Step 2C
├── cluster_polish.py            # Step 3
├── lib/
│   ├── default_params.py        # defaults + EXAMPLE_PATHS (all tunable params)
│   ├── cluster_by_vsearch.py    # helper for cluster_polish.py
│   └── polish_by_racon.py       # helper for cluster_polish.py
├── example/                     # example input files
└── output/                      # pipeline outputs
```

Default paths for the worked example are defined in `lib/default_params.py` → `EXAMPLE_PATHS`.

## Expected inputs

| Input | Used in | Description |
|-------|---------|-------------|
| FASTQ from nanopore basecalling¹ | Step 1 | Multiplexed raw reads |
| `example/barcodes.fasta` | Step 1 | Forward barcodes (headers contain `FWD`) |
| `example/primers_and_y.fasta` | Step 2A | Primer + Y-adapter sequences |
| `example/arg_sequences.fasta` | Step 2B–2C | ARG references for filtering / mapping |

¹ Inputs can be downloaded from the associated SRA — see the publication for SRA accessions.

Outputs go under `output/`. Run commands from the repository root. For an end-to-end run, use `run_pipeline.py` (see below).

## Default parameters vs CLI flags

All tunable defaults live in `lib/default_params.py`:

- `ADAPTER_SPLITTER['pre_demux']` / `ADAPTER_SPLITTER['post_demux']` — selected with `--preset`
- `CUTADAPT_DEMUX`, `ADAPTER_FILTER`, `CLUSTER_POLISH`
- `EXAMPLE_PATHS` — example input/output path strings

**Required to run:** positional / named input-output paths.  
**Optional flags:** override a value from `lib/default_params.py`.  
For exact parameter sets used in the publication, please refer to the publication.

---

## Pipeline overview

```text
Step 1 — Split concatemers and demultiplex
  FASTQ from nanopore basecalling
    → adapter_splitter.py --preset pre_demux  → output/splittr/
    → cutadapt_demux.py                       → output/demux/

Step 2 — Per-barcode split, ARG filter, map, concatenate
  output/demux/<barcode>
    → adapter_splitter.py --preset post_demux → output/split/
    → adapter_filter.py                       → output/arg_filter/
    → arg_filter_and_plot.py                  → output/arg_mapping/
    → cat                                     → output/arg_concat/

Step 3 — Cluster and polish
  output/arg_concat/..._combined.fasta
    → cluster_polish.py                       → clusters/
```

### End-to-end run

```bash
python run_pipeline.py \
  --input /PATH/to/reads.fastq \
  --barcodes example/barcodes.fasta \
  --primers-and-y example/primers_and_y.fasta \
  --arg-sequences example/arg_sequences.fasta \
  --output-dir output
```

Optional: `--barcodes-whitelist barcode01 barcode02 ...` (default: barcode01–barcode06), `--targets ctx-m KPC OXA qnrS`, `--skip-existing-clusters`, `--min-cluster-size 15`. For exact parameter sets used, please refer to the publication.

---

## Step 1 — Split concatemers and demultiplex

### Step 1A. Concatemer / chimera splitting

| | |
|---|---|
| **Script** | `adapter_splitter.py` |
| **Preset** | `--preset pre_demux` |
| **Input** | FASTQ from nanopore basecalling |
| **Reference** | `example/barcodes.fasta` |
| **Output** | `output/splittr/split_reads.fasta` |
| **Defaults** | `ADAPTER_SPLITTER['pre_demux']` |

**Required:**

```bash
mkdir -p output/splittr output/tmp_blast_db

python adapter_splitter.py \
  /PATH/to/reads.fastq \
  example/barcodes.fasta \
  output/splittr/split_reads.fasta \
  --preset pre_demux
```

**Optional flags:** `--blast-db-dir`, `--identity`, `--coverage`, `--evalue`, `--min-primer-overlap`, `--min-fragment-length`, `--verbose`, `--debug`, and related BLAST options. See `python adapter_splitter.py -h`.

### Step 1B. Barcode demultiplexing

| | |
|---|---|
| **Script** | `cutadapt_demux.py` |
| **Input** | `output/splittr/split_reads.fasta` |
| **Reference** | `example/barcodes.fasta` (headers containing `FWD`) |
| **Output** | `output/demux/<FWD_barcode_name>.fasta`, `unsorted.fasta`, `demux_summary.txt` |
| **Defaults** | `CUTADAPT_DEMUX` |

**Required:**

```bash
mkdir -p output/demux

python cutadapt_demux.py \
  --input output/splittr/split_reads.fasta \
  --barcodes example/barcodes.fasta \
  --output-dir output/demux
```

**Optional flags:** `--cutadapt-bin`, `--error-rate`, `--min-overlap`, `--min-len`, `--max-len`, `--internal`, `--no-trim`, `--times`.

Demux output names follow FASTA headers (e.g. `barcode01_FWD.fasta`). Examples below use `barcode01_FWD`; rename paths to match your files.

---

## Step 2 — Per-barcode split, ARG filter, map, concatenate

ARG targets used in the publication include `ctx-m`, `KPC`, `OXA`, and `qnrS`.

### Step 2A. Primer / Y-adapter splitting

| | |
|---|---|
| **Script** | `adapter_splitter.py` |
| **Preset** | `--preset post_demux` |
| **Input** | `output/demux/barcode01_FWD.fasta` |
| **Reference** | `example/primers_and_y.fasta` |
| **Output** | `output/split/barcode01_split.fasta` |
| **Defaults** | `ADAPTER_SPLITTER['post_demux']` |

**Required:**

```bash
mkdir -p output/split

python adapter_splitter.py \
  output/demux/barcode01_FWD.fasta \
  example/primers_and_y.fasta \
  output/split/barcode01_split.fasta \
  --preset post_demux
```

### Step 2B. Remove reads with multiple ARG hits

| | |
|---|---|
| **Script** | `adapter_filter.py` |
| **Input** | `output/split/barcode01_split.fasta` |
| **Reference** | `example/arg_sequences.fasta` |
| **Output** | `output/arg_filter/barcode01_filtered.fasta` (+ summary) |
| **Defaults** | `ADAPTER_FILTER` |

**Required:**

```bash
mkdir -p output/arg_filter

python adapter_filter.py \
  output/split/barcode01_split.fasta \
  example/arg_sequences.fasta \
  output/arg_filter/barcode01_filtered.fasta
```

**Optional flags:** `--max-hits`, `--min-length`, `--identity`, `--coverage`, `--evalue`, `--min-primer-overlap`, `--verbose`, `--debug`.

### Step 2C. ARG mapping / position filtering (per target)

| | |
|---|---|
| **Script** | `arg_filter_and_plot.py` |
| **Input** | `output/arg_filter/barcode01_filtered.fasta` |
| **Reference** | `example/arg_sequences.fasta` |
| **Output** | `output/arg_mapping/<TARGET>/barcode01/start_<TARGET>.fasta` and `end_<TARGET>.fasta` |

```bash
TARGET=ctx-m
mkdir -p "output/arg_mapping/${TARGET}/barcode01"

python arg_filter_and_plot.py \
  output/arg_filter/barcode01_filtered.fasta \
  example/arg_sequences.fasta \
  "output/arg_mapping/${TARGET}/barcode01" \
  --target "${TARGET}"
```

### Step 2D. Concatenate start + end fragments

```bash
TARGET=ctx-m
mkdir -p "output/arg_concat/${TARGET}"

cat \
  "output/arg_mapping/${TARGET}/barcode01/start_${TARGET}.fasta" \
  "output/arg_mapping/${TARGET}/barcode01/end_${TARGET}.fasta" \
  > "output/arg_concat/${TARGET}/barcode01_combined.fasta"
```

Repeat Steps 2C–2D for `KPC`, `OXA`, and `qnrS` as needed.

---

## Step 3 — Cluster and polish

| | |
|---|---|
| **Script** | `cluster_polish.py` |
| **Input** | Combined start+end FASTA |
| **Output** | Cluster directory under `output/` |
| **Defaults** | `CLUSTER_POLISH` |
| **Helpers** | `lib/cluster_by_vsearch.py`, `lib/polish_by_racon.py` |

**Required:**

```bash
mkdir -p output/arg_concat/ctx-m/clusters001

python cluster_polish.py \
  output/arg_concat/ctx-m/barcode01_combined.fasta \
  output/arg_concat/ctx-m/clusters001
```

**Optional flags:** `--coverage`, `--identity`, `--min-cluster-size`, `--polish` / `-p`, `--verbose`, `--debug`, `--remove-redundant-consensus`.

For exact parameter sets used, please refer to the publication.

---

## Example directory layout after a full run

```text
example/
  barcodes.fasta
  primers_and_y.fasta
  arg_sequences.fasta
output/
  splittr/split_reads.fasta
  demux/
    barcode01_FWD.fasta
    demux_summary.txt
  split/barcode01_split.fasta
  arg_filter/barcode01_filtered.fasta
  arg_mapping/ctx-m/barcode01/...
  arg_concat/ctx-m/barcode01_combined.fasta
```

---

## Scripts in this repository

| Script | Role |
|--------|------|
| `run_pipeline.py` | End-to-end orchestrator (Steps 1–3) |
| `adapter_splitter.py` | BLAST concatemer splitter (`--preset pre_demux` or `post_demux`) |
| `cutadapt_demux.py` | Barcode demultiplexing via cutadapt |
| `adapter_filter.py` | Removes reads with too many ARG hits |
| `arg_filter_and_plot.py` | Maps filtered reads to a target ARG; writes start/end FASTAs |
| `cluster_polish.py` | VSEARCH clustering + optional RACON/minimap2 polishing |
| `lib/default_params.py` | All tunable defaults and `EXAMPLE_PATHS` |
| `lib/cluster_by_vsearch.py` | Clustering helper for `cluster_polish.py` |
| `lib/polish_by_racon.py` | Polishing helper for `cluster_polish.py` (requires RACON + minimap2) |

## Citation

This repository was created for the associated publication:

> Megan E. O’Brien, Bradie Ahern, Piper Brase, Sooyeol Kim, Caroline McCormack, Denise Garcia, Erika Keim, Erin Dahl, Matthew Feck, Kelly Kauber, Jorge A. Marchand, Rose S. Kantor, Kara L. Nelson, Amy J. Pickering, Breanna McArdle, Erica R. Fuhrmeister *Capturing genomic contexts of clinically relevant ARGs in wastewater by Ligation-Mediation PCR* (2026).
