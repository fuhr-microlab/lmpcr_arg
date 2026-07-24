#!/usr/bin/env python3
"""
Default parameters for Y-adapter ARG pipeline scripts.

Code autoannotated using Composer v2.5.

Each top-level dict is named for the script / stage that consumes it so shared
parameter names (e.g. identity, evalue) can differ by step. Example input
and output paths used in the README live in EXAMPLE_PATHS.

Author: O'Brien, M.E., Fuhrmeister, E.R., Marchand, J.A.
Date: 2026
"""

# =============================================================================
# Example input / output paths (relative to repository root)
# =============================================================================
EXAMPLE_PATHS = {
    'reads': None,  # FASTQ from nanopore basecalling; see publication for SRA accessions
    'barcodes': 'example/barcodes.fasta',
    'primers_and_y': 'example/primers_and_y.fasta',
    'arg_sequences': 'example/arg_sequences.fasta',
    'output_dir': 'output',
    'splittr_dir': 'output/splittr',
    'demux_dir': 'output/demux',
    'split_dir': 'output/split',
    'arg_filter_dir': 'output/arg_filter',
    'arg_mapping_dir': 'output/arg_mapping',
    'arg_concat_dir': 'output/arg_concat',
    'blast_db_dir': 'output/tmp_blast_db',
}

# =============================================================================
# adapter_splitter.py
# Presets: pre_demux (Step 1A) and post_demux (Step 2A)
# pre_demux values match the publication demux-stage settings
# =============================================================================
ADAPTER_SPLITTER = {
    'pre_demux': {
        'BLAST_PARAMS': {
            'identity': 80,
            'query_coverage': 90,
            'evalue': 1e-5,
            'word_size': 4,
            'dust': 'no',
            'soft_masking': 'false',
            'gapopen': 5,
            'gapextend': 2,
            'penalty': -1,
            'reward': 1,
            'max_target_seqs': 1000,
            'task': 'blastn-short',
        },
        'PARSING_PARAMS': {
            'identity_threshold': 80,
            'min_primer_overlap': 0.8,
        },
        'SPLITTING_PARAMS': {
            'min_fragment_length': 1000,
            'split_before_primer': True,
        },
        'DEBUG_PARAMS': {
            'save_blast_outputs': False,
            'save_round_outputs': True,
            'use_reversed_blast': False,
            'verbose': False,
        },
    },
    'post_demux': {
        'BLAST_PARAMS': {
            'identity': 80,
            'query_coverage': 90,
            'evalue': 1e-5,
            'word_size': 4,
            'dust': 'no',
            'soft_masking': 'false',
            'gapopen': 5,
            'gapextend': 2,
            'penalty': -1,
            'reward': 1,
            'max_target_seqs': 1000,
            'task': 'blastn-short',
        },
        'PARSING_PARAMS': {
            'identity_threshold': 80,
            'min_primer_overlap': 0.9,
        },
        'SPLITTING_PARAMS': {
            'min_fragment_length': 1000,
            'split_before_primer': True,
        },
        'DEBUG_PARAMS': {
            'save_blast_outputs': False,
            'save_round_outputs': True,
            'use_reversed_blast': False,
            'verbose': False,
        },
    },
}

# =============================================================================
# adapter_filter.py  (Step 2B — multi-ARG hit filter)
# =============================================================================
ADAPTER_FILTER = {
    'BLAST_PARAMS': {
        'identity': 80,
        'query_coverage': 90,
        'evalue': 1e-5,
        'word_size': 4,
        'dust': 'no',
        'soft_masking': 'false',
        'gapopen': 5,
        'gapextend': 2,
        'penalty': -1,
        'reward': 1,
        'max_target_seqs': 1000,
        'task': 'blastn-short',
    },
    'PARSING_PARAMS': {
        'identity_threshold': 80,
        'min_primer_overlap': 0.2,
    },
    'FILTERING_PARAMS': {
        'max_hits_per_read': 1,
        'min_read_length': 0,
    },
    'DEBUG_PARAMS': {
        'save_blast_outputs': False,
        'verbose': False,
    },
}

# =============================================================================
# cutadapt_demux.py  (Step 1B — barcode demultiplexing)
# Values match the publication demux-stage settings
# =============================================================================
CUTADAPT_DEMUX = {
    'cutadapt_bin': 'cutadapt',
    'error_rate': 0.2,
    'min_overlap': 14,
    'min_len': 500,
    'max_len': 50000,
    'internal': False,
    'no_trim': True,
    'times': None,
}

# =============================================================================
# arg_filter_and_plot.py  (Step 2C — ARG mapping / position filtering)
# =============================================================================
ARG_FILTER_AND_PLOT = {
    'BLAST_PARAMS': {
        'identity': 80,
        'query_coverage': 60,
        'evalue': 1e-5,
        'word_size': 4,
        'dust': 'no',
        'soft_masking': 'false',
        'gapopen': 5,
        'gapextend': 2,
        'penalty': -1,
        'reward': 1,
        'max_target_seqs': 1000,
        'task': 'blastn-short',
    },
    'PARSING_PARAMS': {
        'identity_threshold': 60,
        'min_primer_overlap': 0.5,
        'max_distance_from_end': 300,
    },
    'FILTER_PARAMS': {
        'min_length': 1000,
        'max_length': 8000,
    },
    'DEBUG_PARAMS': {
        'save_blast_outputs': False,
        'verbose': False,
    },
}

# =============================================================================
# cluster_polish.py + lib/cluster_by_vsearch.py + lib/polish_by_racon.py
# =============================================================================
CLUSTER_POLISH = {
    # CLI / ReadClusterer defaults
    'coverage_threshold': 95,
    'identity_threshold': 90,
    'min_cluster_size': 1,
    'max_clusters': None,
    'max_read_length': None,
    'polish_clusters': False,
    'DEFAULT_CLUSTERING_ALGORITHM': 'vsearch',
    'DEFAULT_POLISHING_ALGORITHM': 'racon',
    # Passed through to VsearchClusterer
    'VSEARCH_PARAMS': {
        # clustering_algorithm: cluster_fast | cluster_smallmem | cluster_size | cluster_unoise
        'clustering_algorithm': 'cluster_fast',
        'cluster_size': 15,
        'identity_threshold': 0.85,
        'query_coverage': 0.9,
        'target_coverage': 0.9,
        'strand': 'both',
        'sizein': True,
        'sizeout': True,
    },
    # Passed through to RaconPolisher
    'RACON_PARAMS': {
        'match': 3,
        'mismatch': -5,
        'gap': -4,
        'window_length': 500,
        'quality_threshold': 10.0,
        'error_threshold': 0.3,
        'threads': 1,
        'include_unpolished': False,
        'split': False,
        'fragments': False,
        'max_reads_for_polishing': 50,
    },
    'DEBUG_PARAMS': {
        'save_blast_outputs': False,
        'verbose': False,
        'save_intermediate': False,
        'output_longest_reads': True,
    },
    'REDUNDANT_REMOVAL_PARAMS': {
        'enabled': False,
        'identity_threshold': 0.95,
        'query_coverage': 0.80,
        'target_coverage': 0.80,
        'strand': 'both',
        'output_file_suffix': '_deduplicated',
    },
    'CLUSTERING_ALGORITHMS': {
        'vsearch': 'cluster_by_vsearch',
    },
    'POLISHING_ALGORITHMS': {
        'racon': 'polish_by_racon',
        'simple': 'polish_simple',
    },
}

# Top-level aliases used by lib/cluster_by_vsearch.py and lib/polish_by_racon.py
VSEARCH_PARAMS = CLUSTER_POLISH['VSEARCH_PARAMS']
RACON_PARAMS = CLUSTER_POLISH['RACON_PARAMS']
DEBUG_PARAMS = CLUSTER_POLISH['DEBUG_PARAMS']
REDUNDANT_REMOVAL_PARAMS = CLUSTER_POLISH['REDUNDANT_REMOVAL_PARAMS']
DEFAULT_CLUSTERING_ALGORITHM = CLUSTER_POLISH['DEFAULT_CLUSTERING_ALGORITHM']
DEFAULT_POLISHING_ALGORITHM = CLUSTER_POLISH['DEFAULT_POLISHING_ALGORITHM']
CLUSTERING_ALGORITHMS = CLUSTER_POLISH['CLUSTERING_ALGORITHMS']
POLISHING_ALGORITHMS = CLUSTER_POLISH['POLISHING_ALGORITHMS']
