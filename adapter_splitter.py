#!/usr/bin/env python3
"""
Nanopore Concatemer Splitter

Code autoannotated using Composer v2.5.

This script processes nanopore sequencing data containing concatemers by identifying
primer / barcode / adapter sequences and splitting reads at those matches. The same
script is used for pre-demux (barcode) splitting and post-demux (primer / Y-adapter)
splitting; choose a parameter preset with --preset.

Author: O'Brien, M.E., Fuhrmeister, E.R., Marchand, J.A.
Date: 2024
"""

import os
import sys
import subprocess
import tempfile
import shutil
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse
import logging

from lib.default_params import ADAPTER_SPLITTER

# Module-level aliases updated when a preset is loaded (used by process helpers)
BLAST_PARAMS = ADAPTER_SPLITTER['pre_demux']['BLAST_PARAMS']
PARSING_PARAMS = ADAPTER_SPLITTER['pre_demux']['PARSING_PARAMS']
SPLITTING_PARAMS = ADAPTER_SPLITTER['pre_demux']['SPLITTING_PARAMS']
DEBUG_PARAMS = ADAPTER_SPLITTER['pre_demux']['DEBUG_PARAMS']


def apply_preset(preset_name):
    """Load ADAPTER_SPLITTER[preset_name] into module-level param dicts."""
    global BLAST_PARAMS, PARSING_PARAMS, SPLITTING_PARAMS, DEBUG_PARAMS
    if preset_name not in ADAPTER_SPLITTER:
        raise ValueError(f"Unknown preset '{preset_name}'. Choose from: {list(ADAPTER_SPLITTER)}")
    cfg = ADAPTER_SPLITTER[preset_name]
    BLAST_PARAMS = cfg['BLAST_PARAMS']
    PARSING_PARAMS = cfg['PARSING_PARAMS']
    SPLITTING_PARAMS = cfg['SPLITTING_PARAMS']
    DEBUG_PARAMS = cfg['DEBUG_PARAMS']
    return cfg

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class NanoporeConcatemerSplitter:
    def __init__(self, input_file, primer_file, output_file, blast_db_dir=None, save_blast_outputs=None, use_reversed_blast=None, verbose=None, 
                 identity=None, coverage=None, evalue=None, word_size=None, dust=None, soft_masking=None, gapopen=None, gapextend=None, 
                 penalty=None, reward=None, max_target_seqs=None, task=None, identity_threshold=None, min_primer_overlap=None, 
                 min_fragment_length=None, split_before_primer=None, save_round_outputs=None, preset='pre_demux'):
        """
        Initialize the splitter with input files and parameters.
        
        Args:
            preset (str): Parameter set from lib.default_params.ADAPTER_SPLITTER
                ('pre_demux' or 'post_demux'). CLI flags override preset values.
        """
        apply_preset(preset)

        self.input_file = Path(input_file)
        self.primer_file = Path(primer_file)
        self.output_file = Path(output_file)
        self.blast_db_dir = Path(blast_db_dir) if blast_db_dir else Path.cwd()
        
        # Use provided parameters or defaults from DEBUG_PARAMS
        self.save_blast_outputs = save_blast_outputs if save_blast_outputs is not None else DEBUG_PARAMS['save_blast_outputs']
        self.use_reversed_blast = use_reversed_blast if use_reversed_blast is not None else DEBUG_PARAMS['use_reversed_blast']
        self.verbose = verbose if verbose is not None else DEBUG_PARAMS['verbose']
        self.save_round_outputs = save_round_outputs if save_round_outputs is not None else DEBUG_PARAMS['save_round_outputs']
        
        # Use parameters from the selected preset
        self.blast_params = BLAST_PARAMS.copy()
        self.parsing_params = PARSING_PARAMS.copy()
        self.splitting_params = SPLITTING_PARAMS.copy()
        
        # Override with provided parameters
        if identity is not None:
            self.blast_params['identity'] = identity
        if coverage is not None:
            self.blast_params['query_coverage'] = coverage
        if evalue is not None:
            self.blast_params['evalue'] = evalue
        if word_size is not None:
            self.blast_params['word_size'] = word_size
        if dust is not None:
            self.blast_params['dust'] = dust
        if soft_masking is not None:
            self.blast_params['soft_masking'] = soft_masking
        if gapopen is not None:
            self.blast_params['gapopen'] = gapopen
        if gapextend is not None:
            self.blast_params['gapextend'] = gapextend
        if penalty is not None:
            self.blast_params['penalty'] = penalty
        if reward is not None:
            self.blast_params['reward'] = reward
        if max_target_seqs is not None:
            self.blast_params['max_target_seqs'] = max_target_seqs
        if task is not None:
            self.blast_params['task'] = task
        if identity_threshold is not None:
            self.parsing_params['identity_threshold'] = identity_threshold
        if min_primer_overlap is not None:
            self.parsing_params['min_primer_overlap'] = min_primer_overlap
        if min_fragment_length is not None:
            self.splitting_params['min_fragment_length'] = min_fragment_length
        if split_before_primer is not None:
            self.splitting_params['split_before_primer'] = split_before_primer
        
        # Validate input files
        self._validate_inputs()
        
    def _validate_inputs(self):
        """Validate that input files exist and are readable."""
        if not self.input_file.exists():
            raise FileNotFoundError(f"Input file not found: {self.input_file}")
        if not self.primer_file.exists():
            raise FileNotFoundError(f"Primer file not found: {self.primer_file}")
        
        # Check if input is FASTA or FASTQ
        self.input_format = self._detect_format(self.input_file)
        logger.info(f"Detected input format: {self.input_format}")
        
    def _detect_format(self, file_path):
        """Detect if file is FASTA or FASTQ format."""
        with open(file_path, 'r') as f:
            first_line = f.readline().strip()
            if first_line.startswith('@'):
                return 'fastq'
            elif first_line.startswith('>'):
                return 'fasta'
            else:
                raise ValueError(f"Unknown file format: {file_path}")
    
    def _create_blast_db(self, primer_file, db_name):
        """Create a BLAST database from primer sequences using simple approach."""
        # Create the database in the blast_db_dir with just the base name
        db_base_path = self.blast_db_dir / db_name
        
        logger.info(f"Creating BLAST database: {db_base_path}")
        logger.info(f"Input primer file: {primer_file}")
        logger.info(f"Database directory: {self.blast_db_dir}")
        
        # Check if primer file exists and has content
        if not primer_file.exists():
            raise FileNotFoundError(f"Primer file not found: {primer_file}")
        
        # Read and log primer file content
        with open(primer_file, 'r') as f:
            primer_content = f.read()
            logger.info(f"Primer file content:\n{primer_content}")
        
        # Use simple approach like simple_blast.py
        cmd = [
            'makeblastdb',
            '-in', str(primer_file),
            '-dbtype', 'nucl',
            '-out', str(db_base_path)
        ]
        
        logger.info(f"BLAST command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info(f"BLAST database creation STDOUT: {result.stdout}")
            logger.info(f"BLAST database creation STDERR: {result.stderr}")
            logger.info(f"Created BLAST database: {db_base_path}")
            
            # Check if database files were created
            db_files = list(self.blast_db_dir.glob(f"{db_name}.*"))
            logger.info(f"Database files created: {[f.name for f in db_files]}")
            
            return db_base_path
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to create BLAST database: {e}")
            logger.error(f"STDOUT: {e.stdout}")
            logger.error(f"STDERR: {e.stderr}")
            raise
    
    def _run_blast(self, query_file, db_path, output_file, primer_name=None):
        """Run BLAST search with simple approach (like simple_blast.py)."""
        
        # Check if query file exists and has content
        if not query_file.exists():
            raise FileNotFoundError(f"Query file not found: {query_file}")
        
        # Check query file size and content
        query_size = query_file.stat().st_size
        logger.info(f"Query file size: {query_size} bytes")
        
        # Read first few lines of query file
        with open(query_file, 'r') as f:
            first_lines = [f.readline().strip() for _ in range(6)]
            logger.info(f"Query file first lines:\n{chr(10).join(first_lines)}")
        
        # Check if database files exist
        db_files = list(self.blast_db_dir.glob(f"{db_path.name}.*"))
        logger.info(f"Database files for {db_path.name}: {[f.name for f in db_files]}")
        
        # Build command with BLAST_PARAMS for higher sensitivity on short primers
        p = self.blast_params
        cmd = [
            'blastn',
            '-query', str(query_file),
            '-db', str(db_path),
            '-task', str(p.get('task', 'blastn-short')),
            '-evalue', str(p.get('evalue', 1e-5)),
            '-reward', str(p.get('reward', 1)),
            '-penalty', str(p.get('penalty', -1)),
            '-gapopen', str(p.get('gapopen', 5)),
            '-gapextend', str(p.get('gapextend', 2)),
            '-dust', str(p.get('dust', 'no')),
            '-soft_masking', str(p.get('soft_masking', 'false')),
            '-perc_identity', str(p.get('identity', 80)),
            '-max_target_seqs', str(p.get('max_target_seqs', 1000)),
            '-strand', 'both',
            '-out', str(output_file),
            '-outfmt', '10 qseqid sseqid pident qcovs qstart qend sstart send evalue bitscore sstrand'
        ]
        
        logger.info(f"BLAST command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            logger.info(f"BLAST search STDOUT: {result.stdout}")
            logger.info(f"BLAST search STDERR: {result.stderr}")
            logger.info(f"BLAST search completed: {output_file}")
            
            # Check output file size
            if output_file.exists():
                output_size = output_file.stat().st_size
                logger.info(f"BLAST output file size: {output_size} bytes")
                
                if output_size == 0:
                    logger.warning("BLAST output file is empty - no hits found or error occurred")
                    # Read first few lines to see if there's any content
                    with open(output_file, 'r') as f:
                        content = f.read(100)
                        logger.info(f"BLAST output content (first 100 chars): '{content}'")
                else:
                    # Read first few lines of output
                    with open(output_file, 'r') as f:
                        first_lines = [f.readline().strip() for _ in range(5)]
                        logger.info(f"BLAST output first lines:\n{chr(10).join(first_lines)}")
            else:
                logger.error("BLAST output file was not created")
            
            # Save BLAST output to local directory if requested
            if self.save_blast_outputs and primer_name:
                # Create blast_outputs directory if it doesn't exist
                blast_outputs_dir = Path("blast_outputs")
                blast_outputs_dir.mkdir(exist_ok=True)
                
                # Copy BLAST output to local directory
                local_blast_output = blast_outputs_dir / f"blast_{primer_name}_output.txt"
                shutil.copy2(output_file, local_blast_output)
                logger.info(f"Saved BLAST output to: {local_blast_output}")
                
                # Also save the query file for reference
                local_query_file = blast_outputs_dir / f"blast_{primer_name}_query.fasta"
                shutil.copy2(query_file, local_query_file)
                logger.info(f"Saved query file to: {local_query_file}")
                
                # Create summary file with BLAST parameters
                summary_file = blast_outputs_dir / f"blast_{primer_name}_summary.txt"
                with open(summary_file, 'w') as f:
                    f.write(f"BLAST Search Summary for Primer: {primer_name}\n")
                    f.write("=" * 50 + "\n\n")
                    f.write("BLAST Parameters:\n")
                    f.write(f"  Output format: 10 (CSV)\n")
                    f.write("\nCommand executed:\n")
                    f.write(f"  {' '.join(cmd)}\n\n")
                    f.write("Output format: CSV (format 10)\n")
                    f.write(f"\nQuery file size: {query_size} bytes\n")
                    f.write(f"Output file size: {output_file.stat().st_size if output_file.exists() else 0} bytes\n")
                    f.write(f"\nBLAST return code: {result.returncode}\n")
                    f.write(f"BLAST STDOUT: {result.stdout}\n")
                    f.write(f"BLAST STDERR: {result.stderr}\n")
                logger.info(f"Saved BLAST summary to: {summary_file}")
                
        except subprocess.CalledProcessError as e:
            logger.error(f"BLAST search failed: {e}")
            logger.error(f"STDOUT: {e.stdout}")
            logger.error(f"STDERR: {e.stderr}")
            raise

    def _run_blast_reversed(self, reads_file, primer_file, output_file, primer_name=None):
        """Run BLAST search with reversed query-subject order (reads as DB, primer as query)."""
        
        # Check if files exist and have content
        if not reads_file.exists():
            raise FileNotFoundError(f"Reads file not found: {reads_file}")
        if not primer_file.exists():
            raise FileNotFoundError(f"Primer file not found: {primer_file}")
        
        # Check file sizes
        reads_size = reads_file.stat().st_size
        primer_size = primer_file.stat().st_size
        logger.info(f"Reads file size: {reads_size} bytes")
        logger.info(f"Primer file size: {primer_size} bytes")
        
        # Read first few lines of files
        with open(reads_file, 'r') as f:
            first_lines = [f.readline().strip() for _ in range(4)]
            logger.info(f"Reads file first lines:\n{chr(10).join(first_lines)}")
        
        with open(primer_file, 'r') as f:
            primer_content = f.read()
            logger.info(f"Primer file content:\n{primer_content}")
        
        # Create database from reads using simple approach
        reads_db_path = reads_file.parent / f"{reads_file.stem}_db"
        db_cmd = [
            'makeblastdb',
            '-in', str(reads_file),
            '-dbtype', 'nucl',
            '-out', str(reads_db_path)
        ]
        
        logger.info(f"Creating reads database: {' '.join(db_cmd)}")
        try:
            db_result = subprocess.run(db_cmd, capture_output=True, text=True)
            logger.info(f"Database creation STDOUT: {db_result.stdout}")
            logger.info(f"Database creation STDERR: {db_result.stderr}")
            if db_result.returncode != 0:
                logger.error(f"Database creation failed: {db_result.stderr}")
                raise subprocess.CalledProcessError(db_result.returncode, db_cmd)
            logger.info("✓ Reads database created successfully")
            
            # Check if database files were created
            db_files = list(reads_file.parent.glob(f"{reads_db_path.name}.*"))
            logger.info(f"Database files created: {[f.name for f in db_files]}")
            
        except Exception as e:
            logger.error(f"Database creation error: {e}")
            raise
        
        # Run BLAST search with BLAST_PARAMS in reversed mode
        p = self.blast_params
        cmd = [
            'blastn',
            '-query', str(primer_file),
            '-db', str(reads_db_path),
            '-task', str(p.get('task', 'blastn-short')),
            '-evalue', str(p.get('evalue', 1e-5)),
            '-word_size', str(p.get('word_size', 4)),
            '-reward', str(p.get('reward', 1)),
            '-penalty', str(p.get('penalty', -1)),
            '-gapopen', str(p.get('gapopen', 5)),
            '-gapextend', str(p.get('gapextend', 2)),
            '-dust', str(p.get('dust', 'no')),
            '-soft_masking', str(p.get('soft_masking', 'false')),
            '-perc_identity', str(p.get('identity', 80)),
            '-max_target_seqs', str(p.get('max_target_seqs', 1000)),
            '-strand', 'both',
            '-out', str(output_file),
            '-outfmt', '10 qseqid sseqid pident qcovs qstart qend sstart send evalue bitscore sstrand'
        ]
        
        logger.info(f"BLAST command (reversed): {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            logger.info(f"BLAST search STDOUT: {result.stdout}")
            logger.info(f"BLAST search STDERR: {result.stderr}")
            logger.info(f"BLAST return code: {result.returncode}")
            logger.info(f"BLAST search completed: {output_file}")
            
            # Check output file size
            if output_file.exists():
                output_size = output_file.stat().st_size
                logger.info(f"BLAST output file size: {output_size} bytes")
                
                if output_size == 0:
                    logger.warning("BLAST output file is empty - no hits found or error occurred")
                    # Read first few lines to see if there's any content
                    with open(output_file, 'r') as f:
                        content = f.read(100)
                        logger.info(f"BLAST output content (first 100 chars): '{content}'")
                else:
                    # Read first few lines of output
                    with open(output_file, 'r') as f:
                        first_lines = [f.readline().strip() for _ in range(5)]
                        logger.info(f"BLAST output first lines:\n{chr(10).join(first_lines)}")
            else:
                logger.error("BLAST output file was not created")
            
            # Save BLAST output to local directory if requested
            if self.save_blast_outputs and primer_name:
                # Create blast_outputs directory if it doesn't exist
                blast_outputs_dir = Path("blast_outputs")
                blast_outputs_dir.mkdir(exist_ok=True)
                
                # Copy BLAST output to local directory
                local_blast_output = blast_outputs_dir / f"blast_{primer_name}_output_reversed.txt"
                shutil.copy2(output_file, local_blast_output)
                logger.info(f"Saved BLAST output to: {local_blast_output}")
                
                # Also save the files for reference
                local_primer_file = blast_outputs_dir / f"blast_{primer_name}_primer.fasta"
                shutil.copy2(primer_file, local_primer_file)
                logger.info(f"Saved primer file to: {local_primer_file}")
                
                # Create summary file with BLAST parameters
                summary_file = blast_outputs_dir / f"blast_{primer_name}_summary_reversed.txt"
                with open(summary_file, 'w') as f:
                    f.write(f"BLAST Search Summary for Primer: {primer_name} (Reversed)\n")
                    f.write("=" * 60 + "\n\n")
                    f.write("BLAST Parameters:\n")
                    f.write(f"  Output format: 10 (CSV)\n")
                    f.write("\nCommand executed:\n")
                    f.write(f"  {' '.join(cmd)}\n\n")
                    f.write("Output format: CSV (format 10)\n")
                    f.write(f"\nReads file size: {reads_size} bytes\n")
                    f.write(f"Primer file size: {primer_size} bytes\n")
                    f.write(f"Output file size: {output_file.stat().st_size if output_file.exists() else 0} bytes\n")
                    f.write("\nNote: This uses REVERSED query-subject order (primer as query, reads as database)\n")
                    f.write(f"\nBLAST return code: {result.returncode}\n")
                    f.write(f"BLAST STDOUT: {result.stdout}\n")
                    f.write(f"BLAST STDERR: {result.stderr}\n")
                logger.info(f"Saved BLAST summary to: {summary_file}")
                
        except subprocess.CalledProcessError as e:
            logger.error(f"BLAST search failed: {e}")
            logger.error(f"STDOUT: {e.stdout}")
            logger.error(f"STDERR: {e.stderr}")
            raise
    
    def _calculate_primer_overlap(self, hit_start, hit_end, primer_length):
        """Calculate the fraction of primer that overlaps with the hit."""
        # Convert to 0-based coordinates if needed
        if hit_start > hit_end:
            hit_start, hit_end = hit_end, hit_start
        
        # Calculate overlap length
        overlap_length = hit_end - hit_start + 1
        
        # Calculate fraction of primer that overlaps
        overlap_fraction = overlap_length / primer_length
        
        return overlap_fraction, overlap_length
    
    def _parse_blast_results(self, blast_output, reversed_order=False):
        """Parse BLAST results and group by query sequence."""
        hits = {}
        
        if not os.path.exists(blast_output) or os.path.getsize(blast_output) == 0:
            logger.warning(f"BLAST output file is empty or doesn't exist: {blast_output}")
            return hits
        
        logger.info(f"Parsing BLAST results from: {blast_output}")
        logger.info(f"Reversed order: {reversed_order}")
        
        # Use parameters from the top of the script
        identity_threshold = self.parsing_params['identity_threshold']
        min_primer_overlap = self.parsing_params['min_primer_overlap']
        
        # Get primer lengths for overlap calculation and determine if primers are reverse
        primer_lengths = {}
        reverse_primers = set()
        for record in SeqIO.parse(self.primer_file, 'fasta'):
            primer_lengths[record.id] = len(record.seq)
            # Detect reverse primers by name patterns
            if '(-)' in record.id or record.id.endswith('(-)') or 'reverse' in record.id.lower():
                reverse_primers.add(record.id)
                logger.info(f"Detected reverse primer: {record.id}")
        
        logger.info(f"Primer lengths: {primer_lengths}")
        logger.info(f"Reverse primers: {reverse_primers}")
        
        with open(blast_output, 'r') as f:
            line_count = 0
            for line in f:
                line_count += 1
                line = line.strip()
                if not line:  # Skip empty lines
                    continue
                
                # Try to detect format - check if it's CSV (format 10) or tabular (format 6)
                if ',' in line:
                    # CSV format (format 10)
                    fields = line.split(',')
                    logger.debug(f"Detected CSV format, {len(fields)} fields")
                else:
                    # Tabular format (format 6)
                    fields = line.split('\t')
                    logger.debug(f"Detected tabular format, {len(fields)} fields")
                
                if len(fields) >= 10:  # Now expecting at least 10 fields with strand info
                    if reversed_order:
                        # In reversed order: qseqid=primer, sseqid=read
                        qseqid = fields[0]  # primer name
                        sseqid = fields[1]  # read name
                        pident = float(fields[2])
                        qcovs = float(fields[3])
                        qstart = int(fields[4])  # primer start
                        qend = int(fields[5])    # primer end
                        sstart = int(fields[6])  # read start
                        send = int(fields[7])    # read end
                        evalue = float(fields[8])
                        bitscore = float(fields[9]) if len(fields) > 9 else 0.0
                        sstrand = fields[10] if len(fields) > 10 else 'plus'  # strand info
                        
                        logger.debug(f"Parsed line {line_count}: primer={qseqid}, read={sseqid}, identity={pident}%, coverage={qcovs}%, read_pos={sstart}-{send}, strand={sstrand}")
                        
                        # Calculate primer overlap
                        if sseqid in primer_lengths:
                            primer_length = primer_lengths[sseqid]
                            overlap_fraction, overlap_length = self._calculate_primer_overlap(qstart, qend, primer_length)
                            
                            # Filter by identity and primer overlap
                            if pident >= identity_threshold and overlap_fraction >= min_primer_overlap:
                                if qseqid not in hits:  # Group by read (query)
                                    hits[qseqid] = []
                                # Determine if this is a forward or reverse primer
                                # Consider both primer name AND strand orientation for reverse complement reads
                                is_reverse_primer_by_name = sseqid in reverse_primers
                                
                                # For reverse complement reads, we need to consider the actual strand
                                # If the primer matches on the minus strand, it's actually a reverse complement match
                                is_minus_strand_match = sstrand == 'minus'
                                
                                # A primer is "forward" if:
                                # 1. It's a forward primer by name AND matches on plus strand, OR
                                # 2. It's a reverse primer by name AND matches on minus strand (reverse complement)
                                is_forward = (not is_reverse_primer_by_name and not is_minus_strand_match) or \
                                           (is_reverse_primer_by_name and is_minus_strand_match)
                                
                                is_reverse_primer = not is_forward
                                
                                hits[qseqid].append({
                                    'primer': sseqid,
                                    'start': read_start,  # Use read coordinates (query)
                                    'end': read_end,
                                    'identity': pident,
                                    'overlap_fraction': overlap_fraction,
                                    'overlap_length': overlap_length,
                                    'evalue': evalue,
                                    'strand': sstrand,
                                    'is_forward': is_forward,  # Based on primer name, not strand
                                    'is_reverse_primer': is_reverse_primer
                                })
                                if self.verbose:
                                    logger.debug(f"Added hit for read {qseqid} at position {read_start}-{read_end} (strand: {sstrand}, overlap: {overlap_fraction:.2f})")
                                    logger.debug(f"  Primer: {sseqid}, is_reverse_by_name: {is_reverse_primer_by_name}, is_minus_strand: {is_minus_strand_match}, is_forward: {is_forward}")
                            else:
                                if self.verbose:
                                    logger.debug(f"Hit filtered out: identity={pident}% (threshold={identity_threshold}%), overlap={overlap_fraction:.2f} (threshold={min_primer_overlap})")
                        else:
                            logger.warning(f"Primer {sseqid} not found in primer file")
                    else:
                        # Original order: qseqid=read, sseqid=primer
                        qseqid = fields[0]  # read name
                        sseqid = fields[1]  # primer name
                        pident = float(fields[2])
                        qcovs = float(fields[3])
                        qstart = int(fields[4])  # query start (read start) - CORRECT!
                        qend = int(fields[5])    # query end (read end) - CORRECT!
                        sstart = int(fields[6])  # subject start (primer start)
                        send = int(fields[7])    # subject end (primer end)
                        evalue = float(fields[8])
                        sstrand = fields[10] if len(fields) > 10 else 'plus'  # strand info
                        
                        # Use READ coordinates for split positions (qstart, qend)
                        read_start = qstart  # Use query coordinates for read position
                        read_end = qend
                        
                        logger.debug(f"Parsed line {line_count}: read={qseqid}, primer={sseqid}, identity={pident}%, coverage={qcovs}%, read_pos={read_start}-{read_end}, strand={sstrand}")
                        
                        # Calculate primer overlap
                        if sseqid in primer_lengths:
                            primer_length = primer_lengths[sseqid]
                            overlap_fraction, overlap_length = self._calculate_primer_overlap(qstart, qend, primer_length)
                            
                            # Filter by identity and primer overlap
                            if pident >= identity_threshold and overlap_fraction >= min_primer_overlap:
                                if qseqid not in hits:  # Group by read (query)
                                    hits[qseqid] = []
                                # Determine if this is a forward or reverse primer
                                # Consider both primer name AND strand orientation for reverse complement reads
                                is_reverse_primer_by_name = sseqid in reverse_primers
                                
                                # For reverse complement reads, we need to consider the actual strand
                                # If the primer matches on the minus strand, it's actually a reverse complement match
                                is_minus_strand_match = sstrand == 'minus'
                                
                                # A primer is "forward" if:
                                # 1. It's a forward primer by name AND matches on plus strand, OR
                                # 2. It's a reverse primer by name AND matches on minus strand (reverse complement)
                                is_forward = (not is_reverse_primer_by_name and not is_minus_strand_match) or \
                                           (is_reverse_primer_by_name and is_minus_strand_match)
                                
                                is_reverse_primer = not is_forward
                                
                                hits[qseqid].append({
                                    'primer': sseqid,
                                    'start': read_start,  # Use read coordinates (query)
                                    'end': read_end,
                                    'identity': pident,
                                    'overlap_fraction': overlap_fraction,
                                    'overlap_length': overlap_length,
                                    'evalue': evalue,
                                    'strand': sstrand,
                                    'is_forward': is_forward,  # Based on primer name, not strand
                                    'is_reverse_primer': is_reverse_primer
                                })
                                if self.verbose:
                                    logger.debug(f"Added hit for read {qseqid} at position {read_start}-{read_end} (strand: {sstrand}, overlap: {overlap_fraction:.2f})")
                                    logger.debug(f"  Primer: {sseqid}, is_reverse_by_name: {is_reverse_primer_by_name}, is_minus_strand: {is_minus_strand_match}, is_forward: {is_forward}")
                            else:
                                if self.verbose:
                                    logger.debug(f"Hit filtered out: identity={pident}% (threshold={identity_threshold}%), overlap={overlap_fraction:.2f} (threshold={min_primer_overlap})")
                        else:
                            logger.warning(f"Primer {sseqid} not found in primer file")
                else:
                    logger.warning(f"Invalid line {line_count} in BLAST output: {line.strip()}")
        
        # Sort hits by position for each read
        for read_id in hits:
            hits[read_id].sort(key=lambda x: x['start'])
            if self.verbose:
                logger.info(f"Read {read_id} has {len(hits[read_id])} hits: {[(h['start'], h['end']) for h in hits[read_id]]}")
        
        logger.info(f"Total reads with hits: {len(hits)}")
        return hits
    
    def _split_read_at_primer(self, record, primer_hits, primer_name):
        """Split a read at primer positions and return fragments."""
        fragments = []
        seq = str(record.seq)
        read_len = len(seq)
        
        logger.debug(f"Splitting read {record.id} (length: {read_len}) with {len(primer_hits)} hits")
        
        # Get quality scores if available
        qual = None
        if hasattr(record, 'letter_annotations') and 'phred_quality' in record.letter_annotations:
            qual = record.letter_annotations['phred_quality']
        
        if not primer_hits:
            if self.verbose:
                logger.warning(f"No primer hits found for read {record.id}")
            return fragments
        
        # Sort hits by position to process them in order
        sorted_hits = sorted(primer_hits, key=lambda x: x['start'])
        if self.verbose:
            logger.debug(f"Sorted hits for {record.id}: {[(h['start'], h['end']) for h in sorted_hits]}")
        
        # Start with the full sequence
        current_seq = seq
        current_qual = qual
        current_offset = 0  # Track position offset as we split
        
        for i, hit in enumerate(sorted_hits):
            # Calculate split position relative to the current sequence
            # hit['start'] and hit['end'] are now the READ coordinates (not primer coordinates)
            if hit['is_forward']:
                # Forward primer hit: split BEFORE the primer to keep primer in output
                split_pos = hit['start'] - current_offset
                if self.verbose:
                    logger.debug(f"Forward primer hit: splitting BEFORE primer at position {split_pos} (read pos {hit['start']}, offset {current_offset})")
            else:
                # Reverse primer hit: split AFTER the primer to keep primer in output
                split_pos = hit['end'] - current_offset + 1
                if self.verbose:
                    logger.debug(f"Reverse primer hit: splitting AFTER primer at position {split_pos} (read pos {hit['end']}, offset {current_offset})")
            
            if self.verbose:
                logger.debug(f"Processing hit {i+1}/{len(sorted_hits)}: start={hit['start']}, end={hit['end']}, split_pos={split_pos}, strand={hit['strand']}, is_forward={hit['is_forward']}")
                logger.debug(f"Current sequence length: {len(current_seq)}, current_offset: {current_offset}")
            
            # Validate split position
            if split_pos < 0 or split_pos > len(current_seq):
                if self.verbose:
                    logger.warning(f"Invalid split position {split_pos} for current sequence (length: {len(current_seq)})")
                continue
            
            # Create left fragment (before split position)
            left_seq = current_seq[:split_pos]
            if len(left_seq) >= self.splitting_params['min_fragment_length']:
                left_qual = current_qual[:split_pos] if current_qual else None
                left_record = self._create_fragment_record(
                    record, left_seq, left_qual, primer_name, 0, hit['start']
                )
                fragments.append(left_record)
                if self.verbose:
                    logger.debug(f"Created left fragment: {left_record.id} (length: {len(left_seq)})")
            else:
                if self.verbose:
                    logger.debug(f"Left fragment too short (length: {len(left_seq)}, min: {self.splitting_params['min_fragment_length']}), skipping")
            
            # Update current sequence to be the right fragment (after split position)
            current_seq = current_seq[split_pos:]
            current_qual = current_qual[split_pos:] if current_qual else None
            current_offset = hit['start'] if hit['is_forward'] else hit['end'] + 1
            
            if self.verbose:
                logger.debug(f"Updated current sequence length: {len(current_seq)}")
        
        # Add the final fragment (everything after the last primer)
        if len(current_seq) >= self.splitting_params['min_fragment_length']:
            final_record = self._create_fragment_record(
                record, current_seq, current_qual, primer_name, 1, current_offset
            )
            fragments.append(final_record)
            if self.verbose:
                logger.debug(f"Created final fragment: {final_record.id} (length: {len(current_seq)})")
        else:
            if self.verbose:
                logger.debug(f"Final fragment too short (length: {len(current_seq)}, min: {self.splitting_params['min_fragment_length']}), skipping")
        
        if self.verbose:
            logger.info(f"Successfully split read {record.id} into {len(fragments)} fragments")
        logger.debug(f"Returning {len(fragments)} fragments for read {record.id}")
        return fragments
    
    def _create_fragment_record(self, original_record, seq_str, qual, primer_name, fragment_side, split_pos):
        """Create a new SeqRecord for a fragment."""
        # Create new ID with primer and order information
        start_pos = split_pos
        end_pos = split_pos + len(seq_str)
        
        # Include primer name and fragment order in the header
        if fragment_side == 0:
            # First fragment (before primer)
            new_id = f"{original_record.id}[{start_pos}:{end_pos}][{primer_name}]"
        else:
            # Subsequent fragments (after primer)
            new_id = f"{original_record.id}[{primer_name}][{start_pos}:{end_pos}]"
        
        # Create new record (always FASTA, no quality scores)
        new_record = SeqRecord(
            seq=Seq(seq_str),
            id=new_id,
            description=""  # Empty description to avoid duplicate information
        )
        
        return new_record
    
    def _count_reads(self, file_path):
        """Count the number of reads in a FASTA/FASTQ file."""
        count = 0
        # Detect format dynamically
        file_format = self._detect_format(file_path)
        for _ in SeqIO.parse(file_path, file_format):
            count += 1
        return count
    
    def process(self):
        """Main processing function."""
        logger.info("Starting nanopore concatemer splitting process")
        
        # Create round_outputs under the output directory parent
        output_parent = self.output_file.parent
        
        # Create round_outputs directory for debugging
        round_outputs_dir = output_parent / "round_outputs"
        if self.save_round_outputs:
            round_outputs_dir.mkdir(exist_ok=True)
            logger.info(f"Round outputs will be saved to: {round_outputs_dir}")
        else:
            logger.info("Round outputs disabled")

        # Create a BLAST-safe temp directory without spaces for intermediate files and DBs
        # Prefer user-provided blast_db_dir when available; otherwise use /tmp
        safe_base = self.blast_db_dir if self.blast_db_dir else Path('/tmp')
        temp_dir = safe_base / f"tmp_nanopore_splitter"
        temp_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Using BLAST temp directory: {temp_dir}")
        
        try:
            # Update blast_db_dir to use temp directory if not specified
            if self.blast_db_dir == Path.cwd():
                self.blast_db_dir = temp_dir
            
            # Copy input file to BLAST temp directory and convert to FASTA if needed
            current_input = temp_dir / f"current_input.fasta"
            if self.input_format == 'fastq':
                # Convert FASTQ to FASTA
                with open(current_input, 'w') as f:
                    for record in SeqIO.parse(self.input_file, 'fastq'):
                        SeqIO.write(record, f, 'fasta')
                logger.info(f"Converted FASTQ input to FASTA: {current_input}")
            else:
                # Copy FASTA file
                shutil.copy2(self.input_file, current_input)
                logger.info(f"Copied FASTA input: {current_input}")
            
            # Save initial input to round_outputs for debugging
            initial_round_file = round_outputs_dir / "round_0_initial.fasta"
            shutil.copy2(current_input, initial_round_file)
            logger.info(f"Saved initial reads to: {initial_round_file}")
            
            # Count initial reads
            initial_count = self._count_reads(current_input)
            logger.info(f"Initial read count: {initial_count}")
            
            # Process each primer
            primer_records = list(SeqIO.parse(self.primer_file, 'fasta'))
            total_splits = 0
            
            for i, primer_record in enumerate(primer_records):
                primer_name = primer_record.id
                primer_seq = str(primer_record.seq)
                
                logger.info(f"Processing primer {i+1}/{len(primer_records)}: {primer_name}")
                logger.info(f"Primer sequence: {primer_seq}")
                
                # Save current input before processing this primer
                before_primer_file = round_outputs_dir / f"round_{i+1}_before_{primer_name}.fasta"
                shutil.copy2(current_input, before_primer_file)
                logger.info(f"Saved reads before {primer_name} processing to: {before_primer_file}")
                
                # Create temporary primer file
                temp_primer_file = temp_dir / f"primer_{i}.fasta"
                with open(temp_primer_file, 'w') as f:
                    f.write(f">{primer_name}\n{primer_seq}\n")
                
                # Run BLAST search
                blast_output = temp_dir / f"blast_output_{i}.txt"
                if self.use_reversed_blast:
                    # Use reversed order (reads as DB, primer as query)
                    logger.info(f"Using reversed BLAST approach for {primer_name}")
                    self._run_blast_reversed(current_input, temp_primer_file, blast_output, primer_name)
                    hits = self._parse_blast_results(blast_output, reversed_order=True)
                else:
                    # Use original order (primer as DB, reads as query)
                    logger.info(f"Using original BLAST approach for {primer_name}")
                    db_name = f"primer_db_{i}"
                    db_path = self._create_blast_db(temp_primer_file, db_name)
                    self._run_blast(current_input, db_path, blast_output, primer_name)
                    hits = self._parse_blast_results(blast_output, reversed_order=False)
                
                # Count reads with hits
                reads_with_hits = len(hits)
                logger.info(f"Reads with hits to {primer_name}: {reads_with_hits}")
                
                # Save list of reads that will be split
                if reads_with_hits > 0:
                    reads_to_split_file = round_outputs_dir / f"round_{i+1}_reads_to_split_{primer_name}.txt"
                    with open(reads_to_split_file, 'w') as f:
                        f.write(f"Reads to be split by primer {primer_name}:\n")
                        f.write("=" * 50 + "\n")
                        for read_id in hits:
                            f.write(f"{read_id}\n")
                            for hit in hits[read_id]:
                                f.write(f"  Hit: {hit['start']}-{hit['end']}, identity={hit['identity']}%, overlap={hit['overlap_fraction']:.2f}\n")
                    logger.info(f"Saved list of reads to split to: {reads_to_split_file}")
                
                if reads_with_hits == 0:
                    logger.info(f"No hits found for primer {primer_name}, continuing...")
                    # Save unchanged file for this round
                    after_primer_file = round_outputs_dir / f"round_{i+1}_after_{primer_name}_no_hits.fasta"
                    shutil.copy2(current_input, after_primer_file)
                    logger.info(f"Saved unchanged reads to: {after_primer_file}")
                    continue
                
                # Split reads and create new file
                new_records = []
                splits_this_round = 0
                split_details = []
                
                logger.info(f"Processing {len(list(SeqIO.parse(current_input, 'fasta')))} reads for splitting")
                if self.verbose:
                    logger.info(f"Hits dictionary keys: {list(hits.keys())}")
                
                for record in SeqIO.parse(current_input, 'fasta'):
                    if self.verbose:
                        logger.debug(f"Processing read: {record.id} (length: {len(record.seq)})")
                    
                    if record.id in hits:
                        # Split this read
                        if self.verbose:
                            logger.info(f"Splitting read {record.id} with {len(hits[record.id])} hits")
                            logger.info(f"Hit details for {record.id}: {hits[record.id]}")
                        
                        fragments = self._split_read_at_primer(record, hits[record.id], primer_name)
                        
                        if self.verbose:
                            logger.info(f"Created {len(fragments)} fragments for read {record.id}")
                        
                        if len(fragments) > 0:
                            new_records.extend(fragments)
                            splits_this_round += 1
                            split_details.append({
                                'original_id': record.id,
                                'original_length': len(record.seq),
                                'fragments_created': len(fragments),
                                'fragment_ids': [f.id for f in fragments]
                            })
                            if self.verbose:
                                logger.debug(f"Split read {record.id} into {len(fragments)} fragments")
                        else:
                            if self.verbose:
                                logger.warning(f"No fragments created for read {record.id}, keeping original")
                            new_records.append(record)
                    else:
                        # Keep read unchanged
                        if self.verbose:
                            logger.debug(f"Read {record.id} has no hits, keeping unchanged")
                        new_records.append(record)
                
                # Save splitting details
                split_details_file = round_outputs_dir / f"round_{i+1}_split_details_{primer_name}.txt"
                with open(split_details_file, 'w') as f:
                    f.write(f"Splitting details for primer {primer_name}:\n")
                    f.write("=" * 50 + "\n")
                    f.write(f"Total reads processed: {len(new_records)}\n")
                    f.write(f"Reads split this round: {splits_this_round}\n\n")
                    for detail in split_details:
                        f.write(f"Original read: {detail['original_id']} (length: {detail['original_length']})\n")
                        f.write(f"  Fragments created: {detail['fragments_created']}\n")
                        for frag_id in detail['fragment_ids']:
                            f.write(f"    {frag_id}\n")
                        f.write("\n")
                logger.info(f"Saved splitting details to: {split_details_file}")
                
                logger.info(f"Total reads processed: {len(new_records)}")
                logger.info(f"Reads split this round: {splits_this_round}")
                
                # Write new file (always FASTA)
                new_input = temp_dir / f"round_{i+1}_output.fasta"
                SeqIO.write(new_records, new_input, 'fasta')
                
                # Save after-primer file for debugging
                after_primer_file = round_outputs_dir / f"round_{i+1}_after_{primer_name}.fasta"
                shutil.copy2(new_input, after_primer_file)
                logger.info(f"Saved reads after {primer_name} processing to: {after_primer_file}")
                
                # Update current input for next round
                current_input = new_input
                
                # Count reads after this round
                new_count = len(new_records)
                logger.info(f"Reads split with {primer_name}: {splits_this_round}")
                logger.info(f"Total reads after {primer_name}: {new_count}")
                
                total_splits += splits_this_round
            
            # Copy final result to output (always FASTA)
            shutil.copy2(current_input, self.output_file)
            
            # Save final round file
            final_round_file = round_outputs_dir / f"round_final_output.fasta"
            shutil.copy2(current_input, final_round_file)
            logger.info(f"Saved final output to: {final_round_file}")
            
            # Final statistics
            final_count = self._count_reads(self.output_file)
            logger.info(f"Final read count: {final_count}")
            logger.info(f"Total reads split: {total_splits}")
            logger.info(f"Net increase in reads: {final_count - initial_count}")
            logger.info(f"Output written to: {self.output_file}")
            logger.info(f"Temp files available in: {temp_dir}")
            logger.info(f"Round output files available in: {round_outputs_dir}")
            
        except Exception as e:
            logger.error(f"Error during processing: {e}")
            logger.info(f"Temp files are preserved in: {temp_dir} for debugging")
            logger.info(f"Round output files are preserved in: {round_outputs_dir} for debugging")
            raise


def main():
    """Main function with command line argument parsing."""
    parser = argparse.ArgumentParser(
        description="Split nanopore sequencing reads at primer/barcode positions",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Presets (lib/default_params.py → ADAPTER_SPLITTER):
  pre_demux   Step 1A — split on barcodes before demultiplexing
  post_demux  Step 2A — split on primers / Y-adapters after demultiplexing

Examples:
  python adapter_splitter.py reads.fastq barcodes.fasta out.fasta --preset pre_demux
  python adapter_splitter.py barcode01.fasta primers_and_y.fasta out.fasta --preset post_demux
        """
    )
    
    parser.add_argument('input_file', help='Input FASTA/FASTQ file')
    parser.add_argument('primer_file', help='Reference sequences in FASTA format (barcodes or primers)')
    parser.add_argument('output_file', help='Output file (always FASTA format)')
    parser.add_argument('--preset', choices=list(ADAPTER_SPLITTER.keys()), default='pre_demux',
                        help='Parameter preset (default: pre_demux)')
    parser.add_argument('--blast-db-dir', help='Directory for BLAST databases (default: current directory)')
    parser.add_argument('--identity', type=int, help='BLAST identity percentage (overrides preset)')
    parser.add_argument('--coverage', type=int, help='BLAST query coverage percentage (overrides preset)')
    parser.add_argument('--evalue', type=float, help='BLAST E-value (overrides preset)')
    parser.add_argument('--word-size', type=int, help='BLAST word size (overrides preset)')
    parser.add_argument('--dust', type=str, help='BLAST dust filtering (overrides preset)')
    parser.add_argument('--soft-masking', type=str, help='BLAST soft masking (overrides preset)')
    parser.add_argument('--gapopen', type=int, help='BLAST gap opening penalty (overrides preset)')
    parser.add_argument('--gapextend', type=int, help='BLAST gap extension penalty (overrides preset)')
    parser.add_argument('--penalty', type=int, help='BLAST mismatch penalty (overrides preset)')
    parser.add_argument('--reward', type=int, help='BLAST match reward (overrides preset)')
    parser.add_argument('--max-target-seqs', type=int, help='BLAST max target sequences (overrides preset)')
    parser.add_argument('--task', type=str, help='BLAST task (overrides preset)')
    parser.add_argument('--identity-threshold', type=float, help='Minimum identity for parsing (overrides preset)')
    parser.add_argument('--min-primer-overlap', type=float, help='Minimum primer overlap fraction (overrides preset)')
    parser.add_argument('--min-fragment-length', type=int, help='Minimum fragment length to keep (overrides preset)')
    parser.add_argument('--split-before-primer', action='store_true', help='Force split-before-primer on')
    parser.add_argument('--save-round-outputs', action='store_true', help='Force saving round outputs')
    parser.add_argument('--debug', action='store_true', help='Enable debug logging')
    parser.add_argument('--verbose', action='store_true', help='Print detailed output for each split')
    parser.add_argument('--save-blast-outputs', action='store_true', help='Save BLAST output files for debugging')
    parser.add_argument('--use-reversed-blast', action='store_true', help='Use reversed BLAST (reads as DB, primer as query)')
    
    args = parser.parse_args()
    
    # Set up logging level
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
        logger.setLevel(logging.DEBUG)
    
    try:
        splitter = NanoporeConcatemerSplitter(
            input_file=args.input_file,
            primer_file=args.primer_file,
            output_file=args.output_file,
            blast_db_dir=args.blast_db_dir,
            save_blast_outputs=args.save_blast_outputs or None,
            use_reversed_blast=args.use_reversed_blast or None,
            verbose=args.verbose or None,
            identity=args.identity,
            coverage=args.coverage,
            evalue=args.evalue,
            word_size=args.word_size,
            dust=args.dust,
            soft_masking=args.soft_masking,
            gapopen=args.gapopen,
            gapextend=args.gapextend,
            penalty=args.penalty,
            reward=args.reward,
            max_target_seqs=args.max_target_seqs,
            task=args.task,
            identity_threshold=args.identity_threshold,
            min_primer_overlap=args.min_primer_overlap,
            min_fragment_length=args.min_fragment_length,
            split_before_primer=args.split_before_primer or None,
            save_round_outputs=args.save_round_outputs or None,
            preset=args.preset,
        )
        
        splitter.process()
        
    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
