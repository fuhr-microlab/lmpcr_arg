#!/usr/bin/env python3
"""
Adapter/ARG Primer Filter

Code autoannotated using Composer v2.5.

This script filters nanopore sequencing reads by removing reads with multiple hits
to primer sequences. Useful for removing chimeric reads or reads with multiple ARG hits.

Author: O'Brien, M.E., Fuhrmeister, E.R., Marchand, J.A.
Date: 2026

FIX (260630): Each invocation now creates a unique temporary directory via
tempfile.mkdtemp() instead of writing to the shared /tmp/tmp_adapter_filter path.
Running multiple barcodes in parallel on the same node previously caused BLAST
inputs/outputs from one barcode to overwrite another's, producing impossible summary
statistics (hits >> initial reads, negative "reads with no hits", etc.).
"""

import os
import sys
import subprocess
import shutil
import tempfile
from pathlib import Path
from Bio import SeqIO
import argparse
import logging

from lib.default_params import ADAPTER_FILTER

BLAST_PARAMS = ADAPTER_FILTER['BLAST_PARAMS']
PARSING_PARAMS = ADAPTER_FILTER['PARSING_PARAMS']
FILTERING_PARAMS = ADAPTER_FILTER['FILTERING_PARAMS']
DEBUG_PARAMS = ADAPTER_FILTER['DEBUG_PARAMS']

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class AdapterFilter:
    def __init__(self, input_file, primer_file, output_file, blast_db_dir=None, 
                 save_blast_outputs=None, verbose=None, max_hits=None, min_length=None,
                 identity=None, coverage=None, evalue=None, identity_threshold=None, 
                 min_primer_overlap=None):
        """
        Initialize the filter with input files and parameters.
        
        Args:
            input_file (str): Path to input FASTA/FASTQ file
            primer_file (str): Path to primer FASTA file
            output_file (str): Path to output file
            blast_db_dir (str): Directory for BLAST databases (optional)
            save_blast_outputs (bool): Save BLAST output files for debugging
            verbose (bool): Print detailed output
            max_hits (int): Maximum hits per read (override FILTERING_PARAMS)
            min_length (int): Minimum read length (override FILTERING_PARAMS)
            identity (int): BLAST identity percentage (override BLAST_PARAMS)
            coverage (int): BLAST query coverage (override BLAST_PARAMS)
            evalue (float): BLAST E-value (override BLAST_PARAMS)
            identity_threshold (float): Minimum identity for parsing (override PARSING_PARAMS)
            min_primer_overlap (float): Minimum primer overlap (override PARSING_PARAMS)
        """
        self.input_file = Path(input_file)
        self.primer_file = Path(primer_file)
        self.output_file = Path(output_file)
        self.blast_db_dir = Path(blast_db_dir) if blast_db_dir else Path.cwd()
        
        # Use provided parameters or defaults
        self.save_blast_outputs = save_blast_outputs if save_blast_outputs is not None else DEBUG_PARAMS['save_blast_outputs']
        self.verbose = verbose if verbose is not None else DEBUG_PARAMS['verbose']
        
        # Copy parameter dictionaries
        self.blast_params = BLAST_PARAMS.copy()
        self.parsing_params = PARSING_PARAMS.copy()
        self.filtering_params = FILTERING_PARAMS.copy()
        
        # Override with provided parameters
        if max_hits is not None:
            self.filtering_params['max_hits_per_read'] = max_hits
        if min_length is not None:
            self.filtering_params['min_read_length'] = min_length
        if identity is not None:
            self.blast_params['identity'] = identity
        if coverage is not None:
            self.blast_params['query_coverage'] = coverage
        if evalue is not None:
            self.blast_params['evalue'] = evalue
        if identity_threshold is not None:
            self.parsing_params['identity_threshold'] = identity_threshold
        if min_primer_overlap is not None:
            self.parsing_params['min_primer_overlap'] = min_primer_overlap
        
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
    
    def _create_blast_db(self, primer_file, db_name, db_dir):
        """Create a BLAST database from primer sequences."""
        db_base_path = db_dir / db_name
        
        logger.info(f"Creating BLAST database: {db_base_path}")
        
        cmd = [
            'makeblastdb',
            '-in', str(primer_file),
            '-dbtype', 'nucl',
            '-out', str(db_base_path)
        ]
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info(f"Created BLAST database: {db_base_path}")
            return db_base_path
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to create BLAST database: {e}")
            logger.error(f"STDERR: {e.stderr}")
            raise
    
    def _run_blast(self, query_file, db_path, output_file):
        """Run BLAST search."""
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
        
        logger.info(f"Running BLAST search...")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            logger.info(f"BLAST search completed")
            
            if output_file.exists() and output_file.stat().st_size == 0:
                logger.info("BLAST output is empty - no hits found")
            
        except subprocess.CalledProcessError as e:
            logger.error(f"BLAST search failed: {e}")
            raise
    
    def _calculate_primer_overlap(self, hit_start, hit_end, primer_length):
        """Calculate the fraction of primer that overlaps with the hit."""
        if hit_start > hit_end:
            hit_start, hit_end = hit_end, hit_start
        
        overlap_length = hit_end - hit_start + 1
        overlap_fraction = overlap_length / primer_length
        
        return overlap_fraction, overlap_length
    
    def _parse_blast_results(self, blast_output):
        """Parse BLAST results and count hits per read."""
        hits_per_read = {}
        
        if not os.path.exists(blast_output) or os.path.getsize(blast_output) == 0:
            logger.info(f"BLAST output file is empty or doesn't exist: {blast_output}")
            return hits_per_read
        
        logger.info(f"Parsing BLAST results from: {blast_output}")
        
        identity_threshold = self.parsing_params['identity_threshold']
        min_primer_overlap = self.parsing_params['min_primer_overlap']
        
        # Get primer lengths
        primer_lengths = {}
        for record in SeqIO.parse(self.primer_file, 'fasta'):
            primer_lengths[record.id] = len(record.seq)
        
        with open(blast_output, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                # Parse CSV format
                fields = line.split(',') if ',' in line else line.split('\t')
                
                if len(fields) >= 10:
                    qseqid = fields[0]  # read name
                    sseqid = fields[1]  # primer name
                    pident = float(fields[2])
                    qcovs = float(fields[3])
                    qstart = int(fields[4])
                    qend = int(fields[5])
                    sstart = int(fields[6])
                    send = int(fields[7])
                    evalue = float(fields[8])
                    
                    # Calculate primer overlap
                    if sseqid in primer_lengths:
                        primer_length = primer_lengths[sseqid]
                        overlap_fraction, overlap_length = self._calculate_primer_overlap(
                            sstart, send, primer_length
                        )
                        
                        # Filter by identity and primer overlap
                        if pident >= identity_threshold and overlap_fraction >= min_primer_overlap:
                            if qseqid not in hits_per_read:
                                hits_per_read[qseqid] = []
                            
                            hits_per_read[qseqid].append({
                                'primer': sseqid,
                                'start': qstart,
                                'end': qend,
                                'identity': pident,
                                'overlap_fraction': overlap_fraction,
                                'evalue': evalue
                            })
        
        logger.info(f"Total reads with hits: {len(hits_per_read)}")
        return hits_per_read
    
    def process(self):
        """Main processing function."""
        logger.info("Starting adapter filtering process")
        logger.info(f"Input file: {self.input_file}")
        logger.info(f"Max hits per read: {self.filtering_params['max_hits_per_read']}")
        logger.info(f"Min read length: {self.filtering_params['min_read_length']}")
        
        # --- FIX: use a unique per-run temp directory to avoid cross-job collisions ---
        # Previously this used a fixed path (/tmp/tmp_adapter_filter) shared across all
        # concurrent invocations, causing BLAST outputs from one barcode to corrupt
        # another's counts when jobs ran in parallel on the same node.
        temp_dir = Path(tempfile.mkdtemp(prefix="adapter_filter_"))
        logger.info(f"Using isolated BLAST temp directory: {temp_dir}")
        
        try:
            # Copy input file to temp directory and convert to FASTA if needed
            current_input = temp_dir / "input.fasta"
            if self.input_format == 'fastq':
                with open(current_input, 'w') as f:
                    for record in SeqIO.parse(self.input_file, 'fastq'):
                        SeqIO.write(record, f, 'fasta')
                logger.info(f"Converted FASTQ input to FASTA")
            else:
                shutil.copy2(self.input_file, current_input)
                logger.info(f"Copied FASTA input")
            
            # Count initial reads
            initial_count = sum(1 for _ in SeqIO.parse(current_input, 'fasta'))
            logger.info(f"Initial read count: {initial_count}")
            
            # Create BLAST database from primers
            db_name = "primer_db"
            db_path = self._create_blast_db(self.primer_file, db_name, temp_dir)
            
            # Run BLAST search
            blast_output = temp_dir / "blast_output.txt"
            self._run_blast(current_input, db_path, blast_output)
            
            # Parse BLAST results
            hits_per_read = self._parse_blast_results(blast_output)
            
            # Save BLAST output if requested
            if self.save_blast_outputs:
                blast_outputs_dir = self.output_file.parent / f"blast_outputs_{self.output_file.stem}"
                blast_outputs_dir.mkdir(exist_ok=True)
                
                local_blast_output = blast_outputs_dir / "filter_blast_output.txt"
                shutil.copy2(blast_output, local_blast_output)
                logger.info(f"Saved BLAST output to: {local_blast_output}")
                
                # Create summary
                summary_file = blast_outputs_dir / "filter_blast_summary.txt"
                with open(summary_file, 'w') as f:
                    f.write("BLAST Search Summary - Adapter Filter\n")
                    f.write("=" * 50 + "\n\n")
                    f.write(f"Input file: {self.input_file}\n")
                    f.write(f"Primer file: {self.primer_file}\n")
                    f.write(f"Total reads: {initial_count}\n")
                    f.write(f"Reads with hits: {len(hits_per_read)}\n\n")
                    f.write("Hit distribution:\n")
                    hit_counts = {}
                    for read_id, hits in hits_per_read.items():
                        num_hits = len(hits)
                        hit_counts[num_hits] = hit_counts.get(num_hits, 0) + 1
                    for num_hits in sorted(hit_counts.keys()):
                        f.write(f"  {num_hits} hits: {hit_counts[num_hits]} reads\n")
                logger.info(f"Saved BLAST summary to: {summary_file}")
            
            # Filter reads
            max_hits = self.filtering_params['max_hits_per_read']
            min_length = self.filtering_params['min_read_length']
            
            kept_reads = []
            removed_multi_hit = 0
            removed_too_short = 0
            
            logger.info(f"Filtering reads...")
            
            for record in SeqIO.parse(current_input, 'fasta'):
                # Check hit count
                num_hits = len(hits_per_read.get(record.id, []))
                
                if num_hits > max_hits:
                    removed_multi_hit += 1
                    if self.verbose:
                        logger.debug(f"Removing {record.id}: {num_hits} hits > {max_hits}")
                    continue
                
                # Check length
                if len(record.seq) < min_length:
                    removed_too_short += 1
                    if self.verbose:
                        logger.debug(f"Removing {record.id}: length {len(record.seq)} < {min_length}")
                    continue
                
                # Keep this read
                kept_reads.append(record)
                if self.verbose:
                    logger.debug(f"Keeping {record.id}: {num_hits} hits, length {len(record.seq)}")
            
            # Write output
            SeqIO.write(kept_reads, self.output_file, 'fasta')
            
            # Final statistics
            final_count = len(kept_reads)
            total_removed = removed_multi_hit + removed_too_short
            
            logger.info(f"\n" + "=" * 60)
            logger.info(f"FILTERING SUMMARY")
            logger.info(f"=" * 60)
            logger.info(f"Initial reads: {initial_count}")
            logger.info(f"Reads with hits: {len(hits_per_read)}")
            logger.info(f"Removed (>{max_hits} hits): {removed_multi_hit}")
            logger.info(f"Removed (too short): {removed_too_short}")
            logger.info(f"Total removed: {total_removed}")
            logger.info(f"Final reads: {final_count}")
            logger.info(f"Retention rate: {final_count/initial_count*100:.2f}%")
            logger.info(f"=" * 60)
            
            # Show hit distribution
            hit_counts = {}
            for read_id, hits in hits_per_read.items():
                num_hits = len(hits)
                hit_counts[num_hits] = hit_counts.get(num_hits, 0) + 1
            
            if hit_counts:
                logger.info(f"\nHit distribution:")
                for num_hits in sorted(hit_counts.keys()):
                    logger.info(f"  {num_hits} hits: {hit_counts[num_hits]} reads")
            
            logger.info(f"\nOutput written to: {self.output_file}")
            
            # Write summary report to file
            self._write_summary_report(
                initial_count, final_count, removed_multi_hit, removed_too_short,
                hits_per_read, hit_counts, max_hits, min_length
            )
            
        except Exception as e:
            logger.error(f"Error during processing: {e}")
            logger.info(f"Temp files preserved in: {temp_dir} for debugging")
            raise
        finally:
            # Clean up unique temp directory
            try:
                shutil.rmtree(temp_dir)
                logger.info(f"Cleaned up temp directory: {temp_dir}")
            except Exception:
                pass
    
    def _write_summary_report(self, initial_count, final_count, removed_multi_hit, 
                              removed_too_short, hits_per_read, hit_counts, max_hits, min_length):
        """Write a summary report of the filtering results."""
        report_file = self.output_file.parent / f"{self.output_file.stem}_summary.txt"
        
        with open(report_file, 'w') as f:
            f.write("Adapter/ARG Filter Summary\n")
            f.write("=" * 60 + "\n\n")
            
            # Input information
            f.write("Input Information:\n")
            f.write(f"  Input file: {self.input_file}\n")
            f.write(f"  Primer file: {self.primer_file}\n")
            f.write(f"  Output file: {self.output_file}\n\n")
            
            # Filter criteria
            f.write("Filter Criteria:\n")
            f.write(f"  Maximum hits per read: {max_hits}\n")
            f.write(f"  Minimum read length: {min_length if min_length > 0 else 'None'}\n")
            f.write(f"  Identity threshold: {self.parsing_params['identity_threshold']}%\n")
            f.write(f"  Minimum primer overlap: {self.parsing_params['min_primer_overlap']}\n\n")
            
            # Results summary
            f.write("Results Summary:\n")
            f.write(f"  Initial reads: {initial_count}\n")
            f.write(f"  Reads with hits: {len(hits_per_read)}\n")
            f.write(f"  Reads with no hits: {initial_count - len(hits_per_read)}\n")
            f.write(f"  Removed (>{max_hits} hits): {removed_multi_hit}\n")
            f.write(f"  Removed (too short): {removed_too_short}\n")
            f.write(f"  Total removed: {removed_multi_hit + removed_too_short}\n")
            f.write(f"  Final reads: {final_count}\n")
            f.write(f"  Retention rate: {final_count/initial_count*100:.2f}%\n\n")
            
            # Hit distribution
            if hit_counts:
                f.write("Hit Distribution:\n")
                for num_hits in sorted(hit_counts.keys()):
                    count = hit_counts[num_hits]
                    percentage = (count / initial_count * 100) if initial_count > 0 else 0
                    f.write(f"  {num_hits} hit(s): {count} reads ({percentage:.2f}%)\n")
                f.write("\n")
            
            # Multi-hit reads breakdown
            multi_hit_reads = {read_id: hits for read_id, hits in hits_per_read.items() 
                              if len(hits) > max_hits}
            
            if multi_hit_reads:
                f.write(f"Multi-Hit Reads (>{max_hits} hits) - REMOVED:\n")
                f.write(f"  Total multi-hit reads: {len(multi_hit_reads)}\n")
                
                # Count by number of hits
                multi_hit_counts = {}
                for read_id, hits in multi_hit_reads.items():
                    num_hits = len(hits)
                    multi_hit_counts[num_hits] = multi_hit_counts.get(num_hits, 0) + 1
                
                f.write("  Breakdown by hit count:\n")
                for num_hits in sorted(multi_hit_counts.keys()):
                    f.write(f"    {num_hits} hits: {multi_hit_counts[num_hits]} reads\n")
                
                # Count by ARG/primer type
                primer_counts = {}
                for read_id, hits in multi_hit_reads.items():
                    for hit in hits:
                        primer = hit['primer']
                        primer_counts[primer] = primer_counts.get(primer, 0) + 1
                
                if primer_counts:
                    f.write("\n  ARG/Primer types found in multi-hit reads:\n")
                    for primer, count in sorted(primer_counts.items(), key=lambda x: x[1], reverse=True):
                        f.write(f"    {primer}: {count} hits\n")
                f.write("\n")
            
            # Single-hit reads (kept)
            single_hit_reads = {read_id: hits for read_id, hits in hits_per_read.items() 
                               if len(hits) <= max_hits}
            
            if single_hit_reads:
                f.write(f"Single-Hit Reads (≤{max_hits} hits) - KEPT:\n")
                f.write(f"  Total single-hit reads: {len(single_hit_reads)}\n")
                
                # Count by ARG/primer type
                kept_primer_counts = {}
                for read_id, hits in single_hit_reads.items():
                    for hit in hits:
                        primer = hit['primer']
                        kept_primer_counts[primer] = kept_primer_counts.get(primer, 0) + 1
                
                if kept_primer_counts:
                    f.write("  ARG/Primer types in kept reads:\n")
                    for primer, count in sorted(kept_primer_counts.items(), key=lambda x: x[1], reverse=True):
                        f.write(f"    {primer}: {count} reads\n")
        
        logger.info(f"Summary report written to: {report_file}")


def main():
    """Main function with command line argument parsing."""
    parser = argparse.ArgumentParser(
        description="Filter reads by number of primer hits (remove reads with multiple ARG hits)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Remove reads with 2+ hits (default)
  python adapter_filter.py input.fasta arg_primers.fasta output.fasta
  
  # Remove reads with 3+ hits
  python adapter_filter.py input.fasta arg_primers.fasta output.fasta --max-hits 2
  
  # Also filter by minimum length
  python adapter_filter.py input.fasta arg_primers.fasta output.fasta --min-length 1000
  
  # Adjust BLAST sensitivity
  python adapter_filter.py input.fasta arg_primers.fasta output.fasta --identity 85 --coverage 95
        """
    )
    
    parser.add_argument('input_file', help='Input FASTA/FASTQ file')
    parser.add_argument('primer_file', help='Primer sequences in FASTA format (e.g., ARG primers)')
    parser.add_argument('output_file', help='Output FASTA file (reads with <=max_hits)')
    
    # Filtering parameters (defaults from lib/default_params.py → ADAPTER_FILTER)
    parser.add_argument('--max-hits', type=int, default=None,
                        help=f"Maximum number of hits per read (default: {FILTERING_PARAMS['max_hits_per_read']})")
    parser.add_argument('--min-length', type=int, default=None,
                        help=f"Minimum read length to keep (default: {FILTERING_PARAMS['min_read_length']})")
    
    # BLAST parameters
    parser.add_argument('--blast-db-dir', help='Directory for BLAST databases (default: /tmp)')
    parser.add_argument('--identity', type=int,
                        help=f"BLAST identity percentage (default: {BLAST_PARAMS['identity']})")
    parser.add_argument('--coverage', type=int,
                        help=f"BLAST query coverage percentage (default: {BLAST_PARAMS['query_coverage']})")
    parser.add_argument('--evalue', type=float,
                        help=f"BLAST E-value (default: {BLAST_PARAMS['evalue']})")
    parser.add_argument('--identity-threshold', type=float,
                        help=f"Minimum identity for hit acceptance (default: {PARSING_PARAMS['identity_threshold']})")
    parser.add_argument('--min-primer-overlap', type=float,
                        help=f"Minimum primer overlap fraction (default: {PARSING_PARAMS['min_primer_overlap']})")
    
    # Debug options
    parser.add_argument('--debug', action='store_true', help='Enable debug logging')
    parser.add_argument('--verbose', action='store_true', help='Print detailed output for each read')
    parser.add_argument('--save-blast-outputs', action='store_true', 
                        help='Save BLAST output files for debugging')
    
    args = parser.parse_args()
    
    # Set up logging level
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
        logger.setLevel(logging.DEBUG)
    
    try:
        # Create filter instance
        filter_obj = AdapterFilter(
            input_file=args.input_file,
            primer_file=args.primer_file,
            output_file=args.output_file,
            blast_db_dir=args.blast_db_dir,
            save_blast_outputs=args.save_blast_outputs or None,
            verbose=args.verbose or None,
            max_hits=args.max_hits,
            min_length=args.min_length,
            identity=args.identity,
            coverage=args.coverage,
            evalue=args.evalue,
            identity_threshold=args.identity_threshold,
            min_primer_overlap=args.min_primer_overlap
        )
        
        # Process the data
        filter_obj.process()
        
    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
