#!/usr/bin/env python3
"""
ARG filter and plot (Step 2C).

Filters nanopore reads by length, maps them against ARG reference sequences with
BLAST, and writes start/middle/end FASTA bins by hit position on each read.

Code autoannotated using Composer v2.5.

Author: O'Brien, M.E., Fuhrmeister, E.R., Marchand, J.A.
Date: 2026
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

from lib.default_params import ARG_FILTER_AND_PLOT

BLAST_PARAMS = ARG_FILTER_AND_PLOT['BLAST_PARAMS']
PARSING_PARAMS = ARG_FILTER_AND_PLOT['PARSING_PARAMS']
FILTER_PARAMS = ARG_FILTER_AND_PLOT['FILTER_PARAMS']
DEBUG_PARAMS = ARG_FILTER_AND_PLOT['DEBUG_PARAMS']

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class ReadFilterAndSplitter:
    def __init__(self, input_file, primer_file, output_dir, min_length=None, max_length=None, 
                 save_blast_outputs=None, verbose=None, target_filter=None):
        """
        Initialize the filter and splitter with input files and parameters.

        Args:
            input_file: Input FASTA/FASTQ path
            primer_file: ARG reference FASTA path
            output_dir: Output directory for start/middle/end bins
            min_length: Minimum read length (default from FILTER_PARAMS)
            max_length: Maximum read length (default from FILTER_PARAMS)
            save_blast_outputs: Save BLAST output files for debugging
            verbose: Print detailed output
            target_filter: Optional ARG name substring filter
        """
        if input_file is None or primer_file is None or output_dir is None:
            if len(sys.argv) < 4:
                print("Usage: python3 arg_filter_and_plot.py <input_file> <primer_file> <output_dir>")
                sys.exit(1)
            input_file = sys.argv[1]
            primer_file = sys.argv[2]
            output_dir = sys.argv[3]
            
        self.input_file = Path(input_file)
        self.primer_file = Path(primer_file)
        self.output_dir = Path(output_dir)
        
        # Use provided parameters or defaults
        self.min_length = min_length if min_length is not None else FILTER_PARAMS['min_length']
        self.max_length = max_length if max_length is not None else FILTER_PARAMS['max_length']
        self.save_blast_outputs = save_blast_outputs if save_blast_outputs is not None else DEBUG_PARAMS['save_blast_outputs']
        self.verbose = verbose if verbose is not None else DEBUG_PARAMS['verbose']
        self.target_filter = target_filter
        
        # Use parameters from the top of the script
        self.blast_params = BLAST_PARAMS.copy()
        self.parsing_params = PARSING_PARAMS.copy()
        
        # Validate input files and create output directory
        self._validate_inputs()
        self.output_dir.mkdir(exist_ok=True)
        
    def _validate_inputs(self):
        """Validate that input files exist and are readable."""
        if not self.input_file.exists():
            raise FileNotFoundError(f"Input file not found: {self.input_file}")
        if not self.primer_file.exists():
            raise FileNotFoundError(f"Primer file not found: {self.primer_file}")
        
        # Check if input is FASTA
        self.input_format = self._detect_format(self.input_file)
        if self.input_format != 'fasta':
            raise ValueError(f"Input file must be in FASTA format, detected: {self.input_format}")
        
        logger.info(f"Input format: {self.input_format}")
        
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
        """Create a BLAST database from primer sequences."""
        db_base_path = self.output_dir / db_name
        
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
            raise
    
    def _run_blast(self, query_file, db_path, output_file, primer_name=None):
        """Run BLAST search."""
        cmd = [
            'blastn',
            '-query', str(query_file),
            '-db', str(db_path),
            '-out', str(output_file),
            '-outfmt', '10 qseqid sseqid pident qcovs qstart qend sstart send evalue bitscore sstrand'
        ]
        
        logger.info(f"Running BLAST search: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            logger.info(f"BLAST search completed: {output_file}")
            
            # Save BLAST output to local directory if requested
            if self.save_blast_outputs and primer_name:
                blast_outputs_dir = self.output_dir / "blast_outputs"
                blast_outputs_dir.mkdir(exist_ok=True)
                
                local_blast_output = blast_outputs_dir / f"blast_{primer_name}_output.txt"
                shutil.copy2(output_file, local_blast_output)
                logger.info(f"Saved BLAST output to: {local_blast_output}")
                
        except subprocess.CalledProcessError as e:
            logger.error(f"BLAST search failed: {e}")
            raise
    
    def _parse_blast_results(self, blast_output):
        """Parse BLAST results and identify reads with primers at start/end."""
        hits = {}
        
        if not os.path.exists(blast_output) or os.path.getsize(blast_output) == 0:
            logger.warning(f"BLAST output file is empty or doesn't exist: {blast_output}")
            return hits
        
        logger.info(f"Parsing BLAST results from: {blast_output}")
        
        identity_threshold = self.parsing_params['identity_threshold']
        min_primer_overlap = self.parsing_params['min_primer_overlap']
        max_distance_from_end = self.parsing_params['max_distance_from_end']
        
        # Get primer lengths for overlap calculation
        primer_lengths = {}
        for record in SeqIO.parse(self.primer_file, 'fasta'):
            primer_lengths[record.id] = len(record.seq)
        
        logger.info(f"Primer lengths: {primer_lengths}")
        
        with open(blast_output, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                # Parse CSV format
                fields = line.split(',')
                if len(fields) >= 10:
                    qseqid = fields[0]  # read name
                    sseqid = fields[1]  # primer name
                    pident = float(fields[2])
                    qcovs = float(fields[3])
                    qstart = int(fields[4])  # read start
                    qend = int(fields[5])    # read end
                    sstart = int(fields[6])  # primer start
                    send = int(fields[7])    # primer end
                    evalue = float(fields[8])
                    sstrand = fields[10] if len(fields) > 10 else 'plus'
                    
                    # Calculate primer overlap
                    if sseqid in primer_lengths:
                        primer_length = primer_lengths[sseqid]
                        overlap_length = abs(qend - qstart) + 1
                        overlap_fraction = overlap_length / primer_length
                        
                        # Filter by identity and primer overlap
                        if pident >= identity_threshold and overlap_fraction >= min_primer_overlap:
                            # Determine if this is a start or end match
                            is_start_match = qstart <= max_distance_from_end
                            # For end matches, we'll determine this later when we have the actual read length
                            is_end_match = False  # Will be determined in _categorize_reads_by_primers
                            
                            if qseqid not in hits:
                                hits[qseqid] = []
                            
                            hits[qseqid].append({
                                'primer': sseqid,
                                'start': qstart,
                                'end': qend,
                                'identity': pident,
                                'overlap_fraction': overlap_fraction,
                                'evalue': evalue,
                                'strand': sstrand,
                                'is_start_match': is_start_match,
                                'is_end_match': is_end_match
                            })
                            
                            if self.verbose:
                                logger.debug(f"Hit for read {qseqid}: primer={sseqid}, pos={qstart}-{qend}, "
                                           f"identity={pident}%, overlap={overlap_fraction:.2f}")
        
        # Sort hits by position for each read
        for read_id in hits:
            hits[read_id].sort(key=lambda x: x['start'])
        
        logger.info(f"Total reads with primer hits: {len(hits)}")
        return hits
    
    def _filter_reads_by_length(self, input_file):
        """Filter reads by length criteria."""
        filtered_records = []
        total_reads = 0
        passed_filter = 0
        
        logger.info(f"Filtering reads by length criteria: min={self.min_length}, max={self.max_length}")
        
        for record in SeqIO.parse(input_file, 'fasta'):
            total_reads += 1
            read_length = len(record.seq)
            
            # Apply length filters
            if self.min_length and read_length < self.min_length:
                if self.verbose:
                    logger.debug(f"Read {record.id} too short: {read_length} < {self.min_length}")
                continue
            
            if self.max_length and read_length > self.max_length:
                if self.verbose:
                    logger.debug(f"Read {record.id} too long: {read_length} > {self.max_length}")
                continue
            
            filtered_records.append(record)
            passed_filter += 1
            
            if self.verbose:
                logger.debug(f"Read {record.id} passed filter: length={read_length}")
        
        logger.info(f"Length filtering results: {passed_filter}/{total_reads} reads passed filter")
        return filtered_records
    
    def _categorize_reads_by_primers(self, reads, blast_hits):
        """Categorize reads based on primer matches at start/end."""
        categories = {
            'start_matches': {},  # primer_name -> list of reads
            'end_matches': {},    # primer_name -> list of reads
            'both_matches': {},   # primer_name -> list of reads
            'no_matches': []      # reads with no primer matches
        }
        
        # Get all primer names
        primer_names = set()
        for record in SeqIO.parse(self.primer_file, 'fasta'):
            primer_names.add(record.id)
        
        # Initialize categories for each primer
        for primer_name in primer_names:
            categories['start_matches'][primer_name] = []
            categories['end_matches'][primer_name] = []
            categories['both_matches'][primer_name] = []
        
        # Categorize each read
        for record in reads:
            read_id = record.id
            read_length = len(record.seq)
            
            if read_id not in blast_hits:
                categories['no_matches'].append(record)
                continue
            
            # Check hits for this read
            start_primers = set()
            end_primers = set()
            
            for hit in blast_hits[read_id]:
                primer_name = hit['primer']
                
                # Determine if this is a start or end match
                is_start = hit['start'] <= self.parsing_params['max_distance_from_end']
                is_end = hit['end'] >= (read_length - self.parsing_params['max_distance_from_end'])
                
                if is_start:
                    start_primers.add(primer_name)
                if is_end:
                    end_primers.add(primer_name)
            
            # Categorize the read
            if start_primers and end_primers:
                # Read has both start and end matches
                for primer_name in start_primers & end_primers:
                    categories['both_matches'][primer_name].append(record)
            elif start_primers:
                # Read has only start matches
                for primer_name in start_primers:
                    categories['start_matches'][primer_name].append(record)
            elif end_primers:
                # Read has only end matches
                for primer_name in end_primers:
                    categories['end_matches'][primer_name].append(record)
            else:
                # Read has no valid start/end matches
                categories['no_matches'].append(record)
        
        return categories
    
    def _write_categorized_reads(self, categories):
        """Write categorized reads to separate FASTA files."""
        output_files = {}
        
        # Write start matches
        for primer_name, reads in categories['start_matches'].items():
            if reads:
                filename = f"start_{primer_name}.fasta"
                filepath = self.output_dir / filename
                SeqIO.write(reads, filepath, 'fasta')
                output_files[f"start_{primer_name}"] = filepath
                logger.info(f"Wrote {len(reads)} reads with {primer_name} at start to: {filepath}")
        
        # Write end matches
        for primer_name, reads in categories['end_matches'].items():
            if reads:
                filename = f"end_{primer_name}.fasta"
                filepath = self.output_dir / filename
                SeqIO.write(reads, filepath, 'fasta')
                output_files[f"end_{primer_name}"] = filepath
                logger.info(f"Wrote {len(reads)} reads with {primer_name} at end to: {filepath}")
        
        # Write both matches
        for primer_name, reads in categories['both_matches'].items():
            if reads:
                filename = f"both_{primer_name}.fasta"
                filepath = self.output_dir / filename
                SeqIO.write(reads, filepath, 'fasta')
                output_files[f"both_{primer_name}"] = filepath
                logger.info(f"Wrote {len(reads)} reads with {primer_name} at both ends to: {filepath}")
        
        # Write no matches
        if categories['no_matches']:
            filepath = self.output_dir / "no_primer_matches.fasta"
            SeqIO.write(categories['no_matches'], filepath, 'fasta')
            output_files["no_matches"] = filepath
            logger.info(f"Wrote {len(categories['no_matches'])} reads with no primer matches to: {filepath}")
        
        return output_files
    
    def process(self):
        """Main processing function."""
        logger.info("Starting read filtering and splitting process")
        
        # Create temp directory
        temp_dir = self.output_dir / "tmp"
        temp_dir.mkdir(exist_ok=True)
        
        try:
            # Step 1: Filter reads by length
            logger.info("Step 1: Filtering reads by length criteria")
            filtered_reads = self._filter_reads_by_length(self.input_file)
            
            if not filtered_reads:
                logger.warning("No reads passed the length filter")
                return {}
            
            # Write filtered reads to temp file
            filtered_file = temp_dir / "filtered_reads.fasta"
            SeqIO.write(filtered_reads, filtered_file, 'fasta')
            logger.info(f"Wrote {len(filtered_reads)} filtered reads to: {filtered_file}")
            
            # Step 2: Run BLAST search for each primer
            logger.info("Step 2: Running BLAST searches for primer identification")
            all_hits = {}
            
            # Filter primers by target if specified
            target_filter = getattr(self, 'target_filter', None)
            
            for i, primer_record in enumerate(SeqIO.parse(self.primer_file, 'fasta')):
                primer_name = primer_record.id
                primer_seq = str(primer_record.seq)
                
                # Skip primers that don't match the target filter
                if target_filter and target_filter not in primer_name:
                    logger.info(f"Skipping primer {primer_name} (doesn't match target filter: {target_filter})")
                    continue
                
                logger.info(f"Processing primer {i+1}: {primer_name}")
                
                # Create temporary primer file
                temp_primer_file = temp_dir / f"primer_{i}.fasta"
                with open(temp_primer_file, 'w') as f:
                    f.write(f">{primer_name}\n{primer_seq}\n")
                
                # Create BLAST database and run search
                db_name = f"primer_db_{i}"
                db_path = self._create_blast_db(temp_primer_file, db_name)
                
                blast_output = temp_dir / f"blast_output_{i}.txt"
                self._run_blast(filtered_file, db_path, blast_output, primer_name)
                
                # Parse BLAST results
                hits = self._parse_blast_results(blast_output)
                all_hits.update(hits)
            
            # Step 3: Categorize reads by primer matches
            logger.info("Step 3: Categorizing reads by primer matches")
            categories = self._categorize_reads_by_primers(filtered_reads, all_hits)
            
            # Step 4: Write categorized reads to files
            logger.info("Step 4: Writing categorized reads to output files")
            output_files = self._write_categorized_reads(categories)
            
            # Generate summary report
            self._write_summary_report(categories, output_files)
            
            logger.info("Processing completed successfully")
            return output_files
            
        except Exception as e:
            logger.error(f"Error during processing: {e}")
            raise
        finally:
            # Clean up temp directory
            if not self.save_blast_outputs:
                shutil.rmtree(temp_dir, ignore_errors=True)
                logger.info("Cleaned up temporary files")
    
    def _write_summary_report(self, categories, output_files):
        """Write a summary report of the processing results."""
        report_file = self.output_dir / "processing_summary.txt"
        
        with open(report_file, 'w') as f:
            f.write("Read Filter and Split Processing Summary\n")
            f.write("=" * 50 + "\n\n")
            
            # Input information
            f.write("Input Information:\n")
            f.write(f"  Input file: {self.input_file}\n")
            f.write(f"  Primer file: {self.primer_file}\n")
            f.write(f"  Output directory: {self.output_dir}\n\n")
            
            # Filter criteria
            f.write("Filter Criteria:\n")
            f.write(f"  Minimum length: {self.min_length if self.min_length else 'None'}\n")
            f.write(f"  Maximum length: {self.max_length if self.max_length else 'None'}\n\n")
            
            # Results summary
            f.write("Results Summary:\n")
            
            # Start matches
            total_start = sum(len(reads) for reads in categories['start_matches'].values())
            f.write(f"  Reads with primer at start: {total_start}\n")
            for primer_name, reads in categories['start_matches'].items():
                if reads:
                    f.write(f"    {primer_name}: {len(reads)} reads\n")
            
            # End matches
            total_end = sum(len(reads) for reads in categories['end_matches'].values())
            f.write(f"  Reads with primer at end: {total_end}\n")
            for primer_name, reads in categories['end_matches'].items():
                if reads:
                    f.write(f"    {primer_name}: {len(reads)} reads\n")
            
            # Both matches
            total_both = sum(len(reads) for reads in categories['both_matches'].values())
            f.write(f"  Reads with primer at both ends: {total_both}\n")
            for primer_name, reads in categories['both_matches'].items():
                if reads:
                    f.write(f"    {primer_name}: {len(reads)} reads\n")
            
            # No matches
            f.write(f"  Reads with no primer matches: {len(categories['no_matches'])}\n\n")
            
            # Output files
            f.write("Output Files:\n")
            for category, filepath in output_files.items():
                f.write(f"  {category}: {filepath}\n")
        
        logger.info(f"Summary report written to: {report_file}")


def main():
    """Main function with command line argument parsing."""
    parser = argparse.ArgumentParser(
        description="Filter reads by length and split by primer sequences at start/end",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python read_filter_and_split.py input.fasta primers.fasta output_dir
  python read_filter_and_split.py input.fasta primers.fasta output_dir --min-length 1000 --max-length 5000
  python read_filter_and_split.py input.fasta primers.fasta output_dir --verbose --save-blast-outputs
  python read_filter_and_split.py input.fasta primers.fasta output_dir --target ctx-m
        """
    )
    
    parser.add_argument('input_file', help='Input FASTA file')
    parser.add_argument('primer_file', help='Primer sequences in FASTA format')
    parser.add_argument('output_dir', help='Output directory for split files')
    parser.add_argument('--min-length', type=int, help='Minimum read length filter')
    parser.add_argument('--max-length', type=int, help='Maximum read length filter')
    parser.add_argument('--identity', type=int, default=80, help='Minimum identity percentage for BLAST hits (default: 80)')
    parser.add_argument('--coverage', type=int, default=90, help='Minimum query coverage percentage (default: 90)')
    parser.add_argument('--evalue', type=float, default=1e-5, help='Maximum E-value (default: 1e-5)')
    parser.add_argument('--max-distance', type=int, default=50, help='Maximum distance from read start/end to consider primer match (default: 50)')
    parser.add_argument('--debug', action='store_true', help='Enable debug logging')
    parser.add_argument('--verbose', action='store_true', help='Print detailed output')
    parser.add_argument('--save-blast-outputs', action='store_true', help='Save BLAST output files for debugging')
    parser.add_argument('--target', type=str, help='Specific target to process (e.g., ctx-m, KPC, OXA, qnrS). If not specified, processes all targets.')
    
    args = parser.parse_args()
    
    # Set up logging level
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
        logger.setLevel(logging.DEBUG)
    
    try:
        # Create filter and splitter instance
        filter_splitter = ReadFilterAndSplitter(
            input_file=args.input_file,
            primer_file=args.primer_file,
            output_dir=args.output_dir,
            min_length=args.min_length,
            max_length=args.max_length,
            save_blast_outputs=args.save_blast_outputs,
            verbose=args.verbose,
            target_filter=args.target
        )
        
        # Update parameters if provided via command line
        if args.identity != 80:
            filter_splitter.parsing_params['identity_threshold'] = args.identity
        if args.coverage != 90:
            filter_splitter.blast_params['query_coverage'] = args.coverage
        if args.evalue != 1e-5:
            filter_splitter.blast_params['evalue'] = args.evalue
        if args.max_distance != 50:
            filter_splitter.parsing_params['max_distance_from_end'] = args.max_distance
        
        # Process the data
        output_files = filter_splitter.process()
        
        logger.info("Processing completed successfully!")
        logger.info(f"Output files created in: {args.output_dir}")
        
    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
