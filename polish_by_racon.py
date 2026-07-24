#!/usr/bin/env python3
"""
RACON-based polishing module.

This module implements consensus sequence calculation using RACON, which is
optimized for long-read consensus polishing and error correction.
"""

import os
import subprocess
import shutil
import random
import logging
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import Counter

from lib.default_params import RACON_PARAMS

logger = logging.getLogger(__name__)


class RaconPolisher:
    """RACON-based consensus polishing implementation."""
    
    def __init__(self, racon_params=None):
        """
        Initialize RACON polisher.
        
        Args:
            racon_params (dict): RACON parameters (uses defaults if None)
        """
        self.racon_params = racon_params or RACON_PARAMS.copy()
    
    def check_racon_availability(self):
        """Check if RACON is available and get version."""
        try:
            result = subprocess.run(['racon', '--version'], capture_output=True, text=True, check=True)
            logger.info(f"RACON version: {result.stdout.strip()}")
            return True
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            logger.warning(f"RACON not available: {e}")
            return False
    
    def check_minimap2_availability(self):
        """Check if minimap2 is available (required for RACON)."""
        try:
            result = subprocess.run(['minimap2', '--version'], capture_output=True, text=True, check=True)
            logger.info(f"minimap2 version: {result.stdout.strip()}")
            return True
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            logger.warning(f"minimap2 not available: {e}")
            return False
    
    def run_minimap2_alignment(self, reads_file, reference_file, temp_dir, cluster_id):
        """Run minimap2 to align reads to reference."""
        if not self.check_minimap2_availability():
            return None
        
        # Create PAF output file
        paf_file = temp_dir / f"minimap2_alignment_{cluster_id}.paf"
        
        # Run minimap2 with parameters optimized for similar sequences
        cmd = [
            'minimap2',
            '-x', 'map-ont',  # Optimized for Oxford Nanopore reads
            '-N', '50',  # Maximum number of secondary alignments
            '--secondary=no',  # No secondary alignments
            str(reference_file),
            str(reads_file)
        ]
        
        logger.info(f"Running minimap2 alignment for cluster {cluster_id}: {' '.join(cmd)}")
        
        try:
            with open(paf_file, 'w') as f:
                result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True, check=True)
            
            logger.info(f"minimap2 alignment completed for cluster {cluster_id}")
            return paf_file
            
        except subprocess.CalledProcessError as e:
            logger.error(f"minimap2 alignment failed for cluster {cluster_id}: {e}")
            if e.stderr:
                logger.error(f"minimap2 stderr: {e.stderr}")
            return None
        except Exception as e:
            logger.error(f"Unexpected error during minimap2 alignment for cluster {cluster_id}: {e}")
            return None
    
    def run_racon_polishing(self, reads_file, paf_file, reference_file, temp_dir, cluster_id):
        """Run RACON to polish the consensus sequence."""
        if not self.check_racon_availability():
            return None
        
        # Create output file
        consensus_file = temp_dir / f"racon_consensus_{cluster_id}.fasta"
        
        # Run RACON with configurable parameters
        cmd = [
            'racon',
            str(reads_file),  # Reads
            str(paf_file),    # PAF alignment file
            str(reference_file),  # Reference sequence
            '-m', str(self.racon_params['match']),
            '-x', str(self.racon_params['mismatch']),
            '-g', str(self.racon_params['gap']),
            '-w', str(self.racon_params['window_length']),
            '-q', str(self.racon_params['quality_threshold']),
            '-e', str(self.racon_params['error_threshold']),
            '-t', str(self.racon_params['threads'])
        ]
        
        # Add optional flags
        if self.racon_params['include_unpolished']:
            cmd.append('-u')
        if self.racon_params['split']:
            cmd.append('-s')
        if self.racon_params['fragments']:
            cmd.append('-f')
        
        logger.info(f"Running RACON polishing for cluster {cluster_id}: {' '.join(cmd)}")
        
        try:
            with open(consensus_file, 'w') as f:
                result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True, check=True)
            
            logger.info(f"RACON polishing completed for cluster {cluster_id}")
            return consensus_file
            
        except subprocess.CalledProcessError as e:
            logger.error(f"RACON polishing failed for cluster {cluster_id}: {e}")
            if e.stderr:
                logger.error(f"RACON stderr: {e.stderr}")
            return None
        except Exception as e:
            logger.error(f"Unexpected error during RACON polishing for cluster {cluster_id}: {e}")
            return None
    
    def create_simple_consensus(self, cluster_file, cluster_id):
        """Create a simple consensus as fallback when RACON is not available."""
        try:
            reads = list(SeqIO.parse(cluster_file, 'fasta'))
            if len(reads) < 2:
                logger.warning(f"Cluster {cluster_id} has only {len(reads)} sequence(s), cannot calculate consensus")
                return None
            
            logger.info(f"Creating simple consensus for cluster {cluster_id} with {len(reads)} reads")
            
            # Find the longest sequence length
            max_length = max(len(read.seq) for read in reads)
            logger.info(f"Cluster {cluster_id}: Max length: {max_length}")
            
            # Pad shorter sequences to match the longest (add N's instead of gaps)
            padded_reads = []
            for read in reads:
                if len(read.seq) < max_length:
                    # Pad with N's to match longest sequence
                    padding = 'N' * (max_length - len(read.seq))
                    padded_seq = str(read.seq) + padding
                    padded_reads.append(padded_seq)
                else:
                    padded_reads.append(str(read.seq))
            
            # Calculate consensus position by position
            consensus_seq = ""
            for pos in range(max_length):
                bases = [seq[pos] for seq in padded_reads]
                
                # Count N's and gaps
                n_count = sum(1 for base in bases if base.upper() == 'N')
                gap_count = sum(1 for base in bases if base.upper() in ['-', '.'])
                total_gap_like = n_count + gap_count
                gap_fraction = total_gap_like / len(bases)
                
                # Filter out N's and gaps for base counting
                valid_bases = [base for base in bases if base.upper() not in ['N', '-', '.']]
                
                # Check if gap threshold is met
                if gap_fraction >= 0.3:  # 30% gap threshold
                    consensus_seq += '-'
                elif not valid_bases:
                    consensus_seq += 'N'
                else:
                    # Count bases
                    base_counts = Counter(valid_bases)
                    most_common_base, count = base_counts.most_common(1)[0]
                    
                    # Calculate coverage
                    total_sequences = len(bases)
                    coverage = count / total_sequences
                    
                    # Use simple majority voting
                    if coverage >= 0.5:  # 50% majority
                        consensus_seq += most_common_base.upper()
                    else:
                        consensus_seq += most_common_base.upper()
            
            # Trim trailing N's
            consensus_seq = consensus_seq.rstrip('N')
            
            if not consensus_seq:
                logger.warning(f"Cluster {cluster_id} consensus is empty after trimming")
                return None
            
            # Create consensus sequence record
            consensus_record = SeqRecord(
                seq=Seq(consensus_seq),
                id=f"cluster_{cluster_id:04d}_simple_consensus",
                description=f"Simple consensus for cluster {cluster_id} ({len(consensus_seq)} bp)"
            )
            
            logger.info(f"Created simple consensus for cluster {cluster_id}: {len(consensus_seq)} bp")
            return consensus_record
            
        except Exception as e:
            logger.error(f"Error creating simple consensus for cluster {cluster_id}: {e}")
            return None
    
    def polish_cluster(self, cluster_file, cluster_id, temp_dir, output_dir):
        """
        Polish a cluster by calculating its consensus sequence using RACON.
        
        Args:
            cluster_file: Path to cluster FASTA file
            cluster_id: Cluster identifier
            temp_dir: Temporary directory for intermediate files
            output_dir: Output directory for consensus files
            
        Returns:
            SeqRecord: Consensus sequence record, or None if failed
        """
        logger.info(f"Polishing cluster {cluster_id} with RACON consensus calculation")
        
        # Debug: Check cluster file content
        try:
            reads = list(SeqIO.parse(cluster_file, 'fasta'))
            logger.info(f"Cluster {cluster_id} contains {len(reads)} sequences")
            for i, read in enumerate(reads):
                logger.debug(f"  Read {i+1}: {read.id} ({len(read.seq)} bp)")
        except Exception as e:
            logger.error(f"Error reading cluster file {cluster_file}: {e}")
            raise RuntimeError(f"Error reading cluster file {cluster_file}: {e}")
        
        # Check if we have enough reads
        if len(reads) < 2:
            logger.warning(f"Cluster {cluster_id} has only {len(reads)} sequence(s), cannot calculate consensus")
            return None
        
        # Randomly sample reads if cluster is too large
        if len(reads) > self.racon_params['max_reads_for_polishing']:
            random.seed(42)  # For reproducible results
            sampled_reads = random.sample(reads, self.racon_params['max_reads_for_polishing'])
            logger.info(f"Cluster {cluster_id}: Sampling {len(sampled_reads)} reads from {len(reads)} total reads for RACON")
            reads_to_use = sampled_reads
        else:
            reads_to_use = reads
            logger.info(f"Cluster {cluster_id}: Using all {len(reads)} reads for RACON")
        
        # Create temporary files
        reads_file = temp_dir / f"reads_{cluster_id}.fasta"
        reference_file = temp_dir / f"reference_{cluster_id}.fasta"
        
        # Write reads to temporary file
        SeqIO.write(reads_to_use, reads_file, 'fasta')
        
        # Use the longest read as reference
        reference_read = max(reads_to_use, key=lambda x: len(x.seq))
        SeqIO.write([reference_read], reference_file, 'fasta')
        
        logger.info(f"Using {reference_read.id} ({len(reference_read.seq)} bp) as reference for cluster {cluster_id}")
        
        # Try RACON polishing
        try:
            # Step 1: Run minimap2 alignment
            paf_file = self.run_minimap2_alignment(reads_file, reference_file, temp_dir, cluster_id)
            if not paf_file:
                logger.warning(f"minimap2 alignment failed for cluster {cluster_id}, falling back to simple consensus")
                consensus = self.create_simple_consensus(cluster_file, cluster_id)
            else:
                # Step 2: Run RACON polishing
                consensus_file = self.run_racon_polishing(reads_file, paf_file, reference_file, temp_dir, cluster_id)
                if consensus_file and consensus_file.exists():
                    # Read the RACON consensus
                    consensus_records = list(SeqIO.parse(consensus_file, 'fasta'))
                    if consensus_records:
                        consensus = consensus_records[0]
                        consensus.id = f"cluster_{cluster_id:04d}_racon_consensus"
                        consensus.description = f"RACON consensus for cluster {cluster_id} ({len(consensus.seq)} bp)"
                        logger.info(f"RACON consensus created for cluster {cluster_id}: {len(consensus.seq)} bp")
                    else:
                        logger.warning(f"RACON produced empty consensus for cluster {cluster_id}, falling back to simple consensus")
                        consensus = self.create_simple_consensus(cluster_file, cluster_id)
                else:
                    logger.warning(f"RACON polishing failed for cluster {cluster_id}, falling back to simple consensus")
                    consensus = self.create_simple_consensus(cluster_file, cluster_id)
            
            if consensus:
                # Write consensus to cluster-specific file
                consensus_output_file = output_dir / f"cluster_{cluster_id:04d}_consensus.fasta"
                SeqIO.write([consensus], consensus_output_file, 'fasta')
                
                logger.info(f"RACON consensus written to: {consensus_output_file}")
                return consensus
            else:
                logger.error(f"Failed to create consensus for cluster {cluster_id} using all methods")
                raise RuntimeError(f"Failed to create consensus for cluster {cluster_id} using all methods")
                
        except Exception as e:
            logger.error(f"Error during RACON polishing for cluster {cluster_id}: {e}")
            # Fallback to simple consensus
            logger.info(f"Falling back to simple consensus for cluster {cluster_id}")
            consensus = self.create_simple_consensus(cluster_file, cluster_id)
            if consensus:
                consensus_output_file = output_dir / f"cluster_{cluster_id:04d}_consensus.fasta"
                SeqIO.write([consensus], consensus_output_file, 'fasta')
                logger.info(f"Simple consensus written to: {consensus_output_file}")
                return consensus
            else:
                raise RuntimeError(f"Failed to create consensus for cluster {cluster_id} using all methods")
