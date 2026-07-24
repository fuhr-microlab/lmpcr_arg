#!/usr/bin/env python3
"""
VSEARCH-based clustering module.

Fast clustering with VSEARCH (typically much faster than BLAST-based clustering).
Cluster representatives are suitable for downstream RACON polishing (via minimap2).
"""

import subprocess
import logging
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from lib.default_params import VSEARCH_PARAMS

logger = logging.getLogger(__name__)


class VsearchClusterer:
    """VSEARCH-based clustering implementation."""
    
    def __init__(self, vsearch_params=None):
        """
        Initialize VSEARCH clusterer.
        
        Args:
            vsearch_params (dict): VSEARCH parameters (uses defaults if None)
        """
        self.vsearch_params = vsearch_params or VSEARCH_PARAMS.copy()
    
    def check_vsearch_availability(self):
        """Check if VSEARCH is available and get version."""
        try:
            result = subprocess.run(['vsearch', '--version'], capture_output=True, text=True, check=True)
            logger.info(f"VSEARCH version: {result.stdout.strip()}")
            return True
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            logger.error(f"VSEARCH is required for fast clustering but not found: {e}")
            logger.error("Please install VSEARCH: conda install -c bioconda vsearch")
            raise RuntimeError("VSEARCH is required for fast clustering. Please install it first.")
    
    def run_vsearch_clustering(self, reads, temp_dir):
        """Run VSEARCH clustering on the reads."""
        logger.info("Using VSEARCH for fast clustering")
        
        # Create temporary input file
        input_file = temp_dir / "reads_for_vsearch.fasta"
        SeqIO.write(reads, input_file, 'fasta')
        
        # VSEARCH clustering output files
        output_prefix = temp_dir / "vsearch_clusters"
        uc_file = temp_dir / "clusters.uc"
        
        # Run VSEARCH clustering with configurable parameters for better quality
        cmd = [
            'vsearch', f'--{self.vsearch_params["clustering_algorithm"]}', str(input_file),
            '--id', str(self.vsearch_params['identity_threshold']),
            '--strand', self.vsearch_params['strand'],
            '--uc', str(uc_file),
            '--centroids', str(output_prefix) + '.fasta'
        ]
        
        # Add coverage thresholds if specified
        if self.vsearch_params['query_coverage']:
            cmd.extend(['--query_cov', str(self.vsearch_params['query_coverage'])])
        if self.vsearch_params['target_coverage']:
            cmd.extend(['--target_cov', str(self.vsearch_params['target_coverage'])])
        
        # Add size annotations if enabled
        if self.vsearch_params['sizein']:
            cmd.append('--sizein')
        if self.vsearch_params['sizeout']:
            cmd.append('--sizeout')
        
        logger.info(f"Running VSEARCH clustering: {' '.join(cmd)}")
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info("VSEARCH clustering completed")
        
        # Parse VSEARCH output
        clusters = self.parse_vsearch_output(uc_file, reads)
        
        return clusters
    
    def parse_vsearch_output(self, uc_file, reads):
        """Parse VSEARCH UC output format and return clusters."""
        clusters = []
        cluster_id_to_index = {}  # Map VSEARCH cluster ID to cluster index
        
        # Read UC file
        with open(uc_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                fields = line.split('\t')
                if len(fields) < 10:
                    continue
                
                record_type = fields[0]
                cluster_id = fields[1]  # VSEARCH cluster identifier
                read_id = fields[8]     # Read identifier
                
                if record_type == 'S':  # Seed sequence (representative)
                    # Create new cluster
                    read_obj = next((r for r in reads if r.id == read_id), None)
                    if read_obj:
                        clusters.append([read_obj])
                        cluster_id_to_index[cluster_id] = len(clusters) - 1  # Map cluster ID to cluster index
                
                elif record_type == 'H':  # Hit sequence (member)
                    # Add to existing cluster
                    read_obj = next((r for r in reads if r.id == read_id), None)
                    if read_obj and cluster_id in cluster_id_to_index:  # Use cluster ID to find cluster
                        cluster_idx = cluster_id_to_index[cluster_id]
                        clusters[cluster_idx].append(read_obj)
        
        # Debug: Show cluster size distribution
        if clusters:
            cluster_sizes = [len(cluster) for cluster in clusters]
            logger.info(f"Original cluster size distribution: min={min(cluster_sizes)}, max={max(cluster_sizes)}, mean={sum(cluster_sizes)/len(cluster_sizes):.1f}")
        
        logger.info(f"VSEARCH created {len(clusters)} clusters")
        
        return clusters
    
    def filter_cluster_quality(self, cluster):
        """Filter VSEARCH cluster to remove sequences that are too divergent for consensus polishing."""
        if len(cluster) <= 2:
            return cluster  # Too small to filter
        
        # Use the longest sequence as reference
        reference = max(cluster, key=lambda x: len(x.seq))
        reference_length = len(reference.seq)
        
        # Calculate pairwise similarities and filter out divergent sequences
        filtered_cluster = [reference]  # Always keep reference
        
        for read in cluster:
            if read.id == reference.id:
                continue
            
            # Simple length-based similarity check
            length_diff = abs(len(read.seq) - reference_length)
            length_similarity = 1.0 - (length_diff / max(len(read.seq), reference_length))
            
            # Keep sequences that are reasonably similar in length (within 20%)
            if length_similarity >= 0.8:
                filtered_cluster.append(read)
            else:
                logger.debug(f"Filtered out divergent sequence {read.id}: length similarity {length_similarity:.2f}")
        
        logger.debug(f"Cluster quality filtering: {len(cluster)} -> {len(filtered_cluster)} sequences")
        return filtered_cluster
    
    def cluster_reads(self, reads, temp_dir, output_dir, min_cluster_size, use_quality_filtering=True, output_longest_reads=False):
        """
        Perform VSEARCH-based clustering.
        
        Args:
            reads: List of SeqRecord objects
            temp_dir: Temporary directory for intermediate files
            output_dir: Output directory for cluster files
            min_cluster_size: 15
            use_quality_filtering: Whether to filter clusters for quality
            output_longest_reads: Whether to collect longest reads from clusters
            
        Returns:
            tuple: (cluster_files_dict, cluster_list, longest_reads_list)
        """
        # Check VSEARCH availability
        self.check_vsearch_availability()
        
        # Run VSEARCH clustering
        clusters = self.run_vsearch_clustering(reads, temp_dir)
        
        # Filter clusters by minimum size and quality
        filtered_clusters = []
        for cluster in clusters:
            if len(cluster) >= min_cluster_size:
                # For VSEARCH clusters, filter out sequences that are too divergent
                if use_quality_filtering:
                    filtered_cluster = self.filter_cluster_quality(cluster)
                    if len(filtered_cluster) >= min_cluster_size:
                        filtered_clusters.append(filtered_cluster)
                else:
                    filtered_clusters.append(cluster)
        
        # Debug: Show filtered cluster size distribution
        if filtered_clusters:
            filtered_sizes = [len(cluster) for cluster in filtered_clusters]
            logger.info(f"Filtered cluster size distribution: min={min(filtered_sizes)}, max={max(filtered_sizes)}, mean={sum(filtered_sizes)/len(filtered_sizes):.1f}")
            logger.info(f"Clusters with size >= {min_cluster_size}: {len(filtered_clusters)}/{len(clusters)}")
            
            # Show some examples of small clusters
            small_clusters = [cluster for cluster in clusters if len(cluster) < min_cluster_size]
            if small_clusters:
                logger.info(f"Examples of small clusters (size < {min_cluster_size}):")
                for i, cluster in enumerate(small_clusters[:5]):  # Show first 5 small clusters
                    logger.info(f"  Small cluster {i+1}: {len(cluster)} reads")
        
        logger.info(f"VSEARCH created {len(clusters)} clusters, {len(filtered_clusters)} meet minimum size and quality requirements")
        
        # Process clusters and create cluster files
        cluster_files = self._process_clusters(filtered_clusters, output_dir, output_longest_reads)
        
        return cluster_files
    
    def _process_clusters(self, clusters, output_dir, output_longest_reads=False):
        """Process VSEARCH clusters and create cluster files."""
        cluster_files = {}
        cluster_id = 1
        cluster_list = []
        longest_reads_list = []
        growing_cluster_file = output_dir / "cluster_list_growing.fasta"
        
        # Initialize growing cluster list file
        if growing_cluster_file.exists():
            growing_cluster_file.unlink()
        
        for cluster in clusters:
            if len(cluster) >= 1:  # At least one read
                # Use first read as representative (they're already sorted by length)
                representative = cluster[0]
                cluster_members = cluster[1:]
                
                # Create cluster file
                cluster_filepath = self._create_cluster(representative, cluster_members, cluster_id, output_dir, cluster_list, growing_cluster_file)
                cluster_files[f"cluster_{cluster_id:04d}"] = cluster_filepath
                
                # Collect longest read if requested
                if output_longest_reads:
                    longest_read = max(cluster, key=lambda x: len(x.seq))
                    longest_read_copy = SeqRecord(
                        seq=longest_read.seq,
                        id=f"cluster_{cluster_id:04d}_longest",
                        description=f"Longest read from cluster {cluster_id} (original: {longest_read.id}, {len(longest_read.seq)} bp)"
                    )
                    longest_reads_list.append(longest_read_copy)
                
                logger.info(f"VSEARCH cluster {cluster_id} completed: {len(cluster)} reads")
                cluster_id += 1
        
        return cluster_files, cluster_list, longest_reads_list
    
    def _create_cluster(self, representative_read, cluster_reads, cluster_id, output_dir, cluster_list, growing_cluster_file):
        """Create a cluster file with the representative and cluster reads."""
        cluster_filename = f"cluster_{cluster_id:04d}.fasta"
        cluster_filepath = output_dir / cluster_filename
        
        # Write representative read first, then cluster reads
        all_reads = [representative_read] + cluster_reads
        SeqIO.write(all_reads, cluster_filepath, 'fasta')
        
        # Add representative to cluster list for tracking
        representative_copy = SeqRecord(
            seq=representative_read.seq,
            id=f"cluster_{cluster_id:04d}_representative",
            description=f"Representative for cluster {cluster_id} (original: {representative_read.id})"
        )
        cluster_list.append(representative_copy)
        
        # Write to growing cluster file immediately
        with open(growing_cluster_file, 'a') as f:
            SeqIO.write(representative_copy, f, 'fasta')
        
        logger.info(f"Added cluster {cluster_id} representative to growing list: {growing_cluster_file}")
        logger.info(f"Created cluster {cluster_id}: {len(all_reads)} reads -> {cluster_filepath}")
        
        return cluster_filepath
