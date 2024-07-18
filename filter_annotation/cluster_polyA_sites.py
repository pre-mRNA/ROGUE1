import pandas as pd
import numpy as np
import argparse
import logging
from collections import defaultdict

def setup_logger():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    return logging.getLogger(__name__)

def find_clusters(positions, window_size, min_cluster_size):
    positions.sort()
    clusters = []
    current_cluster = [positions[0]]
    
    for pos in positions[1:]:
        if pos - current_cluster[0] <= window_size:
            current_cluster.append(pos)
        else:
            if len(current_cluster) >= min_cluster_size:
                clusters.append(current_cluster)
            current_cluster = [pos]
    
    if len(current_cluster) >= min_cluster_size:
        clusters.append(current_cluster)
    
    return clusters

def process_gene(gene_id, gene_data, window_size, min_cluster_size, logger, co_trans_cutoff, post_trans_cutoff, no_polyA_filtering):
    results = []
    cluster_numbers = defaultdict(int)
    
    if no_polyA_filtering:
        categories = [('all', slice(None))]
    else:
        categories = [
            ('co_transcriptional', gene_data['polya_length'] <= co_trans_cutoff),
            ('post_transcriptional', gene_data['polya_length'] >= post_trans_cutoff)
        ]
    
    for category, condition in categories:
        gene_data_filtered = gene_data[condition]
        
        if len(gene_data_filtered) < min_cluster_size:
            logger.info(f"Skipping {gene_id} for {category}: fewer than {min_cluster_size} reads")
            continue
        
        clusters = find_clusters(gene_data_filtered['read_end_position'].values, window_size, min_cluster_size)
        
        for cluster in clusters:
            cluster_df = gene_data_filtered[gene_data_filtered['read_end_position'].isin(cluster)]
            median_index = cluster_df['read_end_position'].argmin()
            median_row = cluster_df.iloc[median_index]
            min_distance = min(abs(median_row['polyA_downstream_distance']), 
                               abs(median_row['polyA_upstream_distance']))
            
            cluster_numbers[category] += 1
            cluster_name = f"{gene_id}_{category}_{cluster_numbers[category]}"
            
            results.append({
                'gene_id': gene_id,
                'category': category,
                'median_position': median_row['read_end_position'],
                'nearest_distance': min_distance,
                'median_chromosome': median_row['read_end_chromosome'],
                'median_strand': median_row['read_end_strand'],
                'read_count': len(cluster),
                'cluster_name': cluster_name
            })
    
    return pd.DataFrame(results) if results else None

def clusters_to_bed(clusters_df, output_file):
    bed_df = pd.DataFrame({
        'chrom': clusters_df['median_chromosome'],
        'chromStart': clusters_df['median_position'] - 1,
        'chromEnd': clusters_df['median_position'],
        'name': clusters_df['cluster_name'],
        'score': clusters_df['read_count'],
        'strand': clusters_df['median_strand'],
        'nearest_distance': clusters_df['nearest_distance']
    })

    bed_df['nearest_distance'] = bed_df['nearest_distance'].fillna(-1)
    bed_df['nearest_distance'] = bed_df['nearest_distance'].replace([np.inf, -np.inf], -1)
    bed_df['nearest_distance'] = bed_df['nearest_distance'].round().astype(int)

    bed_df = bed_df.astype({
        'chrom': str,
        'chromStart': int,
        'chromEnd': int,
        'name': str,
        'score': int,
        'strand': str,
        'nearest_distance': int
    })

    bed_df.to_csv(output_file, sep='\t', header=False, index=False)

def main(args):
    logger = setup_logger()
    
    logger.info(f"Reading input file: {args.input}")
    try:
        df = pd.read_csv(args.input, sep='\t')
    except Exception as e:
        logger.error(f"Error reading input file: {e}")
        return

    logger.info("Filtering for protein-coding genes")
    df_protein_coding = df[df['biotype'] == 'protein_coding']

    if args.no_polyA_filtering:
        logger.info("Poly(A) filtering disabled. Using all reads for clustering.")
        co_trans_cutoff = post_trans_cutoff = None
    elif args.use_polyA_cutoffs:
        co_trans_cutoff, post_trans_cutoff = map(int, args.use_polyA_cutoffs.split(','))
        logger.info(f"Using custom poly(A) cutoffs: co-transcriptional <= {co_trans_cutoff}, post-transcriptional >= {post_trans_cutoff}")
    else:
        co_trans_cutoff, post_trans_cutoff = 50, 100
        logger.info(f"Using default poly(A) cutoffs: co-transcriptional <= {co_trans_cutoff}, post-transcriptional >= {post_trans_cutoff}")

    logger.info("Processing genes and finding clusters")
    grouped = df_protein_coding.groupby('gene_id')
    all_clusters = []

    for gene_id, group in grouped:
        try:
            gene_result = process_gene(gene_id, group, args.window_size, args.min_cluster_size, logger, 
                                       co_trans_cutoff, post_trans_cutoff, args.no_polyA_filtering)
            if gene_result is not None:
                all_clusters.append(gene_result)
        except Exception as e:
            logger.error(f"Error processing gene {gene_id}: {e}")

    if not all_clusters:
        logger.warning("No clusters found")
        return

    all_clusters_df = pd.concat(all_clusters, ignore_index=True)

    logger.info("Adding cluster information to original data")
    df['cluster'] = ''
    for _, row in all_clusters_df.iterrows():
        mask = (df['gene_id'] == row['gene_id']) & \
               (df['read_end_position'] == row['median_position']) & \
               (df['read_end_chromosome'] == row['median_chromosome']) & \
               (df['read_end_strand'] == row['median_strand'])
        df.loc[mask, 'cluster'] = row['cluster_name']

    logger.info(f"Writing output file: {args.output}")
    df.to_csv(args.output, sep='\t', index=False)

    if args.bed:
        logger.info(f"Writing BED file: {args.bed}")
        clusters_to_bed(all_clusters_df, args.bed)

    logger.info("Processing completed successfully")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Perform clustering on gene data from nanopore direct RNA sequencing.",
        epilog="""
Examples:
  # Using default settings
  python clustering_script.py input.tsv output.tsv

  # Disabling poly(A) filtering
  python clustering_script.py input.tsv output.tsv --no_polyA_filtering

  # Using custom poly(A) cutoffs
  python clustering_script.py input.tsv output.tsv --use_polyA_cutoffs 40,120

  # Generating a BED file and customizing clustering parameters
  python clustering_script.py input.tsv output.tsv --bed output.bed --window_size 50 --min_cluster_size 15
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("input", help="Path to input TSV file containing nanopore sequencing data")
    parser.add_argument("output", help="Path to output TSV file with added cluster information")
    parser.add_argument("--bed", help="Path to output BED file summarizing cluster information (optional)")
    parser.add_argument("--window_size", type=int, default=30, 
                        help="Window size for clustering in base pairs (default: 30)")
    parser.add_argument("--min_cluster_size", type=int, default=10, 
                        help="Minimum number of reads to form a cluster (default: 10)")
    
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--no_polyA_filtering", action="store_true", 
                       help="Disable poly(A) filtering and use all reads for clustering")
    group.add_argument("--use_polyA_cutoffs", metavar="INT1,INT2",
                       help="Manually select poly(A) cutoffs. Format: co_transcriptional_max,post_transcriptional_min")
    
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')

    args = parser.parse_args()

    main(args)