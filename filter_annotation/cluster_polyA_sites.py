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

def process_gene(gene_id, gene_data, window_size, min_cluster_size, logger, cluster_column, co_trans_cutoff, post_trans_cutoff, no_polyA_filtering):
    results = []
    cluster_numbers = defaultdict(int)
    
    if no_polyA_filtering:
        categories = [('all', slice(None))]
    else:
        categories = [
            ('co_transcriptional', gene_data[cluster_column] <= co_trans_cutoff),
            ('post_transcriptional', gene_data[cluster_column] >= post_trans_cutoff)
        ]
    
    for category, condition in categories:
        gene_data_filtered = gene_data[condition]
        
        if len(gene_data_filtered) < min_cluster_size:
            logger.info(f"skipping {gene_id} for {category}: fewer than {min_cluster_size} reads")
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
    
    logger.info(f"reading input file: {args.input}")
    try:
        df = pd.read_csv(args.input, sep='\t')
    except Exception as e:
        logger.error(f"error reading input file: {e}")
        return

    logger.info("filtering for protein-coding genes")
    df_protein_coding = df[df['biotype'] == 'protein_coding']

    # filter based on minimum alignment length if provided
    if args.min_alignment_length is not None:
        if 'alignment_length' not in df_protein_coding.columns:
            logger.error("alignment_length column is missing from the input data")
            return
        df_protein_coding = df_protein_coding[df_protein_coding['alignment_length'] >= args.min_alignment_length]
        logger.info(f"filtered reads with alignment_length >= {args.min_alignment_length}")

    if args.no_polyA_filtering:
        logger.info("poly(A) filtering disabled. using all reads for clustering.")
        co_trans_cutoff = post_trans_cutoff = None
    else:
        if not args.cluster_cutoffs:
            logger.error("cluster_cutoffs are required unless --no_polyA_filtering is specified.")
            return
        
        try:
            co_trans_cutoff, post_trans_cutoff = map(float, args.cluster_cutoffs.split(','))
        except ValueError:
            logger.error("cluster_cutoffs must be two comma-separated numbers.")
            return
        
        logger.info(f"using custom cluster cutoffs: co_trans_cutoff <= {co_trans_cutoff}, post_trans_cutoff >= {post_trans_cutoff}")
    
    logger.info("processing genes and finding clusters")
    grouped = df_protein_coding.groupby('gene_id')
    all_clusters = []

    for gene_id, group in grouped:
        try:
            gene_result = process_gene(
                gene_id, 
                group, 
                args.window_size, 
                args.min_cluster_size, 
                logger, 
                args.cluster_by,
                co_trans_cutoff, 
                post_trans_cutoff, 
                args.no_polyA_filtering
            )
            if gene_result is not None:
                all_clusters.append(gene_result)
        except Exception as e:
            logger.error(f"error processing gene {gene_id}: {e}")

    if not all_clusters:
        logger.warning("no clusters found")
        return

    all_clusters_df = pd.concat(all_clusters, ignore_index=True)

    # logger.info("adding cluster information to original data")
    # df['cluster'] = ''
    # for _, row in all_clusters_df.iterrows():
    #     mask = (df['gene_id'] == row['gene_id']) & \
    #            (df['read_end_position'] == row['median_position']) & \
    #            (df['read_end_chromosome'] == row['median_chromosome']) & \
    #            (df['read_end_strand'] == row['median_strand'])
    #     df.loc[mask, 'cluster'] = row['cluster_name']

    # logger.info(f"writing output file: {args.output}")
    # df.to_csv(args.output, sep='\t', index=False)

    if args.bed:
        logger.info(f"writing BED file: {args.bed}")
        clusters_to_bed(all_clusters_df, args.bed)

    logger.info("processing completed successfully")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="perform clustering on gene data from nanopore direct RNA sequencing.",
        epilog="""
Examples:
  # using default settings
  python clustering_script.py input.tsv output.tsv --cluster_by polya_length --cluster_cutoffs 50,100

  # disabling poly(A) filtering
  python clustering_script.py input.tsv output.tsv --no_polyA_filtering --cluster_by polya_z_score --cluster_cutoffs 1.5,2.5

  # generating a BED file and customizing clustering parameters
  python clustering_script.py input.tsv output.tsv --bed output.bed --window_size 50 --min_cluster_size 15 --cluster_by polya_length --cluster_cutoffs 40,120

  # specifying minimum alignment length
  python clustering_script.py input.tsv output.tsv --min_alignment_length 500 --cluster_by polya_length --cluster_cutoffs 50,100
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("input", help="path to input TSV file containing nanopore sequencing data")
    parser.add_argument("output", help="path to output TSV file with added cluster information")
    parser.add_argument("--bed", help="path to output BED file summarizing cluster information (optional)")
    parser.add_argument("--window_size", type=int, default=30, 
                        help="window size for clustering in base pairs (default: 30)")
    parser.add_argument("--min_cluster_size", type=int, default=10, 
                        help="minimum number of reads to form a cluster (default: 10)")
    
    parser.add_argument("--cluster_by", required=True, choices=['polya_length', 'polya_z_score'],
                        help="column to cluster by: 'polya_length' or 'polya_z_score'")
    parser.add_argument("--cluster_cutoffs", help="comma-separated co_trans_cutoff and post_trans_cutoff (e.g., 50,100)", required=False)
    
    parser.add_argument("--min_alignment_length", type=int,
                        help="minimum alignment length to include reads in clustering")
    
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--no_polyA_filtering", action="store_true", 
                       help="disable poly(A) filtering and use all reads for clustering")
    group.add_argument("--use_polyA_cutoffs", metavar="INT1,INT2",
                       help="manually select poly(A) cutoffs. format: co_transcriptional_max,post_transcriptional_min")
    
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')

    args = parser.parse_args()

    main(args)
