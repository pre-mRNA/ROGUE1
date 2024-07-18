import pandas as pd
import numpy as np
import argparse
from tqdm import tqdm
import time

def import_process_data(file_path):
    """Read and process read_id, gene_id, splice_count and intronic_alignment from ROGUE1 table."""
    col_types = {
        'read_id': str,
        'gene_id': str,
        'splice_count': int,
        'intronic_alignment': float
    }
    use_columns = list(col_types.keys())  
    
    return pd.read_csv(file_path, sep='\t', dtype=col_types, usecols=use_columns, na_values=['.'], low_memory=False)

def read_gtf(gtf_file):
    """Efficiently read GTF file to extract gene_ids for protein-coding genes with at least two exons."""
    gtf_columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    gtf_data = []
    with open(gtf_file, 'r') as file:
        for line in tqdm(file, desc="Filtering GTF file"):
            if line.startswith('#') or 'gene_biotype "protein_coding"' not in line or 'exon_number "2";' not in line:
                continue
            values = line.strip().split('\t')
            if values[2] == 'exon':
                gtf_data.append(values)

    gtf = pd.DataFrame(gtf_data, columns=gtf_columns)
    gtf['gene_id'] = gtf['attribute'].str.extract(r'gene_id "([^"]+)"')
    return gtf['gene_id'].dropna().unique()

def main(input_file, output_file, gtf_file, splice_threshold, intron_positive_threshold, intron_negative_threshold):
    
    print("Reading GTF file")
    start_time = time.time()
    protein_coding_intron_containing_genes = read_gtf(gtf_file)
    end_time = time.time()
    print(f"Done reading GTF file in {end_time - start_time:.2f} seconds")
    
    print("Reading ROGUE1 table")
    start_time = time.time()
    data = import_process_data(input_file)
    end_time = time.time()
    print(f"Done reading ROGUE1 table in {end_time - start_time:.2f} seconds")

    initial_read_count = len(data)

    data = data[data['gene_id'].isin(protein_coding_intron_containing_genes)]
    filtered_read_count = len(data)

    # data['splicing_classification'] = np.where(
    #     (data['splice_count'] >= splice_threshold) & 
    #     (data['intronic_alignment'] <= intron_negative_threshold), 'spliced',
    #     np.where(data['intronic_alignment'] > intron_positive_threshold, 'intron_retained', 'ambiguous')
    # )

    # classify as fully-spliced, intron-retained (partially spliced is IR), or ambigous 
    data['splicing_classification'] = np.where(
        (data['splice_count'] >= splice_threshold) & 
        (data['intronic_alignment'] <= intron_negative_threshold), 'spliced',
        np.where(
            data['intronic_alignment'] > intron_positive_threshold, 'intron_retained', 'ambiguous'
        )
    )

    print(f"Using thresholds - Splice: {splice_threshold}, Intron Positive: {intron_positive_threshold}, Intron Negative: {intron_negative_threshold}")
    print(f"Initial read count: {initial_read_count}")
    print(f"Filtered read count (protein-coding with introns): {filtered_read_count}")
    print(f"Spliced reads: {len(data[data['splicing_classification'] == 'spliced'])}")
    print(f"Intron retained reads: {len(data[data['splicing_classification'] == 'intron_retained'])}")
    print(f"Ambiguous reads: {len(data[data['splicing_classification'] == 'ambiguous'])}")
    print(f"Discarded reads: {initial_read_count - filtered_read_count}")

    output_data = data[['read_id', 'splicing_classification']]
    output_data.columns = ['read_id', 'classification']
    output_data.to_csv(output_file, sep='\t', header=True, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process ROGUE tables and classify reads.')
    parser.add_argument('input_file', type=str, help='Path to the input ROGUE table.')
    parser.add_argument('output_file', type=str, help='Path to the output file for saving results.')
    parser.add_argument('gtf_file', type=str, help='Path to the GTF file for annotation.')
    parser.add_argument('--splice_threshold', type=int, default=1, help='Threshold for splice count to classify as spliced.')
    parser.add_argument('--intron_positive_threshold', type=int, default=75, help='Threshold for intronic alignment to classify as intron retained.')
    parser.add_argument('--intron_negative_threshold', type=int, default=25, help='Threshold for intronic alignment to classify as spliced.')

    args = parser.parse_args()

    main(args.input_file, args.output_file, args.gtf_file, args.splice_threshold, args.intron_positive_threshold, args.intron_negative_threshold)
