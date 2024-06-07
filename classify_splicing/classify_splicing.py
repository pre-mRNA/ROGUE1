import pandas as pd
import numpy as np
import argparse
from tqdm import tqdm

def import_process_data(file_path):
    """Read and process ROGUE tables."""
    col_types = {
        'read_id': str,
        'gene_id': str,
        'read_length': int,
        'alignment_length': int,
        'splice_count': int,
        'gene_base_overlap': float,
        'exon_base_overlap': float,
        'intronic_alignment': float,
        'DOG_overlap': float,
        'unclassified_length': float,
        'unaligned_length': int,
        'read_end_chromosome': str,
        'read_end_position': int,
        'read_end_strand': str,
        'strand_sign': str,
        'downstream_distance': float,
        'upstream_distance': float
    }
    
    return pd.read_csv(file_path, sep='\t', dtype=col_types, na_values=['.'], low_memory=False)


def read_gtf(gtf_file):
    """Read GTF file into pandas DataFrame."""
    gtf_columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    gtf_data = []
    with open(gtf_file, 'r') as file:
        for line in tqdm(file, desc="Reading GTF file"):
            if not line.startswith('#'):
                values = line.strip().split('\t')
                gtf_data.append(values)
    gtf = pd.DataFrame(gtf_data, columns=gtf_columns)

    # Extract attributes into separate columns
    def extract_attributes(attribute_str):
        attributes = {}
        for attribute in attribute_str.split('; '):
            key_value = attribute.split(' ')
            if len(key_value) == 2:
                key, value = key_value
                attributes[key] = value.replace('"', '')
        return attributes

    attributes_df = gtf['attribute'].apply(extract_attributes).apply(pd.Series)
    gtf = pd.concat([gtf, attributes_df], axis=1)

    gtf['exon_number'] = gtf['exon_number'].astype(float)
    return gtf

def main(input_file, output_file, gtf_file, splice_threshold, intron_positive_threshold, intron_negative_threshold):
    
    # Load GTF file
    annotation = read_gtf(gtf_file)

    # Filter for protein-coding genes that contain introns
    protein_coding_intron_containing_genes = annotation[
        (annotation['feature'] == 'exon') & 
        (annotation['exon_number'] == 2) & 
        (annotation['gene_biotype'] == 'protein_coding')
    ]['gene_id'].unique()

    # Load dataset
    data = import_process_data(input_file)

    initial_read_count = len(data)

    # Calculate intronic alignment percentage
    data['percent_intronic_alignment'] = data['intronic_alignment'] / data['alignment_length']
    data['percent_intronic_alignment'].replace([np.inf, -np.inf], np.nan, inplace=True)  # Replace infinite values with NaN
    data.dropna(subset=['percent_intronic_alignment'], inplace=True)  # Drop rows with NaN values

    # Filter for protein-coding genes with introns
    data = data[data['gene_id'].isin(protein_coding_intron_containing_genes)]
    filtered_read_count = len(data)

    # Classify reads
    data['splicing_classification'] = np.where(
        (data['splice_count'] >= splice_threshold) & 
        (data['intronic_alignment'] <= intron_negative_threshold), 'spliced',
        np.where(data['intronic_alignment'] > intron_positive_threshold, 'intron_retained', 'ambiguous')
    )

    # Print summary
    print(f"Using thresholds - Splice: {splice_threshold}, Intron Positive: {intron_positive_threshold}, Intron Negative: {intron_negative_threshold}")
    print(f"Initial read count: {initial_read_count}")
    print(f"Filtered read count (protein-coding with introns): {filtered_read_count}")
    print(f"Spliced reads: {len(data[data['splicing_classification'] == 'spliced'])}")
    print(f"Intron retained reads: {len(data[data['splicing_classification'] == 'intron_retained'])}")
    print(f"Ambiguous reads: {len(data[data['splicing_classification'] == 'ambiguous'])}")
    print(f"Discarded reads: {initial_read_count - filtered_read_count}")

    # Save classifications to a single output file
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
