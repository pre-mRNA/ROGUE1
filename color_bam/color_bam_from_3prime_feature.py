import argparse
import pysam
import pandas as pd
import os
import sys
from bam_header_utils import update_pg_header

# map 3' features to colors
FEATURE_COLOR_MAP = {
    'exon': '51,153,255',  # blue
    'intron': '255,0,112',  # pink
    'DOG': '0,255,0',  # green
    'other': '128,128,128'  # grey
}

# fetch R1 version for header 
def get_version():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(script_dir)
    version_file = os.path.join(parent_dir, 'version.txt')
    try:
        with open(version_file, 'r') as f:
            return f.read().strip()
    except FileNotFoundError:
        return "unknown"

# classify feature type from 'three_prime_feature' column
def get_feature_type(feature):
    if pd.isna(feature) or feature == '':
        return 'other'
    if feature.startswith('exon_') or feature.startswith('intron_'):
        return feature.split('_')[0]
    if feature.endswith('DOG'):
        return 'DOG'
    return 'other'

def load_features(file_path):
    
    # read R1 table 
    df = pd.read_csv(file_path, sep='\t')
    df = df[['read_id', 'three_prime_feature']]

    # map features to colors
    df['feature_type'] = df['three_prime_feature'].apply(get_feature_type)
    df['color'] = df['feature_type'].map(FEATURE_COLOR_MAP)

    return df.set_index('read_id')[['color', 'three_prime_feature']].to_dict('index')

def main():
    parser = argparse.ArgumentParser(description="Set RGB colors and EF tag for read IDs in a BAM file based on 3' features.")
    parser.add_argument("-ibam", required=True, help="Input BAM file")
    parser.add_argument("-obam", required=True, help="Output BAM file")
    parser.add_argument("-table", required=True, help="Output table containing read IDs and their 3' features")

    args = parser.parse_args()

    # load 3' features from R1 table 
    features = load_features(args.table)

    modified_count = 0
    missing_reads = set(features.keys())

    # prepare bam header 
    command_line = ' '.join(['python'] + sys.argv)
    description = 'Colored reads based on 3\' features and added EF tag'

    with pysam.AlignmentFile(args.ibam, 'rb', threads=48) as infile:
        header = infile.header.to_dict()
        header = update_pg_header(header, command_line, description)

        with pysam.AlignmentFile(args.obam, 'wb', header=header, threads=48) as outfile:
            for read in infile:
                if read.query_name in features:
                    color = features[read.query_name]['color']
                    end_feature = features[read.query_name]['three_prime_feature']
                    read.set_tag("YC", color) # update YC tag for color 
                    read.set_tag("EF", end_feature if pd.notna(end_feature) else "other", value_type='Z') # update EF tag for classification 
                    modified_count += 1
                    missing_reads.discard(read.query_name)
                outfile.write(read)

    if missing_reads:
        print(f"Warning: {len(missing_reads)} reads from the table were not found in the BAM file.")

    print(f"Indexing output BAM file")
    pysam.index(f"{args.obam}", threads=48)

    print(f"Total reads modified: {modified_count}")
    print("Added tags:")
    print("  YC: RGB color based on 3' feature")
    print("  EF: Exact 3' end feature")

if __name__ == "__main__":
    main()