import argparse
import pysam
import pandas as pd
import os
import sys
from bam_header_utils import update_pg_header

# mapping classifications to colors
CLASSIFICATION_COLOR_MAP = {
    'intron': '255,0,112',
    'spliced': '51,153,255',
    'ambiguous': '255,102,0',
    'intron_retained': '255,0,112',  # added mapping for intron_retained
    'fully-unspliced': '255,0,112' # fully unspliced is the same as intron
}

def get_version():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(script_dir)
    version_file = os.path.join(parent_dir, 'version.txt')
    try:
        with open(version_file, 'r') as f:
            return f.read().strip()
    except FileNotFoundError:
        return "unknown"

def validate_color(color):
    # validate that the color is a comma-separated string of three integers between 0-255
    if isinstance(color, float) and pd.isna(color):
        raise ValueError("color value is missing or not a string.")
    parts = color.split(',')
    if len(parts) != 3:
        raise ValueError("color must be a comma-separated string of three integers between 0-255.")
    for part in parts:
        if not (0 <= int(part) <= 255):
            raise ValueError("each RGB value must be between 0 and 255.")
    return color

def load_splicing_classifications(file_path):

    # read classification file into a dataframe
    df = pd.read_csv(file_path, sep='\t')
    df = df[['read_id', 'splicing_classification']]

    # map splicing_classifications to colors and validate
    df['color'] = df['splicing_classification'].map(CLASSIFICATION_COLOR_MAP)
    
    # check for any NaN values in the 'color' column
    if df['color'].isna().any():
        missing_classes = df[df['color'].isna()]['splicing_classification'].unique()
        raise ValueError(f"Unrecognized splicing_classifications found: {', '.join(missing_classes)}")

    # validate colors
    df['color'] = df['color'].apply(validate_color)
    
    return df.set_index('read_id')[['color', 'splicing_classification']].to_dict('index')

def main():
    parser = argparse.ArgumentParser(description="Set specific RGB colors for read IDs in a BAM file based on splicing_classifications.")
    parser.add_argument("-ibam", required=True, help="Input BAM file")
    parser.add_argument("-obam", required=True, help="Output BAM file")
    parser.add_argument("-class_file", required=True, help="File containing read IDs and their splicing_classifications")

    args = parser.parse_args()

    # load splicing_classifications from the provided file
    splicing_classifications = load_splicing_classifications(args.class_file)

    modified_count = 0  # count modified reads
    missing_reads = set(splicing_classifications.keys())  # initialize missing reads with all read_ids from splicing_classification file

    # prepare bam header entry 
    command_line = ' '.join(['python'] + sys.argv)
    description = 'Colored reads based on splicing splicing_classifications and added SC tag'
    
    with pysam.AlignmentFile(args.ibam, 'rb', threads=48) as infile:
        header = infile.header.to_dict()
        header = update_pg_header(header, command_line, description)

        with pysam.AlignmentFile(args.obam, 'wb', header=header, threads=48) as outfile:
            for read in infile:
                if read.query_name in splicing_classifications:
                    color = splicing_classifications[read.query_name]['color']
                    splicing_classification = splicing_classifications[read.query_name]['splicing_classification']
                    read.set_tag("YC", color)
                    read.set_tag("SC", splicing_classification, value_type='Z')
                    modified_count += 1
                    missing_reads.discard(read.query_name)
                outfile.write(read)

    if missing_reads:
        print(f"Warning: the following reads were not found in the BAM file: {', '.join(missing_reads)}")

    print(f"Indexing output bam file")

    # index the output bam
    pysam.index(f"{args.obam}", threads=48)

    print(f"Total reads modified: {modified_count}")
    print("Added tags:")
    print("  YC: RGB color based on splicing classification")
    print("  SC: Splicing classification")

if __name__ == "__main__":
    main()