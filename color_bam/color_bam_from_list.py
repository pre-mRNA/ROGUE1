import argparse
import pysam
import pandas as pd

# mapping of classifications to colors
CLASSIFICATION_COLOR_MAP = {
    'intron': '255,0,112',
    'spliced': '51,153,255',
    'ambiguous': '255,102,0',
    'intron_retained': '255,0,112'  # added mapping for intron_retained
}

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

def load_classifications(file_path):
    # read classification file into a dataframe
    df = pd.read_csv(file_path, sep='\t')
    df = df[['read_id', 'classification']]

    # map classifications to colors and validate
    df['color'] = df['classification'].map(CLASSIFICATION_COLOR_MAP)
    
    # check for any NaN values in the 'color' column
    if df['color'].isna().any():
        missing_classes = df[df['color'].isna()]['classification'].unique()
        raise ValueError(f"Unrecognized classifications found: {', '.join(missing_classes)}")

    # validate colors
    df['color'] = df['color'].apply(validate_color)
    
    return df.set_index('read_id')['color'].to_dict()

def main():
    parser = argparse.ArgumentParser(description="Set specific RGB colors for read IDs in a BAM file based on classifications.")
    parser.add_argument("-ibam", required=True, help="Input BAM file")
    parser.add_argument("-obam", required=True, help="Output BAM file")
    parser.add_argument("-class_file", required=True, help="File containing read IDs and their classifications")

    args = parser.parse_args()

    # load classifications from the provided file
    classifications = load_classifications(args.class_file)

    modified_count = 0  # count modified reads
    missing_reads = set(classifications.keys())  # initialize missing reads with all read_ids from classification file

    with pysam.AlignmentFile(args.ibam, 'rb', threads=8) as infile, \
         pysam.AlignmentFile(args.obam, 'wb', template=infile, threads=8) as outfile:
        for read in infile:
            if read.query_name in classifications:
                read.set_tag("YC", classifications[read.query_name])
                modified_count += 1
                missing_reads.discard(read.query_name)
            outfile.write(read)

    if missing_reads:
        print(f"Warning: the following reads were not found in the BAM file: {', '.join(missing_reads)}")

    print(f"sorting and indexing output bam file")

    # sort and index the output BAM file
    pysam.sort("-o", f"{args.obam}.sorted.bam", args.obam, threads=8)
    pysam.index(f"{args.obam}.sorted.bam", threads=8)

    print(f"total reads modified: {modified_count}")

if __name__ == "__main__":
    main()
