import os 
import sys 
import argparse
import pysam
import pandas as pd
import numpy as np
from bam_header_utils import update_pg_header

def get_version():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(script_dir)
    version_file = os.path.join(parent_dir, 'version.txt')
    try:
        with open(version_file, 'r') as f:
            return f.read().strip()
    except FileNotFoundError:
        return "unknown"

def calculate_color_sliding(z_score):
    low_threshold = 1
    high_threshold = 2.33
    gray_color = "128,128,128"  # RGB for grey

    if pd.isna(z_score) or z_score == -1:
        return gray_color
    if z_score <= low_threshold:
        return "255,0,0"  # red
    elif z_score >= high_threshold:
        return "0,0,255"  # blue
    else:
        # scale z_score to a 0-1 range
        scaled_value = (z_score - low_threshold) / (high_threshold - low_threshold)
        
        # interpolate between red (255,0,0) and blue (0,0,255)
        red = int(255 * (1 - scaled_value))
        blue = int(255 * scaled_value)
        return f"{red},0,{blue}"

def calculate_color_quantized(z_score):
    gray_color = "128,128,128"  # RGB for grey

    if pd.isna(z_score) or z_score == -1:
        return gray_color, 2
    if z_score < 1:
        return "255,0,0", 0  # red
    elif z_score > 2.33:    
        return "0,0,255", 1  # blue
    else:
        return "128,0,128", 2  # purple for values between thresholds

def main():
    parser = argparse.ArgumentParser(description="Color reads in a BAM file based on polyA z-scores.")
    parser.add_argument("-ibam", required=True, help="Input BAM file")
    parser.add_argument("-obam", required=True, help="Output BAM file")
    parser.add_argument("-table", required=True, help="Table containing read_id and polya_zscore columns")
    parser.add_argument("-q", action="store_true", help="Use quantized colors instead of sliding scale")

    args = parser.parse_args()

    # Read only the necessary columns from the table
    df = pd.read_csv(args.table, sep='\t', usecols=['read_id', 'polya_zscore'])
    z_scores = dict(zip(df['read_id'], df['polya_zscore']))

    modified_count = 0  # count of modified reads
    
    # prepare bam header entry 
    command_line = ' '.join(['python'] + sys.argv)
    description = 'Colored reads based on polyA z-scores and added QP tag for quantized polyA'
    
    with pysam.AlignmentFile(args.ibam, 'rb', threads=8) as infile:
        header = infile.header.to_dict()
        header = update_pg_header(header, command_line, description)
        
        # write to outfile 
        with pysam.AlignmentFile(args.obam, 'wb', header=header, threads=8) as outfile:

            for read in infile:
                z_score = z_scores.get(read.query_name, np.nan)
                
                if args.q:
                    color, qp_value = calculate_color_quantized(z_score)
                    read.set_tag("YC", color)
                    read.set_tag("QP", qp_value, value_type='i')
                else:
                    color = calculate_color_sliding(z_score)
                    read.set_tag("YC", color)
                
                if not pd.isna(z_score) and z_score != -1:
                    read.set_tag("ZS", f"{z_score:.2f}")
                
                outfile.write(read)
                modified_count += 1

    print(f"Sorting and indexing output BAM file")
    pysam.index(f"{args.obam}", threads=8)

    print(f"Total reads modified: {modified_count}")
    if args.q:
        print("Quantized polyA color: Red < 1, 1-2.33 Purple, Blue > 2.33, Gray for missing/invalid")
        print("QP tag: 0 for z-score < 1, 1 for z-score > 2.33, 2 for 1 <= z-score <= 2.33 or missing/invalid")
    else:
        print("Scaled polyA color: Red <= 1, 1-2.33 intermediate, Blue >= 2.33, Gray for missing/invalid")

if __name__ == "__main__":
    main()