import os 
import sys 
import argparse
import pysam
import pandas as pd
from bam_header_utils import update_pg_header

# fetch ROGUE1 version for BAM header
# TODO: move this to a python module 
def get_version():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(script_dir)
    version_file = os.path.join(parent_dir, 'version.txt')
    try:
        with open(version_file, 'r') as f:
            return f.read().strip()
    except FileNotFoundError:
        return "unknown"

# quantize splice distance color for different thresholds
# returns YC tag and QS tag 
def calculate_color_quantized(splice_distance, low_threshold, high_threshold):
    if splice_distance == -1:
        return "128,128,128", -1  # grey for unknown
    elif splice_distance <= low_threshold:
        return "255,0,0", 0  # red
    elif splice_distance <= high_threshold:
        return "255,255,0", 1  # yellow
    else:
        return "0,0,255", 2  # blue

def load_splice_distances(tsv_file):
    df = pd.read_csv(tsv_file, sep='\t', usecols=['read_id', 'splice_donors_downstream_distance', 'splice_donors_upstream_distance', 'biotype'])
    df['min_splice_distance'] = df.apply(lambda row: min(abs(row['splice_donors_downstream_distance']), 
                                                         abs(row['splice_donors_upstream_distance'])), axis=1)
    return df.set_index('read_id')[['min_splice_distance', 'biotype']].to_dict('index')

def main():
    parser = argparse.ArgumentParser(description="Color reads in a BAM file based on minimum distance to splice sites.")
    parser.add_argument("-ibam", required=True, help="Input BAM file")
    parser.add_argument("-obam", required=True, help="Output BAM file")
    parser.add_argument("-q", type=str, required=True, help="Comma-separated low and high thresholds for quantized splice distance color")
    parser.add_argument("--class-file", required=True, help="TSV file containing read_id and splice distance information")

    args = parser.parse_args()

    try:
        low_threshold, high_threshold = map(int, args.q.split(','))
    except ValueError:
        print("Error: -q argument should be two comma-separated integers")
        return
        
    # load splice distances from TSV file
    splice_distances = load_splice_distances(args.class_file)
    
    modified_count = 0  
    
    # prepare bam header entry 
    command_line = ' '.join(['python'] + sys.argv)
    description = 'Colored reads based on minimum distance to splice sites and added QS tag for quantized splice distance'
    
    with pysam.AlignmentFile(args.ibam, 'rb', threads=8) as infile:
        header = infile.header.to_dict()
        header = update_pg_header(header, command_line, description)
        
        # write to outfile 
        with pysam.AlignmentFile(args.obam, 'wb', header=header, threads=8) as outfile:

            for read in infile:
                read_info = splice_distances.get(read.query_name, {'min_splice_distance': -1, 'biotype': 'unknown'})
                min_splice_distance = read_info['min_splice_distance']
                biotype = read_info['biotype']
                
                color, qs_value = calculate_color_quantized(min_splice_distance, low_threshold, high_threshold)
                read.set_tag("YC", color)
                read.set_tag("QS", qs_value, value_type='i')
                read.set_tag("BT", biotype, value_type='Z')  # Add biotype as a tag
                
                outfile.write(read)
                modified_count += 1

    print(f"Sorting and indexing output BAM file")
    pysam.index(f"{args.obam}", threads=8)

    print(f"Total reads modified: {modified_count}")
    print(f"Quantized splice distance color: Red <= {low_threshold}, Yellow {low_threshold+1}-{high_threshold}, Blue > {high_threshold}, Grey for unknown")
    print(f"QS tag: 0 for distance <= {low_threshold}, 1 for {low_threshold} < distance <= {high_threshold}, 2 for distance > {high_threshold}, -1 for unknown")

if __name__ == "__main__":
    main()