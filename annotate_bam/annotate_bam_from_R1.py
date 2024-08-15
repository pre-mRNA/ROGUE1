import argparse
import pysam
import pandas as pd
import os
import sys

sys.path.append("/home/150/as7425/R1/color_bam")
from bam_header_utils import update_pg_header

# color mappings
FEATURE_COLOR_MAP = {
    'exon': '51,153,255',  # blue
    'intron': '255,0,112',  # pink
    'DOG': '0,255,0',  # green
    'other': '128,128,128'  # grey
}

CLASSIFICATION_COLOR_MAP = {
    'intron': '255,0,112',
    'spliced': '51,153,255',
    'ambiguous': '255,102,0',
    'discard': '255,102,0',
    'intron_retained': '255,0,112',
    'fully-unspliced': '255,0,112',
    'fully_unspliced': '255,0,112',
    'partially_spliced': '250,0,255'
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
    if isinstance(color, float) and pd.isna(color):
        raise ValueError("color value is missing or not a string.")
    parts = color.split(',')
    if len(parts) != 3:
        raise ValueError("color must be a comma-separated string of three integers between 0-255.")
    for part in parts:
        if not (0 <= int(part) <= 255):
            raise ValueError("each RGB value must be between 0 and 255.")
    return color

def load_r1_data(file_path):
    df = pd.read_csv(file_path, sep='\t')
    return df.set_index('read_id').to_dict('index')

def calculate_color_quantized(splice_distance, low_threshold, high_threshold):
    if pd.isna(splice_distance) or splice_distance == -1:
        return "128,128,128", -1  # grey for unknown
    elif splice_distance <= low_threshold:
        return "255,0,0", 0  # red
    elif splice_distance <= high_threshold:
        return "255,255,0", 1  # yellow
    else:
        return "0,0,255", 2  # blue

def calculate_polya_color(pt_value, low_threshold, high_threshold):
    if pd.isna(pt_value) or pt_value == -1:
        return "128,128,128", 2  # grey for unknown
    elif pt_value <= low_threshold:
        return "255,0,0", 0  # red
    elif pt_value >= high_threshold:
        return "0,0,255", 1  # blue
    else:
        return "128,128,128", 2  # grey for values between thresholds

def main():
    parser = argparse.ArgumentParser(description="Annotate BAM file with R1.py output data")
    parser.add_argument("-ibam", required=True, help="Input BAM file")
    parser.add_argument("-obam", required=True, help="Output BAM file")
    parser.add_argument("-r1_table", required=True, help="R1.py output table")
    parser.add_argument("-q", type=str, required=True, help="Comma-separated low and high thresholds for quantized splice distance and polyA color")
    args = parser.parse_args()

    low_threshold, high_threshold = map(int, args.q.split(','))
    
    r1_data = load_r1_data(args.r1_table)
    
    command_line = ' '.join(['python'] + sys.argv)
    description = 'Annotated BAM file with R1.py output data'
    
    with pysam.AlignmentFile(args.ibam, 'rb', threads=8) as infile:
        header = infile.header.to_dict()
        header = update_pg_header(header, command_line, description)
        
        with pysam.AlignmentFile(args.obam, 'wb', header=header, threads=8) as outfile:
            for read in infile:
                if read.query_name in r1_data:
                    data = r1_data[read.query_name]
                    
                    # gene id
                    if 'gene_id' in data and not pd.isna(data['gene_id']):
                        read.set_tag("GI", data['gene_id'], value_type='Z')
                    
                    # 3' feature annotation
                    if 'three_prime_feature' in data and not pd.isna(data['three_prime_feature']):
                        feature_type = data['three_prime_feature']
                        color = FEATURE_COLOR_MAP.get(feature_type, FEATURE_COLOR_MAP['other'])
                        read.set_tag("YC", color)
                        read.set_tag("EF", feature_type, value_type='Z')
                    
                    # splice distance annotation
                    if 'splice_donors_downstream_distance' in data and 'splice_donors_upstream_distance' in data:
                        splice_distance = min(abs(data.get('splice_donors_downstream_distance', float('inf'))),
                                              abs(data.get('splice_donors_upstream_distance', float('inf'))))
                        color, qs_value = calculate_color_quantized(splice_distance, low_threshold, high_threshold)
                        read.set_tag("QS", qs_value, value_type='i')
                
                    # biotype annotation
                    if 'biotype' in data and not pd.isna(data['biotype']):
                        read.set_tag("BT", data['biotype'], value_type='Z')
                    
                    # splicing classification annotation
                    if 'canonical_acceptor_splicing_status' in data and not pd.isna(data['canonical_acceptor_splicing_status']):
                        splicing_class = data['canonical_acceptor_splicing_status']
                        class_color = CLASSIFICATION_COLOR_MAP.get(splicing_class, "128,128,128")
                        read.set_tag("SC", splicing_class, value_type='Z')
                        read.set_tag("YC", class_color)
                    
                    # polyA annotation
                    if 'polya_length' in data and not pd.isna(data['polya_length']):
                        pt_value = data['polya_length']
                        polya_color, qp_value = calculate_polya_color(pt_value, low_threshold, high_threshold)
                        read.set_tag("QP", qp_value, value_type='i')
                    
                    # additional fields
                    if 'read_length' in data and not pd.isna(data['read_length']):
                        read.set_tag("RL", int(data['read_length']), value_type='i')
                    if 'alignment_length' in data and not pd.isna(data['alignment_length']):
                        read.set_tag("AL", int(data['alignment_length']), value_type='i')
                    if 'splice_count' in data and not pd.isna(data['splice_count']):
                        read.set_tag("SP", int(data['splice_count']), value_type='i')
                    if 'gene_base_overlap' in data and not pd.isna(data['gene_base_overlap']):
                        read.set_tag("GO", int(data['gene_base_overlap']), value_type='i')
                    if 'exon_base_overlap' in data and not pd.isna(data['exon_base_overlap']):
                        read.set_tag("EO", int(data['exon_base_overlap']), value_type='i')
                    if 'intron_base_overlap' in data and not pd.isna(data['intron_base_overlap']):
                        read.set_tag("IO", int(data['intron_base_overlap']), value_type='i')
                    if 'dog_base_overlap' in data and not pd.isna(data['dog_base_overlap']):
                        read.set_tag("DO", int(data['dog_base_overlap']), value_type='i')
                    if 'gene_name' in data and not pd.isna(data['gene_name']):
                        read.set_tag("GN", data['gene_name'], value_type='Z')
                    if 'exon_count' in data and not pd.isna(data['exon_count']):
                        read.set_tag("EC", int(data['exon_count']), value_type='i')
                
                outfile.write(read)

    print(f"Indexing output BAM file")
    pysam.index(f"{args.obam}", threads=8)

    print("Added tags (when information is available):")
    print("  GI: Gene ID")
    print("  YC: RGB color based on priority (splicing class > 3' feature > splice distance > polyA)")
    print("  EF: 3' end feature")
    print("  QS: Quantized splice distance")
    print("  BT: Biotype")
    print("  SC: Splicing classification (canonical_acceptor_splicing_status)")
    print("  QP: Quantized polyA length")
    print("  RL: Read length")
    print("  AL: Alignment length")
    print("  SP: Splice count")
    print("  GO: Gene base overlap")
    print("  EO: Exon base overlap")
    print("  IO: Intron base overlap")
    print("  DO: Downstream of gene (DOG) base overlap")
    print("  GN: Gene name")
    print("  EC: Exon count")
    print(f"Quantized thresholds: <= {low_threshold}, {low_threshold+1}-{high_threshold}, > {high_threshold}")

if __name__ == "__main__":
    main()