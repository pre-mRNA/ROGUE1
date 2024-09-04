import argparse
import pysam
import pandas as pd
import os
import sys
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'color_bam'))
from bam_header_utils import update_pg_header

FEATURE_COLOR_MAP = {
    'exon': '51,153,255', 'intron': '255,0,112', 'DOG': '0,255,0', 'other': '128,128,128'
}

CLASSIFICATION_COLOR_MAP = {
    'intron': '255,0,112', 'spliced': '51,153,255', 'ambiguous': '255,102,0',
    'discard': '255,102,0', 'intron_retained': '255,0,112',
    'fully-unspliced': '255,0,112', 'fully_unspliced': '255,0,112',
    'partially_spliced': '250,0,255'
}

def load_r1_data(file_path):
    logging.info("Loading R1 data...")
    df = pd.read_csv(file_path, sep='\t', low_memory=False)
    r1_data = df.set_index('read_id').to_dict('index')
    logging.info(f"Loaded {len(r1_data)} entries from R1 data.")
    return r1_data

def calculate_color_quantized(splice_distance, low_threshold, high_threshold):
    if pd.isna(splice_distance) or splice_distance == -1:
        return "128,128,128", -1
    elif splice_distance <= low_threshold:
        return "255,0,0", 0
    elif splice_distance <= high_threshold:
        return "255,255,0", 1
    else:
        return "0,0,255", 2

def calculate_polya_color(pt_value, low_threshold, high_threshold):
    if pd.isna(pt_value) or pt_value == -1:
        return "128,128,128", 2
    elif pt_value <= low_threshold:
        return "255,0,0", 0
    elif pt_value >= high_threshold:
        return "0,0,255", 1
    else:
        return "128,128,128", 2

def calculate_ends(read):
    if read.is_reverse:
        return -read.reference_end, -read.reference_start - 1 
    else:
        return read.reference_start + 1, read.reference_end

def calculate_aligned_bases(read):
    return sum([length for operation, length in read.cigartuples if operation in [0, 7, 8]])

def process_read(read, r1_data, low_threshold, high_threshold):
    added_tags = set()
    try:
        if read.query_name in r1_data:
            data = r1_data[read.query_name]
            
            if 'gene_id' in data and not pd.isna(data['gene_id']):
                read.set_tag("GI", str(data['gene_id']), value_type='Z')
                added_tags.add("GI")
            
            if 'three_prime_feature' in data and not pd.isna(data['three_prime_feature']):
                feature_type = data['three_prime_feature']
                color = FEATURE_COLOR_MAP.get(feature_type, FEATURE_COLOR_MAP['other'])
                read.set_tag("YC", color)
                read.set_tag("EF", feature_type, value_type='Z')
                added_tags.add("YC")
                added_tags.add("EF")
            
            if 'splice_donors_downstream_distance' in data and 'splice_donors_upstream_distance' in data:
                splice_distance = min(abs(data.get('splice_donors_downstream_distance', float('inf'))),
                                      abs(data.get('splice_donors_upstream_distance', float('inf'))))
                color, qs_value = calculate_color_quantized(splice_distance, low_threshold, high_threshold)
                read.set_tag("QS", qs_value, value_type='i')
                added_tags.add("QS")
        
            if 'biotype' in data and not pd.isna(data['biotype']):
                read.set_tag("BT", data['biotype'], value_type='Z')
                added_tags.add("BT")
            
            if 'canonical_acceptor_splicing_status' in data and not pd.isna(data['canonical_acceptor_splicing_status']):
                splicing_class = data['canonical_acceptor_splicing_status']
                class_color = CLASSIFICATION_COLOR_MAP.get(splicing_class, "128,128,128")
                read.set_tag("SC", splicing_class, value_type='Z')
                read.set_tag("YC", class_color)
                added_tags.add("SC")
                added_tags.add("YC")
            
            if 'polya_length' in data and not pd.isna(data['polya_length']):
                pt_value = data['polya_length']
                polya_color, qp_value = calculate_polya_color(pt_value, low_threshold, high_threshold)
                read.set_tag("QP", qp_value, value_type='i')
                added_tags.add("QP")
            
            for field, tag in [('read_length', 'RL'), ('alignment_length', 'AL'), ('splice_count', 'SP'),
                               ('gene_base_overlap', 'GO'), ('exon_base_overlap', 'EO'), ('intron_base_overlap', 'IO'),
                               ('dog_base_overlap', 'DO'), ('exon_count', 'EC')]:
                if field in data and not pd.isna(data[field]):
                    read.set_tag(tag, int(data[field]), value_type='i')
                    added_tags.add(tag)
            
            if 'gene_name' in data and not pd.isna(data['gene_name']):
                read.set_tag("GN", str(data['gene_name']), value_type='Z')
                added_tags.add("GN")
        
        first_position, last_position = calculate_ends(read)
        read.set_tag("RE", last_position, value_type='i')
        read.set_tag("RS", first_position, value_type='i')
        aligned_bases = calculate_aligned_bases(read)
        read.set_tag("AL", aligned_bases, value_type='i')
        tt_value = f"{last_position}.{first_position}"
        read.set_tag("TT", tt_value, value_type='Z')
        ta_value = f"{last_position}.{aligned_bases}"
        read.set_tag("TA", ta_value, value_type='Z')
        added_tags.update(["RE", "RS", "AL", "TT", "TA"])
    except Exception as e:
        logging.error(f"Error processing read {read.query_name}: {str(e)}")
    
    return read, added_tags

def main():
    parser = argparse.ArgumentParser(description="Annotate BAM file with R1.py output data")
    parser.add_argument("-ibam", required=True, help="Input BAM file")
    parser.add_argument("-obam", required=True, help="Output BAM file")
    parser.add_argument("-r1_table", required=True, help="R1.py output table")
    parser.add_argument("-q", type=str, required=True, help="Comma-separated low and high thresholds for quantized splice distance and polyA color")
    parser.add_argument("-t", type=int, default=8, help="Number of threads to use for BAM operations")
    args = parser.parse_args()

    low_threshold, high_threshold = map(int, args.q.split(','))
    
    r1_data = load_r1_data(args.r1_table)
    
    command_line = ' '.join(['python'] + sys.argv)
    description = 'Annotated BAM file with R1.py output data'
    
    logging.info(f"Processing BAM file using {args.t} threads")
    try:
        with pysam.AlignmentFile(args.ibam, 'rb', threads=args.t) as infile:
            header = infile.header.to_dict()
            header = update_pg_header(header, command_line, description)
            
            with pysam.AlignmentFile(args.obam, 'wb', header=header, threads=args.t) as outfile:
                processed_reads = 0
                all_added_tags = set()
                for read in infile:
                    processed_read, added_tags = process_read(read, r1_data, low_threshold, high_threshold)
                    outfile.write(processed_read)
                    all_added_tags.update(added_tags)
                    processed_reads += 1
                    if processed_reads % 1000000 == 0:
                        logging.info(f"Processed {processed_reads} reads")

        logging.info(f"Total reads processed: {processed_reads}")
        logging.info(f"Indexing output BAM file")
        pysam.index(f"{args.obam}", threads=args.t)

    except Exception as e:
        logging.error(f"An error occurred: {str(e)}")
        raise

    logging.info("Processing completed successfully")
    
    print("Added tags (when information is available):")
    tag_descriptions = {
        "GI": "Gene ID",
        "YC": "RGB color based on priority (splicing class > 3' feature > splice distance > polyA)",
        "EF": "3' end feature",
        "QS": "Quantized splice distance",
        "BT": "Biotype",
        "SC": "Splicing classification (canonical_acceptor_splicing_status)",
        "QP": "Quantized polyA length",
        "RL": "Read length",
        "AL": "Alignment length",
        "SP": "Splice count",
        "GO": "Gene base overlap",
        "EO": "Exon base overlap",
        "IO": "Intron base overlap",
        "DO": "Downstream of gene (DOG) base overlap",
        "GN": "Gene name",
        "EC": "Exon count",
        "RE": "Read end position",
        "RS": "Read start position",
        "TT": "Total transcript position (end.start)",
        "TA": "Total transcript end position and aligned length (end.aligned_length)"
    }
    for tag in sorted(all_added_tags):
        print(f"  {tag}: {tag_descriptions.get(tag, 'Unknown tag')}")
    print(f"Quantized thresholds: <= {low_threshold}, {low_threshold+1}-{high_threshold}, > {high_threshold}")

if __name__ == "__main__":
    main()