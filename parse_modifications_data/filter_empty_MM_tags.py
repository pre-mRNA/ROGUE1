import argparse
import pysam
import os
import sys
import logging
import array
import re

sys.path.append("/home/150/as7425/R1/color_bam/")
from bam_header_utils import update_pg_header

def setup_logging():
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

def get_version():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(script_dir)
    version_file = os.path.join(parent_dir, 'version.txt')
    try:
        with open(version_file, 'r') as f:
            return f.read().strip()
    except FileNotFoundError:
        return "unknown"

def parse_mm_tag(mm_tag):
    parts = mm_tag.split(';')
    parsed_tags = []
    for part in parts:
        if part:
            match = re.match(r'([ACGTUN])([-+])([a-z]+|\d+)([.?]?)((?:,\d+)*)', part)
            if match:
                base, strand, mod, flag, positions = match.groups()
                positions = [int(p) for p in positions.split(',') if p]
                parsed_tags.append((base, strand, mod, flag, positions))
    return parsed_tags

def process_read(read):
    if read.has_tag('MM') and read.has_tag('ML'):
        mm_tag = read.get_tag('MM')
        ml_tag = read.get_tag('ML')
        
        mm_parts = mm_tag.split(';')
        for part in mm_parts:
            if part.startswith('A+a'):
                base_mods = list(map(int, part.split(',')[1:]))
                break
        else:
            return read  
        
        ml_values = list(ml_tag)
        
        new_base_mods = []
        new_ml_values = []
        cumulative_sum = 0
        removed_count = 0
        
        for i, (mod, ml) in enumerate(zip(base_mods, ml_values)):
            if ml == 255:
                if cumulative_sum > 0 or removed_count > 0:
                    new_base_mods.append(mod + cumulative_sum + removed_count)
                    cumulative_sum = 0
                    removed_count = 0
                else:
                    new_base_mods.append(mod)
                new_ml_values.append(ml)
            else:
                cumulative_sum += mod
                removed_count += 1
        
        new_mm_tag = f"A+a?,{','.join(map(str, new_base_mods))};"
        read.set_tag('MM', new_mm_tag)
        
        new_ml_tag = array.array('B', new_ml_values)
        read.set_tag('ML', new_ml_tag)
    
    return read

def main():
    parser = argparse.ArgumentParser(description="Process BAM file and modify MM/ML tags.")
    parser.add_argument("-ibam", required=True, help="Input BAM file")
    parser.add_argument("-obam", required=True, help="Output BAM file")
    parser.add_argument("-t", "--threads", type=int, default=8, help="Number of threads to use")

    args = parser.parse_args()

    setup_logging()

    logging.info(f"Starting BAM processing with {args.threads} threads")
    logging.info(f"Input BAM file: {args.ibam}")
    logging.info(f"Output BAM file: {args.obam}")

    if not os.path.exists(args.ibam):
        logging.error(f"Input BAM file does not exist: {args.ibam}")
        sys.exit(1)

    modified_count = 0
    total_count = 0

    # update bam header with version and command line
    command_line = ' '.join(['python'] + sys.argv)
    description = 'Modified MM and ML tags: filtered for ML=255 and adjusted MM values'

    try:
        with pysam.AlignmentFile(args.ibam, 'rb', threads=args.threads) as infile:
            header = infile.header.to_dict()
            header = update_pg_header(header, command_line, description)

            with pysam.AlignmentFile(args.obam, 'wb', header=header, threads=args.threads) as outfile:
                for read in infile:
                    total_count += 1
                    modified_read = process_read(read)
                    outfile.write(modified_read)
                    if read.has_tag('MM') and read.has_tag('ML'):
                        modified_count += 1

                    if total_count % 100000 == 0:
                        logging.info(f"Processed {total_count} reads, modified {modified_count}")

        logging.info(f"Indexing output BAM file: {args.obam}")
        pysam.index(args.obam, threads=args.threads)

        logging.info(f"Total reads processed: {total_count}")
        logging.info(f"Reads with modified MM/ML tags: {modified_count}")

    except Exception as e:
        logging.error(f"An error occurred during processing: {str(e)}")
        logging.error(f"Error details: {str(e)}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()