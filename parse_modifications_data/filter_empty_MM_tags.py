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
        
        logging.info("filtering only for 'A+a' modifications")
        
        mm_parts = mm_tag.split(';')
        cumulative_index = 0
        a_plus_a_indices = []
        a_plus_a_nums = []
        
        for part in mm_parts:
            if not part:
                continue
            mod_split = part.split(',')
            mod_name = mod_split[0]
            nums = mod_split[1:]
            num_count = len(nums)
            if mod_name.startswith('A+a'):
                a_plus_a_indices.extend(range(cumulative_index, cumulative_index + num_count))
                a_plus_a_nums.extend([int(n) for n in nums])
            cumulative_index += num_count
        
        if not a_plus_a_indices:
            return read
        
        ml_values = list(ml_tag)
        
        if len(ml_values) < cumulative_index:
            logging.error(f"ml tag has fewer values ({len(ml_values)}) than expected ({cumulative_index})")
            return read
        
        a_plus_a_ml_values = [ml_values[i] for i in a_plus_a_indices]
        filtered_a_plus_a_nums = []
        new_ml_values_filtered = []
        cumulative_sum = 0
        removed_count = 0
        
        for mod, ml in zip(a_plus_a_nums, a_plus_a_ml_values):
            if ml == 255:
                if cumulative_sum > 0 or removed_count > 0:
                    filtered_a_plus_a_nums.append(mod + cumulative_sum + removed_count)
                    cumulative_sum = 0
                    removed_count = 0
                else:
                    filtered_a_plus_a_nums.append(mod)
                new_ml_values_filtered.append(ml)
            else:
                cumulative_sum += mod
                removed_count += 1
        
        if filtered_a_plus_a_nums:
            new_mm_tag = f"A+a?," + ",".join(map(str, filtered_a_plus_a_nums)) + ";"
            new_ml_tag = array.array('B', new_ml_values_filtered)
        else:
            new_mm_parts = [part for part in mm_parts if not part.startswith('A+a')]
            new_mm_tag = ';'.join(new_mm_parts) if new_mm_parts else ''
            new_ml_values = [ml for i, ml in enumerate(ml_values) if i not in a_plus_a_indices]
            new_ml_tag = array.array('B', new_ml_values)
        
        if new_mm_tag:
            read.set_tag('MM', new_mm_tag)
        else:
            read.remove_tag('MM')
        
        if new_ml_tag:
            read.set_tag('ML', new_ml_tag)
        else:
            read.remove_tag('ML')
        
    return read

def main():
    parser = argparse.ArgumentParser(description="Process BAM file and modify MM/ML tags.")
    parser.add_argument("-ibam", required=True, help="input bam file")
    parser.add_argument("-obam", required=True, help="output bam file")
    parser.add_argument("-t", "--threads", type=int, default=8, help="number of threads to use")

    args = parser.parse_args()

    setup_logging()

    logging.info(f"starting BAM processing with {args.threads} threads")
    logging.info(f"input BAM file: {args.ibam}")
    logging.info(f"output BAM file: {args.obam}")

    if not os.path.exists(args.ibam):
        logging.error(f"input BAM file does not exist: {args.ibam}")
        sys.exit(1)

    modified_count = 0
    total_count = 0

    command_line = ' '.join(['python'] + sys.argv)
    description = 'modified MM and ML tags: filtered for A+a with ML=255'

    try:
        with pysam.AlignmentFile(args.ibam, 'rb', threads=args.threads) as infile:
            header = infile.header.to_dict()
            header = update_pg_header(header, command_line, description)

            with pysam.AlignmentFile(args.obam, 'wb', header=header, threads=args.threads) as outfile:
                for read in infile:
                    total_count += 1
                    modified_read = process_read(read)
                    outfile.write(modified_read)
                    if modified_read.has_tag('MM') and modified_read.has_tag('ML'):
                        modified_count += 1

                    if total_count % 100000 == 0:
                        logging.info(f"processed {total_count} reads, modified {modified_count}")

        logging.info(f"indexing output BAM file: {args.obam}")
        pysam.index(args.obam, threads=args.threads)

        logging.info(f"total reads processed: {total_count}")
        logging.info(f"reads with modified MM/ML tags: {modified_count}")

    except Exception as e:
        logging.error(f"an error occurred during processing: {str(e)}")
        logging.error(f"error details: {str(e)}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
