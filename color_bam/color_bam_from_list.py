import argparse
import pysam
from pysam import AlignedSegment

def validate_color(color):
    parts = color.split(',')
    if len(parts) != 3:
        raise ValueError("Colour must be a comma-separated string of three integers between 0-255.")
    for part in parts:
        if not (0 <= int(part) <= 255):
            raise ValueError("Each RGB value must be between 0 and 255.")
    return color

def read_read_ids(file_path):
    with open(file_path, 'r') as f:
        return {line.strip() for line in f if line.strip()}

def main():
    parser = argparse.ArgumentParser(description="Set a specific RGB color for a list of read IDs in a BAM file.")
    parser.add_argument("-ibam", required=True, help="Input BAM file")
    parser.add_argument("-obam", required=True, help="Output BAM file")
    parser.add_argument("-rids", required=True, help="File containing list of read IDs")
    parser.add_argument("-color", required=True, help="RGB color to set for the specified read IDs")

    args = parser.parse_args()

    read_ids = read_read_ids(args.rids)
    color = validate_color(args.color)
    
    modified_count = 0  # count modified reads
    
    with pysam.AlignmentFile(args.ibam, 'rb', threads=8) as infile, \
         pysam.AlignmentFile(args.obam, 'wb', template=infile, threads=8) as outfile:
        for read in infile:
            if read.query_name in read_ids:
                read.set_tag("YC", color)
                outfile.write(read)
                modified_count += 1
            else:
                outfile.write(read)

    print(f"Sorting and indexing output bam file")
    
    # pysam.sort("-o", f"{args.obam}.sorted.bam", args.obam, threads=8)
    # pysam.index(f"{args.obam}.sorted.bam", threads=8)

    print(f"Total reads modified: {modified_count}")

if __name__ == "__main__":
    main()
