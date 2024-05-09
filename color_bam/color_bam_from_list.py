import argparse
import pysam
from pysam import AlignedSegment

# validate the rgb color provided bt the user 
def validate_color(color):
    parts = color.split(',')
    if len(parts) != 3:
        raise ValueError("colour must be a comma-separated string of three integers between 0-255.")
    for part in parts:
        if not (0 <= int(part) <= 255):
            raise ValueError("Each RGB value must be between 0 and 255.")
    return color

# read list of read IDs 
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
    
    with pysam.AlignmentFile(args.ibam, 'rb') as infile, pysam.AlignmentFile(args.obam, 'wb', template=infile) as outfile:
        for read in infile:
            if read.query_name in read_ids:
                read.set_tag("YC", color)
            outfile.write(read)

    pysam.sort("-o", f"{args.obam}.sorted.bam", args.obam)
    pysam.index(f"{args.obam}.sorted.bam")

if __name__ == "__main__":
    main()
