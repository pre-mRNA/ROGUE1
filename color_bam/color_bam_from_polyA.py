import argparse
import pysam

def calculate_color(pt_value):
    
    # define color thresholds
    low_threshold = 30
    high_threshold = 120
    
    # clamping the pt_value between the thresholds
    if pt_value <= low_threshold:
        return "255,0,0"  # red
    elif pt_value >= high_threshold:
        return "0,0,255"  # blue
    else:

        # scale pt_value to a 0-1 range
        scaled_value = (pt_value - low_threshold) / (high_threshold - low_threshold)
        
        # interpolate between red (255,0,0) and blue (0,0,255)
        red = int(255 * (1 - scaled_value))
        blue = int(255 * scaled_value)
        return f"{red},0,{blue}"

def main():
    parser = argparse.ArgumentParser(description="Color reads in a BAM file based on PT tag values.")
    parser.add_argument("-ibam", required=True, help="Input BAM file")
    parser.add_argument("-obam", required=True, help="Output BAM file")

    args = parser.parse_args()
    
    modified_count = 0  # count of modified reads
    
    with pysam.AlignmentFile(args.ibam, 'rb', threads=8) as infile, \
         pysam.AlignmentFile(args.obam, 'wb', template=infile, threads=8) as outfile:
        for read in infile:
            pt_value = read.get_tag("PT", default=0)
            color = calculate_color(pt_value)
            read.set_tag("YC", color)
            outfile.write(read)
            modified_count += 1

    # print(f"Sorting and indexing output BAM file")
    
    # pysam.sort("-o", f"{args.obam}.sorted.bam", args.obam, threads=8)
    # pysam.index(f"{args.obam}.sorted.bam", threads=8)

    print(f"Total reads modified: {modified_count}")

if __name__ == "__main__":
    main()
