import os
import argparse
import sys
import pysam

def main(args):

    os.makedirs(args.output_folder, exist_ok=True)

    output_base_path = os.path.join(args.output_folder, args.sample_name)

    try:
        bam_file = pysam.AlignmentFile(args.input_bam_path, "rb")
    except FileNotFoundError:
        print(f"Error: The file {args.input_bam_path} does not exist.")
        sys.exit(1)
    except IOError:
        print(f"Error: Could not open {args.input_bam_path}. Ensure it is a valid BAM file.")
        sys.exit(1)

    # create output BAM files for spliced and intron-containing reads
    spliced_bam = pysam.AlignmentFile(f"{output_base_path}_spliced.bam", "wb", header=bam_file.header)
    intron_bam = pysam.AlignmentFile(f"{output_base_path}_intron.bam", "wb", header=bam_file.header)

    # statistics 
    read_count = 0
    spliced_count = 0
    intron_count = 0
    ambiguous_count = 0
    unclassified_count = 0

    # process each read in the input BAM
    for read in bam_file:
        read_count += 1
        yc_tag = read.get_tag('YC') if read.has_tag('YC') else None

        if yc_tag == "51,153,255":  # tag for spliced reads
            spliced_bam.write(read)
            spliced_count += 1
        elif yc_tag == "255,0,112":  # tag for intron-retained reads
            intron_bam.write(read)
            intron_count += 1
        elif yc_tag == "255,102,0":  # tag for ambiguous splicing
            ambiguous_count += 1
        else:
            unclassified_count += 1

    # close files
    bam_file.close()
    spliced_bam.close()
    intron_bam.close()

    # index the BAM files
    pysam.index(f"{output_base_path}_spliced.bam")
    pysam.index(f"{output_base_path}_intron.bam")

    # print statistics
    print(f"Processed {read_count} reads in total.")
    print(f"{spliced_count} reads classified as spliced.")
    print(f"{intron_count} reads classified as intron-retained.")
    print(f"{ambiguous_count} reads classified as ambiguous (not written out).")
    print(f"{unclassified_count} reads with unclassified or missing YC tags.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Split BAM files based on YC tag classifications.')
    parser.add_argument('-i', '--input_bam_path', required=True, type=str, help='Path to the input BAM file')
    parser.add_argument('-o', '--output_folder', required=True, type=str, help='Path to the output folder where BAM files will be stored')
    parser.add_argument('-s', '--sample_name', required=True, type=str, help='Sample name to use as prefix for output files')
    args = parser.parse_args()
    
    main(args)
