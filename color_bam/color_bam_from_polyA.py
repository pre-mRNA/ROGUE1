import os 
import sys 
import argparse
import pysam
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
    
# sliding scale for poly(A) color 
def calculate_color_sliding(pt_value):
    
    low_threshold = 30
    high_threshold = 120
    gray_color = "128,128,128"  # RGB for grey

    if pt_value == -1:
        return gray_color
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

# quantize poly(A) color for nascent and mature thresholds
# returns YC tag and QP tag 
def calculate_color_quantized(pt_value, low_threshold, high_threshold):
    
    gray_color = "128,128,128"  # RGB for grey

    if pt_value == -1:
        return gray_color, 1
    if pt_value <= low_threshold:
        return "255,0,0", 0  # red
    elif pt_value >= high_threshold:    
        return "0,0,255", 1  # blue
    else:
        return gray_color, 2  # grey for values between low and high thresholds

def main():
    parser = argparse.ArgumentParser(description="Color reads in a BAM file based on PT tag values.")
    parser.add_argument("-ibam", required=True, help="Input BAM file")
    parser.add_argument("-obam", required=True, help="Output BAM file")
    parser.add_argument("-q", type=str, help="Comma-separated low and high thresholds for quantized poly(A) color")

    args = parser.parse_args()

    # check if quantized thresholds should be used 
    use_quantized_polya_thresholds = False
    if args.q:
        try:
            low_threshold, high_threshold = map(int, args.q.split(','))
            use_quantized_polya_thresholds = True
        except ValueError:
            print("Error: -q argument should be two comma-separated integers")
            return
        
    modified_count = 0  # count of modified reads
    
    # prepare bam header entry 
    command_line = ' '.join(['python'] + sys.argv)
    description = 'Colored reads based on PT tag values and added QP tag for quantized polyA'
    
    with pysam.AlignmentFile(args.ibam, 'rb', threads=8) as infile:
        header = infile.header.to_dict()
        header = update_pg_header(header, command_line, description)
        
        # write to outfile 
        with pysam.AlignmentFile(args.obam, 'wb', header=header, threads=8) as outfile:

            for read in infile:
                try:
                    pt_value = read.get_tag("pt")
                except KeyError:
                    pt_value = -1  # default PT value if not present
                
                if use_quantized_polya_thresholds:
                    color, qp_value = calculate_color_quantized(pt_value, low_threshold, high_threshold)
                    read.set_tag("YC", color)

                    # set quantized polyA QP tag 
                    if qp_value is not None:
                        read.set_tag("QP", qp_value, value_type='i')
                else:
                    color = calculate_color_sliding(pt_value)
                
                read.set_tag("YC", color)
                outfile.write(read)
                modified_count += 1

    print(f"Sorting and indexing output BAM file")
    pysam.index(f"{args.obam}", threads=8)

    print(f"Total reads modified: {modified_count}")
    if use_quantized_polya_thresholds:
        print(f"Quantized polyA color: Red <= {low_threshold}, {low_threshold+1}-{high_threshold-1} Grey, Blue >= {high_threshold}")
        print(f"QP tag: 0 for polyA <= {low_threshold}, 1 for polyA >= {high_threshold}, 2 for {low_threshold} < polyA < {high_threshold}")
    else:
        print("Scaled polyA color: Red <= 30, 31-119 intermediate, Blue >= 120")


if __name__ == "__main__":
    main()
