import argparse
import pysam

def calculate_ends(read):
    """
    Calculate the last mapped position of a read.
    For reverse strand reads, this is the leftmost position (start).
    For forward strand reads, this is the rightmost position (end).
    
    :param read: A pysam AlignedSegment object
    :return: The last mapped position as an integer
    """
    if read.is_reverse:
        return -read.reference_end , -read.reference_start - 1 
    else:
        return read.reference_start + 1, read.reference_end

def calculate_aligned_bases(read):
    """
    Calculate the total number of aligned bases for a read.
    
    :param read: A pysam AlignedSegment object
    :return: The total number of aligned bases as an integer
    """
    return sum([length for operation, length in read.cigartuples if operation in [0, 7, 8]])

def main():
    parser = argparse.ArgumentParser(description="Add RE and AL tags to reads in a BAM file.")
    parser.add_argument("-ibam", required=True, help="Input BAM file")
    parser.add_argument("-obam", required=True, help="Output BAM file")

    args = parser.parse_args()
    
    modified_count = 0  # count of modified reads

    with pysam.AlignmentFile(args.ibam, 'rb', threads=48) as infile, \
         pysam.AlignmentFile(args.obam, 'wb', template=infile, threads=48) as outfile:
        for read in infile:

            first_position, last_position = calculate_ends(read)
            
            # add RE tag with the last mapped position
            # RE tag allows sorting plus strand alignments by read end position in IGV, which is useful for visualizing pol II progression
            read.set_tag("RE", last_position, value_type='i')

            # add RS tag with first mapped position
            # RS tag allows sorting for minus strand alignments by pol II position 
            read.set_tag("RS", first_position, value_type='i')

            # Add AL tag with the total number of aligned bases     
            aligned_bases = calculate_aligned_bases(read)
            read.set_tag("AL", aligned_bases, value_type='i')
            
            outfile.write(read)
            modified_count += 1

    print(f"Indexing output BAM file")
    
    # Index BAM file
    pysam.index(f"{args.obam}", threads=8)

    print(f"Total reads modified: {modified_count}")

if __name__ == "__main__":
    main()