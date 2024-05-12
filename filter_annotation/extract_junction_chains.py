import logging
import pandas as pd
import pysam

def setup_logging(verbosity=0):
    levels = {0: logging.ERROR, 1: logging.INFO, 2: logging.DEBUG}
    logging.basicConfig(level=levels.get(verbosity, logging.ERROR),
                        format='%(asctime)s - %(levelname)s - %(message)s')

def parse_bam_for_junctions(bam_file, verbosity=0):
    setup_logging(verbosity)
    
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
    except FileNotFoundError:
        logging.error(f"BAM file not found: {bam_file}")
        return pd.DataFrame()  # Return empty DataFrame on error
    except ValueError:
        logging.error(f"Error reading BAM file: {bam_file}")
        return pd.DataFrame()  # Return empty DataFrame on error

    junctions = []
    read_count = 0

    try:
        for read in bam:
            if not read.is_unmapped and not read.is_secondary and not read.is_qcfail:
                chrom = bam.get_reference_name(read.reference_id)
                pos = read.reference_start

                for op, length in read.cigar:
                    if op == 3:  # N CIGAR operation (splicing)
                        junction_start = pos
                        junction_end = pos + length
                        junctions.append([chrom, junction_start, junction_end])
                    if op in [0, 1, 2, 3]:  # M, I, D, N operations affect the reference position
                        pos += length
            read_count += 1
    finally:
        bam.close()

    logging.info(f"Processed {read_count} reads from BAM file, identified splice junctions in {len(junctions)} chromosomal locations.")
    
    if junctions:
        return pd.DataFrame(junctions, columns=['chromosome', 'start', 'end'])
    else:
        return pd.DataFrame(columns=['chromosome', 'start', 'end'])

if __name__ == "__main__":
    junctions_df = parse_bam_for_junctions("path_to_bam_file.bam", verbosity=1)
    print(junctions_df)
