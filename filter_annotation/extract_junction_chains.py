import argparse
import logging
import pysam

def setup_logging(verbosity):
    level = logging.INFO if verbosity == 1 else logging.DEBUG if verbosity > 1 else logging.ERROR
    logging.basicConfig(level=level, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_bam_for_junctions(bam_file):
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
    except FileNotFoundError:
        logging.error(f"BAM file not found: {bam_file}")
        raise SystemExit
    except ValueError:
        logging.error(f"Error reading BAM file: {bam_file}")
        raise SystemExit

    junctions = {}
    read_count = 0

    try:
        for read in bam:
            if not read.is_unmapped and not read.is_secondary and not read.is_qcfail:
                chrom = bam.get_reference_name(read.reference_id)
                pos = read.reference_start

                for op, length in read.cigar:
                    if op == 3:  # N CIGAR operation 
                        junction_start = pos
                        junction_end = pos + length
                        junctions.setdefault(chrom, []).append((junction_start, junction_end))
                    if op in [0, 1, 2, 3]:  # M, I, D, N operations affect the reference position
                        pos += length
            read_count += 1
    finally:
        bam.close()

    logging.info(f"Processed {read_count} reads from BAM file, identified splice junctions in {len(junctions)} locations.")
    return junctions

def main(args):
    setup_logging(args.verbosity)
    junctions = parse_bam_for_junctions(args.bam_file)

    # Optionally, save the junctions to a file
    if args.output_file:
        with open(args.output_file, 'w') as file:
            for chrom, junction_list in junctions.items():
                for start, end in junction_list:
                    file.write(f"{chrom}\t{start}\t{end}\n")
        logging.info(f"Splice junctions saved to {args.output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract splice junctions from a BAM file.")
    parser.add_argument("--bam_file", required=True, help="Input BAM file path.")
    parser.add_argument("--output_file", help="Output file path for saving junction coordinates.")
    parser.add_argument("-v", "--verbosity", type=int, default=0, help="Increase output verbosity.")
    args = parser.parse_args()
    main(args)
