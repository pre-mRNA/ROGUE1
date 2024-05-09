import argparse
import logging
import os
import pysam
import gffutils

def setup_logging(verbosity):
    level = logging.INFO if verbosity == 1 else logging.DEBUG if verbosity > 1 else logging.ERROR
    logging.basicConfig(level=level, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_bam_for_ends_and_junctions(bam_file):
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
    except FileNotFoundError:
        logging.error(f"BAM file not found: {bam_file}")
        raise SystemExit
    except ValueError:
        logging.error(f"Error reading BAM file: {bam_file}")
        raise SystemExit

    ends = {}
    junctions = {}
    read_count = 0

    try:
        for read in bam:
            if not read.is_unmapped and not read.is_secondary and not read.is_qcfail:
                chrom = bam.get_reference_name(read.reference_id)
                three_prime_end = read.reference_end if not read.is_reverse else read.reference_start
                ends.setdefault(chrom, []).append(three_prime_end)

                if read.cigar is not None:
                    pos = read.reference_start
                    for op, length in read.cigar:
                        if op == 3:
                            junctions.setdefault(chrom, []).append((pos, pos + length))
                        if op in [0, 2, 3]:
                            pos += length
            read_count += 1
    finally:
        bam.close()

    logging.info(f"Processed {read_count} reads from BAM file, found junctions in {len(junctions)} locations.")
    return ends, junctions

def load_gtf_db(gtf_file, force_rebuild):
    db_file = gtf_file + ".db"
    if force_rebuild or not os.path.exists(db_file):
        try:
            db = gffutils.create_db(gtf_file, db_file, force=True, keep_order=True)
            logging.info(f"Created new GTF database.")
        except Exception as e:
            logging.error(f"Error creating GTF database: {str(e)}")
            raise SystemExit
    else:
        db = gffutils.FeatureDB(db_file)
        logging.info(f"Loaded existing GTF database.")

    gene_count = sum(1 for _ in db.features_of_type('gene'))
    logging.info(f"GTF database loaded with {gene_count} genes.")
    return db

def filter_transcripts(ends, junctions, db):
    filtered_transcripts = {}
    total_transcripts = 0
    filtered_out = 0

    for gene in db.features_of_type('gene'):
        gene_id = gene.id
        for transcript in db.children(gene, featuretype='transcript'):
            total_transcripts += 1
            transcript_end = None
            transcript_junctions = []
            has_end_support = False
            has_junction_support = True

            for exon in sorted(db.children(transcript, featuretype='exon'), key=lambda x: x.start):
                if transcript.strand == '+':
                    transcript_end = exon.end
                else:
                    transcript_end = exon.start

                if transcript_end is not None:
                    last_end = exon.end if transcript.strand == '+' else exon.start
                    if last_end != transcript_end:  # Ensure it's a junction between exons
                        transcript_junctions.append((last_end, exon.start))

            if transcript_end:
                end_range = range(transcript_end - 20, transcript_end + 21)
                has_end_support = any(end in end_range for end in ends.get(gene.chrom, []))
            if transcript_junctions:
                has_junction_support = all(any(j == t_j for j in junctions.get(gene.chrom, [])) for t_j in transcript_junctions)

            if not (has_end_support and has_junction_support):
                filtered_out += 1
                logging.debug(f"Transcript {transcript.id} of gene {gene_id} filtered out due to lack of support.")
                continue

            filtered_transcripts[transcript.id] = transcript

    logging.info(f"Filtered out {filtered_out} out of {total_transcripts} transcripts due to no support.")
    return filtered_transcripts

def main(args):
    setup_logging(args.verbosity)
    ends, junctions = parse_bam_for_ends_and_junctions(args.bam_file)
    db = load_gtf_db(args.gtf_file, args.force_rebuild)
    filtered_transcripts = filter_transcripts(ends, junctions, db)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter transcripts based on BAM and GTF files.")
    parser.add_argument("--bam_file", required=True, help="Input BAM file path.")
    parser.add_argument("--gtf_file", required=True, help="Input GTF file path.")
    parser.add_argument("-v", "--verbosity", type=int, default=0, help="Increase output verbosity (0 = errors only, 1 = info, 2 = debug).")
    parser.add_argument("--force_rebuild", action='store_true', help="Force rebuild of GTF database if one exists.")
    args = parser.parse_args()
    main(args)
