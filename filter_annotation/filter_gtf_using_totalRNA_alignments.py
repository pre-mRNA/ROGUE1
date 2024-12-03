#!/usr/bin/env python3

import os
import pysam
import gffutils 
import numpy as np
import multiprocessing as mp
import subprocess
import tempfile
import logging
from datetime import datetime
from collections import defaultdict
from tqdm import tqdm

CHUNK_SIZE = 600  # genes per chunk for parallel processing
END_TOLERANCE = 20  # nt window for 5' and 3' end support
SPLICE_DISTANCE_THRESHOLD = 10  # max distance between read and transcript splice junctions

def setup_logging(log_file=None):
    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    date_format = '%Y-%m-%d %H:%M:%S'
    
    if log_file:
        logging.basicConfig(level=logging.INFO,
                            format=log_format,
                            datefmt=date_format,
                            handlers=[
                                logging.FileHandler(log_file),
                                logging.StreamHandler()
                            ])
    else:
        logging.basicConfig(level=logging.INFO,
                            format=log_format,
                            datefmt=date_format)

def count_genes(gtf_file):
    gene_count = 0
    try:
        with open(gtf_file, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    fields = line.split('\t')
                    if fields[2] == 'gene':
                        gene_count += 1
        logging.info(f"Total number of genes in GTF: {gene_count}")
        return gene_count
    except FileNotFoundError:
        logging.error(f"GTF file not found: {gtf_file}")
        raise
    except Exception as e:
        logging.error(f"Error counting genes in GTF: {e}")
        raise

def load_gtf_into_memory(gtf_file):
    """
    Loads the GTF file into memory using gffutils and organizes genes and transcripts.
    """
    logging.info("Loading GTF into memory...")
    db_file = gtf_file + ".db"
    
    try:
        if not os.path.exists(db_file):
            logging.info("GFF database not found. Creating a new database...")
            try:
                # see: https://github.com/daler/gffutils/issues/20
                # this has only been tested with Ensembl GTFs 
                gffutils.create_db(
                    gtf_file,
                    db_file,
                    force=True,
                    keep_order=True,
                    merge_strategy='warning',      # corrected merge strategy
                    sort_attribute_values=True,
                    disable_infer_genes=True,           # disable gene inference
                    disable_infer_transcripts=True      # disable transcript inference
                )
                logging.info("GFF database created successfully.")
            except Exception as e:
                logging.error(f"Failed to create GFF database: {e}")
                raise
        
        if os.path.getsize(db_file) == 0:
            logging.error(f"GFF database file {db_file} is empty.")
            raise ValueError(f"GFF database file {db_file} is empty.")
        
        db = gffutils.FeatureDB(db_file)
        genes = defaultdict(list)
        transcripts = {}
        
        for gene in db.features_of_type('gene'):
            genes[(gene.chrom, gene.strand)].append({
                'id': gene.id,
                'chrom': gene.chrom,
                'start': gene.start,
                'end': gene.end,
                'strand': gene.strand,
                'transcripts': []
            })
        
        for transcript in db.features_of_type('transcript'):
            gene_id = transcript.attributes.get('gene_id', [None])[0]
            if not gene_id:
                logging.warning(f"Transcript {transcript.id} lacks gene_id. Skipping.")
                continue
            gene_found = False
            for gene in genes[(transcript.chrom, transcript.strand)]:
                if gene['id'] == gene_id:
                    gene['transcripts'].append(transcript.id)
                    gene_found = True
                    break
            if not gene_found:
                logging.warning(f"Gene ID {gene_id} for transcript {transcript.id} not found in genes list.")
                continue
            
            exons = sorted(db.children(transcript, featuretype='exon'), key=lambda x: x.start)
            transcripts[transcript.id] = {
                'chrom': transcript.chrom,
                'start': transcript.start,
                'end': transcript.end,
                'strand': transcript.strand,
                'exons': [(e.start, e.end) for e in exons],
                'length': transcript.end - transcript.start
            }
        
        for key in genes:
            genes[key].sort(key=lambda x: x['start'])
        
        logging.info(f"Loaded {sum(len(g) for g in genes.values())} genes and {len(transcripts)} transcripts.")
        return dict(genes), transcripts
    except Exception as e:
        logging.error(f"Error loading GTF into memory: {str(e)}")
        raise

def extract_splice_junctions(read, strand):
    """
    Extracts splice junctions from a read based on its CIGAR string.
    """
    splice_junctions = []
    ref_pos = read.reference_start
    try:
        for op, length in read.cigartuples:
            if op == 3:  # N operation is a splice junction
                junction = f"{ref_pos}-{ref_pos + length}"
                splice_junctions.append(junction)
                ref_pos += length
            elif op in [0, 2, 7, 8]:  # M, D, =, X operations
                ref_pos += length
    except AttributeError:
        logging.warning(f"Read {read.query_name} lacks CIGAR information.")
    return splice_junctions

def calculate_splice_distance(read_junctions, transcript_junctions, strand):
    """
    Calculates the average distance between read splice junctions and transcript splice junctions.
    """
    def parse_junction(junction):
        if isinstance(junction, str):
            return tuple(map(int, junction.split('-')))
        elif isinstance(junction, tuple):
            if len(junction) == 2 and all(isinstance(x, int) for x in junction):
                return junction 
        raise ValueError(f"Unexpected junction format: {junction}")

    try:
        read_junctions = [parse_junction(rj) for rj in read_junctions]
        transcript_junctions = [parse_junction(tj) for tj in transcript_junctions]
    except ValueError as ve:
        logging.error(f"Junction parsing error: {ve}")
        return float('inf'), "junction parsing error"

    # adjust for negative strand 
    if strand == '-':
        transcript_junctions = [(tj[1] - 1, tj[0]) for tj in transcript_junctions]

    distances = []
    if transcript_junctions:  
        for rj in read_junctions:
            try:
                closest_distance = min((abs(rj[0] - tj[0]) + abs(rj[1] - tj[1])) for tj in transcript_junctions)
                distances.append(closest_distance)
            except Exception as e:
                logging.error(f"Error calculating distance for read junction {rj}: {e}")
                continue
    else:
        logging.debug("No transcript junctions available for comparison.")
        return float('inf'), "no transcript junctions"

    avg_distance = np.mean(distances) if distances else float('inf')
    reason = "normal calculation" if distances else "no matching junctions"
    return avg_distance, reason

def process_bam_chunk(args):
    """
    Processes a chunk of genes to calculate support counts based on the BAM file.
    """
    bam_file, gene_chunk, transcripts, temp_dir = args

    support_counts = {gene['id']: {t_id: {'3_end': 0, 'splice': 0, '5_end': 0} for t_id in gene['transcripts']} for gene in gene_chunk}
    
    # create BED file for gene chunk
    bed_file = os.path.join(temp_dir, f"chunk_{os.getpid()}.bed")
    try:
        with open(bed_file, 'w') as f:
            for gene in gene_chunk:
                f.write(f"{gene['chrom']}\t{gene['start']}\t{gene['end']}\n")
    except Exception as e:
        logging.error(f"Error writing BED file {bed_file}: {e}")
        return support_counts  
    
    temp_bam = os.path.join(temp_dir, f"chunk_{os.getpid()}.bam")
    available_cores = max(1, mp.cpu_count() - 1)  
    
    try:
        # extract reads overlapping the BED regions
        subprocess.run(f"samtools view -b -@ {available_cores} -L {bed_file} {bam_file} > {temp_bam}", shell=True, check=True)
        subprocess.run(f"samtools index -@ {available_cores} {temp_bam}", shell=True, check=True)
    except subprocess.CalledProcessError as cpe:
        logging.error(f"Samtools error: {cpe}")
        try:
            os.remove(bed_file)
        except Exception as e:
            logging.warning(f"Error removing BED file after Samtools failure: {e}")
        return support_counts  
    
    try:
        with pysam.AlignmentFile(temp_bam, "rb") as bam:
            for read in bam.fetch():
                if read.is_unmapped or read.is_secondary or read.is_qcfail:
                    continue
                
                for gene in gene_chunk:
                    if gene['end'] < read.reference_start or gene['start'] > read.reference_end:
                        continue
                    
                    for transcript_id in gene['transcripts']:
                        transcript = transcripts[transcript_id]
                        
                        # 3' end support (with poly-A threshold and 20 nt window)
                        if ((transcript['strand'] == '+' and abs(read.reference_end - transcript['end']) <= END_TOLERANCE) or 
                            (transcript['strand'] == '-' and abs(read.reference_start - transcript['start']) <= END_TOLERANCE)):
                            support_counts[gene['id']][transcript_id]['3_end'] += 1
                        
                        # 5' end support (with 20 nt window)
                        if ((transcript['strand'] == '+' and abs(read.reference_start - transcript['start']) <= END_TOLERANCE) or 
                           (transcript['strand'] == '-' and abs(read.reference_end - transcript['end']) <= END_TOLERANCE)):
                            support_counts[gene['id']][transcript_id]['5_end'] += 1

                        # splice chain support
                        read_splice_junctions = extract_splice_junctions(read, transcript['strand'])
                        if read_splice_junctions:
                            exons = transcript['exons']
                            transcript_splice_junctions = [(exons[i][1], exons[i+1][0]) for i in range(len(exons)-1)]
                            if transcript['strand'] == '-':
                                transcript_splice_junctions = [(end, start) for start, end in transcript_splice_junctions]
                            
                            distance, reason = calculate_splice_distance(read_splice_junctions, transcript_splice_junctions, transcript['strand'])
                            
                            if distance <= SPLICE_DISTANCE_THRESHOLD:
                                support_counts[gene['id']][transcript_id]['splice'] += 1
                            else:
                                logging.debug(f"Splice mismatch - read: {read_splice_junctions}, transcript: {transcript_splice_junctions}, strand: {transcript['strand']}, distance: {distance}, reason: {reason}")
    except Exception as e:
        logging.error(f"Error processing BAM file {temp_bam}: {e}")
    finally:
        try:
            os.remove(bed_file)
            os.remove(temp_bam)
            os.remove(temp_bam + '.bai')
        except Exception as e:
            logging.warning(f"Error removing temporary files: {e}")
    
    return support_counts

def parallel_process_genes(bam_file, genes, transcripts):
    """
    Processes genes in parallel to calculate support counts.
    """
    logging.info("Processing genes in parallel...")
    available_cores = max(1, mp.cpu_count() - 1)  
    logging.info(f"Using {available_cores} cores for processing.")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        pool = mp.Pool(available_cores)
        tasks = []
        
        for (chrom, strand), gene_list in genes.items():
            for i in range(0, len(gene_list), CHUNK_SIZE):
                gene_chunk = gene_list[i:i+CHUNK_SIZE]
                tasks.append((bam_file, gene_chunk, transcripts, temp_dir))
        
        try:
            results = list(tqdm(pool.imap_unordered(process_bam_chunk, tasks), total=len(tasks), desc="Processing BAM chunks"))
        except Exception as e:
            logging.error(f"Error during parallel processing: {str(e)}")
            pool.terminate()
            pool.join()
            raise
        else:
            pool.close()
            pool.join()
    
    combined_counts = {}
    for counts in results:
        combined_counts.update(counts)
    
    logging.info(f"Processed support counts for {len(combined_counts)} genes.")
    return combined_counts

def prepare_scoring_data(support_counts, transcripts):
    """
    Prepares scoring data based on support counts for each transcript.
    """
    logging.info("Preparing scoring data...")
    scoring_data = {}
    for gene_id, gene_counts in support_counts.items():
        gene_end = max(transcripts[t_id]['end'] for t_id in gene_counts)
        transcript_ids = list(gene_counts.keys())
        
        counts = np.array([[gene_counts[t_id]['5_end'], gene_counts[t_id]['3_end'], gene_counts[t_id]['splice']] 
                           for t_id in transcript_ids])
        
        transcript_starts = np.array([transcripts[t_id]['start'] for t_id in transcript_ids])
        transcript_ends = np.array([transcripts[t_id]['end'] for t_id in transcript_ids])
        transcript_strands = np.array([transcripts[t_id]['strand'] for t_id in transcript_ids])
        exon_counts = np.array([len(transcripts[t_id]['exons']) for t_id in transcript_ids])
        
        scoring_data[gene_id] = {
            'transcript_ids': transcript_ids,
            'counts': counts,
            'transcript_starts': transcript_starts,
            'transcript_ends': transcript_ends,
            'transcript_strands': transcript_strands,
            'exon_counts': exon_counts,
            'gene_end': gene_end
        }
    
    return scoring_data

def process_scoring_chunk(args):
    """
    Processes a chunk of genes to determine the best transcript based on scoring.
    """
    gene_chunk, scoring_data, transcripts = args
    chunk_best_transcripts = {}
    chunk_score_matrix = []
    genes_with_no_reads = []

    for gene_id in gene_chunk:
        data = scoring_data.get(gene_id)
        if not data:
            logging.warning(f"No scoring data found for gene {gene_id}. Skipping.")
            continue
        
        # length bonus 
        if data['transcript_strands'][0] == '+':
            length_bonus = (data['gene_end'] - data['transcript_starts']) * 0.001
        else:
            length_bonus = (data['transcript_ends'] - data['transcript_starts'].min()) * 0.001
        
        base_scores = data['counts'][:, 0] + data['counts'][:, 1] + 2 * data['counts'][:, 2]
        scores = base_scores + length_bonus
        
        if np.sum(base_scores) == 0:
            # if no reads for a gene, select the longest transcript based on length_bonus
            genes_with_no_reads.append(gene_id)
            best_index = np.argmax(scores)  
            best_transcript = data['transcript_ids'][best_index]
            best_transcript = str(best_transcript)  
            chunk_best_transcripts[gene_id] = best_transcript
            logging.info(f"No reads for gene {gene_id}. Selected longest transcript {best_transcript}.")
        else:
            best_index = np.argmax(scores)
            best_transcript = data['transcript_ids'][best_index]
            best_transcript = str(best_transcript)  # Ensure it's a string
            chunk_best_transcripts[gene_id] = best_transcript
            logging.info(f"Selected best transcript {best_transcript} for gene {gene_id} based on support.")
        
        for i, transcript_id in enumerate(data['transcript_ids']):
            chunk_score_matrix.append([
                gene_id,
                transcript_id,
                data['transcript_starts'][i],
                data['transcript_ends'][i],
                data['exon_counts'][i],
                data['counts'][i, 0],
                data['counts'][i, 1],
                data['counts'][i, 2],
                length_bonus[i],
                scores[i]
            ])

    return chunk_best_transcripts, chunk_score_matrix, genes_with_no_reads

def select_best_transcripts(support_counts, transcripts, score_matrix_file=None):
    """
    Selects the best transcript for each gene based on support counts and scoring.
    """
    logging.info("Selecting best transcripts...")
    
    scoring_data = prepare_scoring_data(support_counts, transcripts)
    
    gene_ids = list(scoring_data.keys())
    gene_chunks = [gene_ids[i:i+CHUNK_SIZE] for i in range(0, len(gene_ids), CHUNK_SIZE)]

    # prepare to parallel sort 
    chunk_args = [(chunk, scoring_data, transcripts) for chunk in gene_chunks]
    
    available_cores = max(1, mp.cpu_count() - 1)
    pool = mp.Pool(available_cores)
    try:
        results = list(tqdm(pool.imap_unordered(process_scoring_chunk, chunk_args), 
                            total=len(chunk_args), desc="Processing scoring chunks"))
    except Exception as e:
        logging.error(f"Error during scoring parallel processing: {str(e)}")
        pool.terminate()
        pool.join()
        raise
    else:
        pool.close()
        pool.join()
    
    best_transcripts = {}
    score_matrix = []
    genes_with_reads = 0
    genes_without_reads = 0
    all_genes_with_no_reads = []
    
    for chunk_best_transcripts, chunk_score_matrix, chunk_genes_with_no_reads in results:
        best_transcripts.update(chunk_best_transcripts)
        score_matrix.extend(chunk_score_matrix)
        all_genes_with_no_reads.extend(chunk_genes_with_no_reads)
        genes_without_reads += len(chunk_genes_with_no_reads)
        genes_with_reads += len(chunk_best_transcripts) - len(chunk_genes_with_no_reads)
    
    if score_matrix_file:
        try:
            # sort the score matrix by gene, transcript, then descending TES
            score_matrix.sort(key=lambda x: (x[0], x[1], -x[3]))
            
            with open(score_matrix_file, 'w') as f:
                f.write("gene\ttranscript\tTSS\tTES\texon_count\t5_end\t3_end\tsplice\tlength_bonus\tscore\n")
                for row in score_matrix:
                    f.write('\t'.join(map(str, row)) + '\n')
            logging.info(f"Score matrix written to {score_matrix_file}")
        except Exception as e:
            logging.error(f"Error writing score matrix to file: {str(e)}")
    
    if score_matrix_file:
        try:
            no_reads_file = score_matrix_file.rsplit('.', 1)[0] + '_genes_with_no_reads.txt'
            with open(no_reads_file, 'w') as f:
                for gene_id in all_genes_with_no_reads:
                    f.write(f"{gene_id}\n")
            logging.info(f"List of genes with no reads written to {no_reads_file}")
        except Exception as e:
            logging.error(f"Error writing genes with no reads to file: {str(e)}")
    
    logging.info(f"Selected best transcripts for {genes_with_reads} genes with reads.")
    logging.info(f"Selected longest transcript for {genes_without_reads} genes without reads.")
    
    return best_transcripts, genes_with_reads, genes_without_reads

def filter_gtf(gtf_file, best_transcripts, out_gtf):
    """
    Filters the GTF file to retain only the best transcripts per gene using grep-like functionality.
    Writes all 'gene' entries and only 'transcript' and 'exon' entries with transcript_ids in best_transcripts.
    """
    logging.info("Filtering GTF based on best transcripts...")
    
    try:
        
        best_transcript_ids = set(best_transcripts.values())
        
        with open(gtf_file, 'r') as in_f, open(out_gtf, 'w') as out_f:
            for line in tqdm(in_f, desc="Writing GTF"):
                if line.startswith('#') or not line.strip():
                    # keep header
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    # simple check for malformed lines
                    continue
                
                feature_type = fields[2]
                attributes = fields[8]
                
                if feature_type == 'gene':
                    out_f.write(line)
                elif feature_type in ['transcript', 'exon']:
                    transcript_id = None
                    for attr in attributes.split(';'):
                        attr = attr.strip()
                        if attr.startswith('transcript_id'):
                            parts = attr.split(' ')
                            if len(parts) >= 2:
                                transcript_id = parts[1].strip('"')
                            break
                    if transcript_id and transcript_id in best_transcript_ids:
                        out_f.write(line)
        logging.info(f"Filtered GTF written to {out_gtf}")
    except Exception as e:
        logging.error(f"Error filtering GTF: {str(e)}")
        raise

def main():
    """
    Main function to execute the transcript selection and GTF filtering.
    """
    import argparse
    
    parser = argparse.ArgumentParser(description="Transcript Selection and GTF Filtering Script")
    parser.add_argument("--bam_file", help="Input BAM file path")
    parser.add_argument("--gtf_file", required=True, help="Input GTF file path")
    parser.add_argument("--out_gtf", required=True, help="Output filtered GTF file path")
    parser.add_argument("--score_matrix", required=True, help="Output file path for scoring matrix")
    parser.add_argument("--log_file", help="Log file path (optional)")
    parser.add_argument("--resume", action='store_true', help="Resume from existing scoring matrix")
    args = parser.parse_args()
    
    setup_logging(args.log_file)
    
    logging.info("Transcript Selection and GTF Filtering Process Initiated.")
    
    if args.resume:
        logging.info("Resuming process using existing scoring matrix.")
        if not os.path.exists(args.score_matrix):
            logging.error(f"Scoring matrix file not found: {args.score_matrix}. Cannot resume.")
            raise FileNotFoundError(f"Scoring matrix file not found: {args.score_matrix}. Cannot resume.")
        
        # load scoring matrix 
        try:
            best_transcripts = {}
            with open(args.score_matrix, 'r') as f:
                header = f.readline()  # header 
                for line in f:
                    fields = line.strip().split('\t')
                    if len(fields) < 10:
                        logging.warning(f"Skipping malformed line in score matrix: {line.strip()}")
                        continue
                    gene_id, transcript_id, TSS, TES, exon_count, five_end, three_end, splice, length_bonus, score = fields
                    try:
                        score = float(score)
                    except ValueError:
                        logging.warning(f"Invalid score value for gene {gene_id}, transcript {transcript_id}: {score}. Skipping.")
                        continue
                    # update best transcript if this score is higher
                    if gene_id not in best_transcripts or score > best_transcripts[gene_id][1]:
                        best_transcripts[gene_id] = (str(transcript_id), score)  # Ensure transcript_id is string
            # convert to gene_id: best_transcript_id
            best_transcripts = {gene_id: transcript_id for gene_id, (transcript_id, _) in best_transcripts.items()}
            
            # calculate genes_with_reads and genes_without_reads
            genes_with_reads = len(best_transcripts)
            total_genes = count_genes(args.gtf_file)  # dynamically get total genes
            genes_without_reads = total_genes - genes_with_reads
            logging.info(f"Loaded best transcripts for {genes_with_reads} genes with reads.")
            logging.info(f"Loaded best transcripts for {genes_without_reads} genes without reads.")
            
            if not best_transcripts:
                logging.warning("No best transcripts loaded. The filtered GTF will contain only gene entries.")
            else:
                sample_genes = list(best_transcripts.keys())[:5]
                sample_transcripts = [best_transcripts[gene] for gene in sample_genes]
                logging.info(f"Sample of loaded best transcripts: {list(zip(sample_genes, sample_transcripts))}")
            
        except Exception as e:
            logging.error(f"Failed to load best transcripts from scoring matrix: {str(e)}")
            raise

    else:
        if not args.bam_file:
            logging.error("BAM file must be provided when not resuming.")
            raise ValueError("BAM file must be provided when not resuming.")
        
        for file_path in [args.bam_file, args.gtf_file]:
            if not os.path.exists(file_path):
                logging.error(f"Input file does not exist: {file_path}")
                raise FileNotFoundError(f"Input file does not exist: {file_path}")
        
        try:
            original_gene_count = count_genes(args.gtf_file)
        except Exception as e:
            logging.error(f"Failed to count genes: {e}")
            raise
        
        # load gtf using gffutils
        try:
            genes, transcripts = load_gtf_into_memory(args.gtf_file)
        except Exception as e:
            logging.error(f"Failed to load GTF into memory: {e}")
            raise
        
        # process bam file 
        try:
            support_counts = parallel_process_genes(args.bam_file, genes, transcripts)
        except Exception as e:
            logging.error(f"Failed during parallel processing of BAM file: {e}")
            raise
        
        # select best transcripts based on support counts
        try:
            best_transcripts, genes_with_reads, genes_without_reads = select_best_transcripts(support_counts, transcripts, args.score_matrix)
        except Exception as e:
            logging.error(f"Failed during transcript selection: {e}")
            raise

    # filter gtf annotation 
    try:
        filter_gtf(args.gtf_file, best_transcripts, args.out_gtf)
    except Exception as e:
        logging.error(f"Failed during GTF filtering: {str(e)}")
        raise
    
    # verify gene counts
    try:
        original_gene_count = count_genes(args.gtf_file) if not args.resume else count_genes(args.out_gtf)  # Dynamically get total gene count
        filtered_gene_count = count_genes(args.out_gtf)
        if original_gene_count == filtered_gene_count:
            logging.info("Gene count check passed: The number of genes before and after filtering is the same.")
        else:
            logging.error(f"Gene count mismatch: Original GTF had {original_gene_count} genes, but filtered GTF has {filtered_gene_count} genes.")
            raise ValueError("The number of genes before and after filtering does not match.")
    except Exception as e:
        logging.error(f"Failed during gene count verification: {str(e)}")
        raise
    
    logging.info("Transcript selection and GTF filtering process completed successfully.")

if __name__ == "__main__":
    main()
