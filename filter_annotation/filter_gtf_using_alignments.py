import os 
import pysam
import gffutils
import numpy as np
import multiprocessing as mp
import subprocess
import os
import tempfile
import logging
from datetime import datetime
from collections import defaultdict
from tqdm import tqdm
from functools import partial

# this code picks the most supported transcript for each gene 
# if no transcript has support, it picks the canonical (longest) transcript 
# it returns all exons for the selected transcript 
# and trims the gene feature to the length of the transcript

CHUNK_SIZE = 600  # genes per chunk for parallel processing

def setup_logging():
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

def load_gtf_into_memory(gtf_file):
    logging.info("Loading GTF into memory...")
    db_file = gtf_file + ".db"
    
    if not os.path.exists(db_file):
        logging.info("GFF database not found. Creating new database...")
        gffutils.create_db(gtf_file, db_file, force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
        logging.info("GFF database created successfully.")
    
    db = gffutils.FeatureDB(db_file)
    genes = defaultdict(list)
    transcripts = {}
    
    for gene in db.features_of_type('gene'):
        genes[(gene.chrom, gene.strand)].append({
            'id': gene.id,
            'start': gene.start,
            'end': gene.end,
            'transcripts': []
        })
    
    for transcript in db.features_of_type('transcript'):
        gene_id = transcript.attributes['gene_id'][0]
        for gene in genes[(transcript.chrom, transcript.strand)]:
            if gene['id'] == gene_id:
                gene['transcripts'].append(transcript.id)
                break
        
        exons = sorted(db.children(transcript, featuretype='exon'), key=lambda x: x.start)
        transcripts[transcript.id] = {
            'chrom': transcript.chrom,
            'start': transcript.start,
            'end': transcript.end,
            'strand': transcript.strand,
            'exons': [(e.start, e.end) for e in exons],
            'length': transcript.end - transcript.start
        }
    
    # sort genes by start
    for key in genes:
        genes[key].sort(key=lambda x: x['start'])
    
    logging.info(f"Loaded {sum(len(g) for g in genes.values())} genes and {len(transcripts)} transcripts")
    return dict(genes), transcripts

def extract_splice_junctions(read, strand):
    splice_junctions = []
    ref_pos = read.reference_start
    for op, length in read.cigartuples:
        if op == 3:  # N operation is a splice junction
            junction = f"{ref_pos}-{ref_pos + length}"
            splice_junctions.append(junction)
            ref_pos += length
        elif op in [0, 2, 7, 8]:  # M, D, =, X operations
            ref_pos += length
    return splice_junctions

def calculate_splice_distance(read_junctions, transcript_junctions, strand):
    
    # logging.info(f"Calculating distances for read junctions {read_junctions} on strand {strand} with transcript junctions {transcript_junctions}")
    
    def parse_junction(junction):
        """Ensure the junction is parsed correctly, whether it's a string or already a tuple."""
        if isinstance(junction, str):
            return tuple(map(int, junction.split('-')))
        elif isinstance(junction, tuple):
            if len(junction) == 2 and all(isinstance(x, int) for x in junction):
                return junction 
        raise ValueError(f"Unexpected junction format: {junction}")

    read_junctions = [parse_junction(rj) for rj in read_junctions]
    transcript_junctions = [parse_junction(tj) for tj in transcript_junctions]

    # adjust junctions for negative strand:
    if strand == '-':
        transcript_junctions = [(tj[1] - 1, tj[0]) for tj in transcript_junctions]

    distances = []
    if transcript_junctions:  
        for rj in read_junctions:
            closest_distance = min((abs(rj[0] - tj[0]) + abs(rj[1] - tj[1])) for tj in transcript_junctions)
            distances.append(closest_distance)
    else:
        # logging.info(f"No transcript junctions available for comparison on strand {strand} for read junctions {read_junctions}.")
        return float('inf'), "No transcript junctions"

    avg_distance = np.mean(distances) if distances else float('inf')
    reason = "Normal calculation" if distances else "No matching junctions"
    return avg_distance, reason

def process_bam_chunk(args):
    bam_file, gene_chunk, transcripts, temp_dir, poly_a_threshold = args
    
    # create a BED file representing the region of the gene chunk
    bed_file = os.path.join(temp_dir, f"chunk_{os.getpid()}.bed")
    with open(bed_file, 'w') as f:
        for gene in gene_chunk:
            f.write(f"{transcripts[gene['transcripts'][0]]['chrom']}\t{gene['start']}\t{gene['end']}\n")
    
    temp_bam = os.path.join(temp_dir, f"chunk_{os.getpid()}.bam")
    available_cores = mp.cpu_count()

    subprocess.run(f"samtools view -b -@ {available_cores} -L {bed_file} {bam_file} > {temp_bam}", shell=True)
    subprocess.run(f"samtools index -@ {available_cores} {temp_bam}", shell=True)
    
    END_TOLERANCE = 20  
    SPLICE_DISTANCE_THRESHOLD = 10  # max distance between read and transcript splice junctions

    splice_chains_observed = []

    support_counts = {gene['id']: {t_id: {'3_end': 0, 'splice': 0, '5_end': 0} for t_id in gene['transcripts']} for gene in gene_chunk}
    
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
                        if poly_a_threshold is None or (read.has_tag('pt') and read.get_tag('pt') >= poly_a_threshold):
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
                        
                        if len(splice_chains_observed) < 10:
                            splice_chains_observed.append((read_splice_junctions, transcript_splice_junctions, transcript['strand']))
                        
                        distance, reason = calculate_splice_distance(read_splice_junctions, transcript_splice_junctions, transcript['strand'])
                        
                        if distance <= SPLICE_DISTANCE_THRESHOLD:
                            support_counts[gene['id']][transcript_id]['splice'] += 1
                        else:
                            logging.debug(f"Splice mismatch - Read: {read_splice_junctions}, Transcript: {transcript_splice_junctions}, Strand: {transcript['strand']}, Distance: {distance}, Reason: {reason}")

#    # print the first 10 observed splice chains and their distances
#     for i, (read_junctions, transcript_junctions, strand) in enumerate(splice_chains_observed[:10], 1):
#         distance, reason = calculate_splice_distance(read_junctions, transcript_junctions, strand)
#         logging.info(f"Observed splice chain {i}:")
#         logging.info(f"  Read junctions: {read_junctions}")
#         logging.info(f"  Closest transcript junctions: {transcript_junctions}")
#         logging.info(f"  Strand: {strand}")
#         logging.info(f"  Average distance: {distance}")
#         logging.info(f"  Reason: {reason}")

    os.remove(bed_file)
    os.remove(temp_bam)
    os.remove(temp_bam + '.bai')
    
    return support_counts

def parallel_process_genes(bam_file, genes, transcripts, poly_a_threshold):
    logging.info("Processing genes in parallel...")
    available_cores = mp.cpu_count()
    logging.info(f"Using all {available_cores} available cores")

    with tempfile.TemporaryDirectory() as temp_dir:
        pool = mp.Pool(available_cores)
        tasks = []
        
        for (chrom, strand), gene_list in genes.items():
            for i in range(0, len(gene_list), CHUNK_SIZE):
                gene_chunk = gene_list[i:i+CHUNK_SIZE]
                tasks.append((bam_file, gene_chunk, transcripts, temp_dir, poly_a_threshold))
        
        results = list(tqdm(pool.imap_unordered(process_bam_chunk, tasks), total=len(tasks), desc="Processing BAM chunks"))
        
        pool.close()
        pool.join()
    
    # combine results
    combined_counts = {}
    for counts in results:
        combined_counts.update(counts)
    
    logging.info(f"Processed {len(combined_counts)} genes")
    return combined_counts

def prepare_scoring_data(support_counts, transcripts):
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
    gene_chunk, scoring_data = args
    chunk_best_transcripts = {}
    chunk_score_matrix = []

    for gene_id in gene_chunk:
        data = scoring_data[gene_id]
        print(data)
        
        # calculate length bonus 
        if data['transcript_strands'][0] == '+':
            length_bonus = (data['gene_end'] - data['transcript_starts']) * 0.001
        else:
            length_bonus = (data['transcript_ends'] - data['transcript_starts'].min()) * 0.001
        
        base_scores = data['counts'][:, 0] + data['counts'][:, 1] + 2 * data['counts'][:, 2]
        scores = base_scores + length_bonus
        
        best_index = np.argmax(scores)
        best_transcript = data['transcript_ids'][best_index]
        
        if scores[best_index] == 0:

            # if no support, choose the transcript with the most upstream end (considering strand)
            if data['transcript_strands'][0] == '+':
                best_transcript = data['transcript_ids'][np.argmin(data['transcript_starts'])]
            else:
                best_transcript = data['transcript_ids'][np.argmax(data['transcript_ends'])]
        
        chunk_best_transcripts[gene_id] = best_transcript
        
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

    return chunk_best_transcripts, chunk_score_matrix

def select_best_transcripts(support_counts, transcripts, score_matrix_file=None):
    logging.info("Selecting best transcripts...")
    
    scoring_data = prepare_scoring_data(support_counts, transcripts)
    
    gene_chunks = [list(scoring_data.keys())[i:i+CHUNK_SIZE] for i in range(0, len(scoring_data), CHUNK_SIZE)]

    chunk_args = [(chunk, scoring_data) for chunk in gene_chunks]

    with mp.Pool(processes=mp.cpu_count()) as pool:
        results = list(tqdm(pool.imap(process_scoring_chunk, chunk_args), 
                            total=len(chunk_args), desc="Processing scoring chunks"))

    best_transcripts = {}
    score_matrix = []
    for chunk_best_transcripts, chunk_score_matrix in results:
        best_transcripts.update(chunk_best_transcripts)
        score_matrix.extend(chunk_score_matrix)

    if score_matrix_file:
        
        score_matrix.sort(key=lambda x: (x[0], x[1], -x[3]))  # sort by gene, transcript, then -TES
        
        with open(score_matrix_file, 'w') as f:
            f.write("gene\ttranscript\tTSS\tTES\texon_count\t5_end\t3_end\tsplice\tlength_bonus\tscore\n")
            for row in score_matrix:
                f.write('\t'.join(map(str, row)) + '\n')
        
        logging.info(f"Score matrix written to {score_matrix_file}")
    
    logging.info(f"Selected {len(best_transcripts)} best transcripts")
    return best_transcripts

def filter_gtf(gtf_file, best_transcripts, out_gtf):
    logging.info("Filtering GTF...")
    db = gffutils.FeatureDB(gtf_file + ".db")
    with open(out_gtf, 'w') as out_f:
        for gene_id, transcript_id in tqdm(best_transcripts.items(), desc="Writing GTF"):
            gene = db[gene_id]
            transcript = db[transcript_id]
            
            # trim gene feature length 
            gene_start = transcript.start
            gene_end = transcript.end
            
            # construct attribute string 
            gene_attrs = '; '.join([f"{k} \"{v[0]}\"" for k, v in gene.attributes.items()])
            gene_line = f"{gene.chrom}\t{gene.source}\tgene\t{gene_start}\t{gene_end}\t.\t{gene.strand}\t.\t{gene_attrs}"
            
            out_f.write(gene_line + '\n')
            out_f.write(str(transcript) + '\n')
            for exon in db.children(transcript, featuretype='exon'):
                out_f.write(str(exon) + '\n')
    
    logging.info(f"Filtered GTF written to {out_gtf}")

def main(bam_file, gtf_file, out_gtf, poly_a_threshold, score_matrix_file):
    setup_logging()
    logging.info("R1 annotation filtering initialised")
    logging.info("Starting transcript selection process")
    
    genes, transcripts = load_gtf_into_memory(gtf_file)
    support_counts = parallel_process_genes(bam_file, genes, transcripts, poly_a_threshold)
    best_transcripts = select_best_transcripts(support_counts, transcripts, score_matrix_file)
    filter_gtf(gtf_file, best_transcripts, out_gtf)
    
    logging.info("Transcript selection process completed")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="ROGUE1 adaptive annotation filtering")
    parser.add_argument("--bam_file", required=True, help="Input BAM file path")
    parser.add_argument("--gtf_file", required=True, help="Input GTF file path")
    parser.add_argument("--out_gtf", required=True, help="Output GTF file path")
    parser.add_argument("--polyA", type=int, help="Minimum poly-A tail length for 3' end support")
    parser.add_argument("--score_matrix", help="Output file path for scoring matrix")
    args = parser.parse_args()

    main(args.bam_file, args.gtf_file, args.out_gtf, args.polyA, args.score_matrix)