import pysam
import gffutils
import numpy as np
import multiprocessing as mp
import subprocess
import os
import tempfile
from collections import defaultdict
from tqdm import tqdm

CHUNK_SIZE = 500  # Number of genes per chunk
AVAILABLE_CORES = 200

def load_gtf_into_memory(gtf_file):
    db = gffutils.FeatureDB(gtf_file + ".db")
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
    
    # Sort genes by start position
    for key in genes:
        genes[key].sort(key=lambda x: x['start'])
    
    return dict(genes), transcripts

def process_gene_chunk(args):
    bam_file, gene_chunk, transcripts, temp_dir = args
    
    # Create a BED file for the gene chunk
    bed_file = os.path.join(temp_dir, f"chunk_{os.getpid()}.bed")
    with open(bed_file, 'w') as f:
        for gene in gene_chunk:
            f.write(f"{transcripts[gene['transcripts'][0]]['chrom']}\t{gene['start']}\t{gene['end']}\n")
    
    # Use samtools to extract reads for the gene chunk
    temp_bam = os.path.join(temp_dir, f"chunk_{os.getpid()}.bam")
    subprocess.run(f"samtools view -b -L {bed_file} {bam_file} > {temp_bam}", shell=True)
    subprocess.run(f"samtools index {temp_bam}", shell=True)
    
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
                    
                    # Check for 3' end support
                    if (transcript['strand'] == '+' and read.reference_end == gene['end']) or \
                       (transcript['strand'] == '-' and read.reference_start == transcript['start']):
                        support_counts[gene['id']][transcript_id]['3_end'] += 1
                    
                    # Check for splice chain support
                    if read.has_tag('XS'):
                        read_splice_junctions = read.get_tag('XS').split(',')
                        transcript_splice_junctions = [f"{e[1]}-{next_e[0]}" for e, next_e in zip(transcript['exons'][:-1], transcript['exons'][1:])]
                        if all(rj in transcript_splice_junctions for rj in read_splice_junctions):
                            support_counts[gene['id']][transcript_id]['splice'] += 1
                    
                    # Check for 5' end support
                    if (transcript['strand'] == '+' and read.reference_start == transcript['start']) or \
                       (transcript['strand'] == '-' and read.reference_end == transcript['end']):
                        support_counts[gene['id']][transcript_id]['5_end'] += 1
    
    # Clean up temporary files
    os.remove(bed_file)
    os.remove(temp_bam)
    os.remove(temp_bam + '.bai')
    
    return support_counts

def parallel_process_genes(bam_file, genes, transcripts):
    with tempfile.TemporaryDirectory() as temp_dir:
        pool = mp.Pool(AVAILABLE_CORES)
        tasks = []
        
        for (chrom, strand), gene_list in genes.items():
            for i in range(0, len(gene_list), CHUNK_SIZE):
                gene_chunk = gene_list[i:i+CHUNK_SIZE]
                tasks.append((bam_file, gene_chunk, transcripts, temp_dir))
        
        results = list(tqdm(pool.imap_unordered(process_gene_chunk, tasks), total=len(tasks), desc="Processing gene chunks"))
        
        pool.close()
        pool.join()
    
    # Combine results
    combined_counts = {}
    for counts in results:
        combined_counts.update(counts)
    
    return combined_counts

def select_best_transcripts(support_counts, transcripts):
    best_transcripts = {}
    for gene_id, gene_counts in support_counts.items():
        best_score = -1
        best_transcript = None
        for transcript_id, counts in gene_counts.items():
            score = sum(counts.values())
            if score > best_score:
                best_score = score
                best_transcript = transcript_id
        
        if best_score == 0:
            # If no support, choose the longest transcript
            best_transcript = max(gene_counts.keys(), key=lambda t: transcripts[t]['length'])
        
        best_transcripts[gene_id] = best_transcript
    
    return best_transcripts

def filter_gtf(gtf_file, best_transcripts, out_gtf):
    db = gffutils.FeatureDB(gtf_file + ".db")
    with open(out_gtf, 'w') as out_f:
        for gene_id, transcript_id in tqdm(best_transcripts.items(), desc="Writing GTF"):
            gene = db[gene_id]
            transcript = db[transcript_id]
            out_f.write(str(gene) + '\n')
            out_f.write(str(transcript) + '\n')
            for exon in db.children(transcript, featuretype='exon'):
                out_f.write(str(exon) + '\n')

def main(bam_file, gtf_file, out_gtf):
    print("Loading GTF into memory...")
    genes, transcripts = load_gtf_into_memory(gtf_file)

    print("Processing genes in parallel...")
    support_counts = parallel_process_genes(bam_file, genes, transcripts)

    print("Selecting best transcripts...")
    best_transcripts = select_best_transcripts(support_counts, transcripts)

    print("Filtering GTF...")
    filter_gtf(gtf_file, best_transcripts, out_gtf)

    print("Done!")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Ultra-fast transcript selector using gene chunks")
    parser.add_argument("--bam_file", required=True, help="Input BAM file path")
    parser.add_argument("--gtf_file", required=True, help="Input GTF file path")
    parser.add_argument("--out_gtf", required=True, help="Output GTF file path")
    args = parser.parse_args()

    main(args.bam_file, args.gtf_file, args.out_gtf)