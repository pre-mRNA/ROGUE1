import pandas as pd
import numpy as np
from collections import defaultdict
import logging
import os
import tempfile
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from parallel_bed_operations import split_bed_file

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

def process_read(intervals, window_size=6):
    intervals.sort(key=lambda x: x[1])  # sort by start 
    total_width = sum(x[2] - x[1] for x in intervals)
    if total_width < window_size:
        return []
    
    strand = intervals[0][5]
    remaining = window_size
    result = []
    
    if strand == '+':
        for interval in reversed(intervals):
            chrom, start, end = interval[0], interval[1], interval[2]
            width = end - start
            if width >= remaining:
                result.append((chrom, end - remaining, end) + interval[3:])
                break
            else:
                result.append(interval)
                remaining -= width
    else:
        for interval in intervals:
            chrom, start, end = interval[0], interval[1], interval[2]
            width = end - start
            if width >= remaining:
                result.append((chrom, start, start + remaining) + interval[3:])
                break
            else:
                result.append(interval)
                remaining -= width
    
    return result

def intersect_and_process_chunk(bed_chunk, feature_files, genome_file, output_dir, gene_id_map, min_overlap=4):
    results = {}

    # intersect reads against exons, introns and DOG regions
    for feature, feature_file in feature_files.items():

        intersect_file = os.path.join(output_dir, f'{os.path.basename(bed_chunk)}_intersect_3prime_{feature}.bed')
        cmd = f"bedtools intersect -a {bed_chunk} -b {feature_file} -wo -s -sorted -g {genome_file} > {intersect_file}"
        subprocess.run(cmd, shell=True, check=True)
        
        with open(intersect_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                read_id = fields[3] # read ID is third col
                feature_info = fields[-3] # feature info is in the third to last col 
                parts = feature_info.split('_')
                if len(parts) >= 3:
                    gene_id = '_'.join(parts[:-2])  
                    feature = '_'.join(parts[-2:])
                else:
                    gene_id = feature_info
                    feature = "Unknown"

                if read_id not in gene_id_map or gene_id_map[read_id] != gene_id:
                    continue

                overlap = int(fields[-1])
                if overlap < min_overlap:
                    continue
                
                if read_id not in results:
                    results[read_id] = {}
                if feature not in results[read_id]:
                    results[read_id][feature] = 0
                results[read_id][feature] += overlap

    return results

def get_end_feature(temp_bed_files, genome_file, output_dir, dog_bed_file, exon_bed_file, intron_bed_file, gene_id_map, min_overlap=4, num_chunks=104):
    logging.info(f"Calculating read end feature for {len(temp_bed_files)} bed files")
    
    all_reads = defaultdict(list)
    
    for bed_chunk in temp_bed_files:
        df = pd.read_csv(bed_chunk, sep='\t', header=None, 
                         names=['chrom', 'start', 'end', 'read_id', 'score', 'strand'])
        for _, row in df.iterrows():
            all_reads[row['read_id']].append(tuple(row))
    
    logging.info(f"Total reads for end feature processing: {len(all_reads)}")
    
    processed_reads = []
    for read_id, intervals in all_reads.items():
        processed = process_read(intervals, window_size=6)
        processed_reads.extend(processed)
    
    with tempfile.NamedTemporaryFile(mode='w+t', dir=output_dir, suffix='_processed_3prime_6nt.bed', delete=False) as temp_file:
        for interval in processed_reads:
            temp_file.write('\t'.join(map(str, interval)) + '\n')
        processed_file = temp_file.name
    
    logging.info(f"Read end positions written to {processed_file}")
    
    # use system sort
    logging.info("Sorting read end positions reads")
    sorted_file = os.path.join(output_dir, 'sorted_processed_3prime_6nt.bed')
    sort_command = f"sort -V -k1,1 -k2,2n -k3,3n --parallel={os.cpu_count()} -S 80% -T {output_dir} {processed_file} > {sorted_file}"

    try:
        subprocess.run(sort_command, shell=True, check=True)
        logging.info(f"Sorted reads written to {sorted_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error occurred while sorting: {e}")
        raise

    processed_file = sorted_file
    
    # split the 3' end positions into chunks for parallel processing
    temp_bed_chunks = split_bed_file(processed_file, output_dir, num_chunks)
    
    # refer to exon, intron and dog files to assign gene positions 
    feature_files = {
        "exon": exon_bed_file,
        "intron": intron_bed_file,
        "dog": dog_bed_file
    }
    
    # use multiprocessing to intersect the read end positions against exons, introns and DOG regions 
    with ThreadPoolExecutor(max_workers=num_chunks) as executor:
        future_to_chunk = {executor.submit(intersect_and_process_chunk, chunk, feature_files, genome_file, output_dir, gene_id_map, min_overlap): chunk for chunk in temp_bed_chunks}
        chunk_results = []
        for future in as_completed(future_to_chunk):
            chunk_results.append(future.result())
    
    # merge chunks
    merged_results = {}
    for chunk_result in chunk_results:
        for read_id, features in chunk_result.items():
            if read_id not in merged_results:
                merged_results[read_id] = features
            else:
                for feature, overlap in features.items():
                    if feature not in merged_results[read_id]:
                        merged_results[read_id][feature] = overlap
                    else:
                        merged_results[read_id][feature] += overlap

    # retain the best feature per read_id based on maximum overlap
    final_data = []
    for read_id, features in merged_results.items():
        if features:
            best_feature = max(features.items(), key=lambda x: x[1])
            final_data.append({'read_id': read_id, 'gene_id': gene_id_map[read_id], 'feature': best_feature[0]})

    df = pd.DataFrame(final_data)
    if not df.empty:
        df.set_index('read_id', inplace=True)
        logging.info(f"DataFrame created with {len(df)} rows, indexed by read_id.")
    else:
        logging.warning("No valid data processed into DataFrame; DataFrame is empty.")

    return df