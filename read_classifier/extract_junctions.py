import pandas as pd 
import pysam
import multiprocessing
from functools import partial
import logging 

# if we're using a second-pass index
# instead of de novo calculating introns
# we can just use all the introns in the refined annotation 
# this function converts the updated_intron.bed file to a format that can be used by the rest of the pipeline
def process_intron_to_junctions(intron_bed):
    df = pd.read_csv(intron_bed, sep='\t', header=None, 
                     names=['Chromosome', 'Start', 'End', 'Gene_ID', 'Intron_ID', 'Strand'])
    df['Counts'] = -1  # default count for annotated junction 
    df = df[['Chromosome', 'Start', 'End', 'Gene_ID', 'Counts', 'Strand']]
    return df

def extract_region_splice_junctions(bam_file, region, gene_id_map):
    """
    Extract splice junctions for a specific region from a BAM file.
    """

    bam = pysam.AlignmentFile(bam_file, "rb")
    results = []
    for read in bam.fetch(region=region):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        gene_id = gene_id_map.get(read.query_name, "unknown")
        current_pos = read.reference_start
        # print(f"{read.cigartuples}")
        for operation, length in read.cigartuples:
            if operation in {0, 7, 8}:  # match
                current_pos += length
            elif operation == 3:  # splice
                # print(f"match, N, {length}")
                start = current_pos
                end = current_pos + length
                strand = '-' if read.is_reverse else '+'
                # print(f"{read.reference_name}, {start}, {end}, {strand}")
                results.append([read.reference_name, start, end, strand, gene_id])
                current_pos += length
            elif operation == 2:  # deletion or softclip
                current_pos += length
            elif operation == 1:  # insertion 
                continue
    bam.close()
    return results

def parallel_extract_splice_junctions(bam_file, num_cpus, gene_id_map):
    """
    Extract splice junctions from a BAM file in parallel using multiple CPUs.
    """
    logging.info("Extracting splice junctions")
    bam = pysam.AlignmentFile(bam_file, "rb")

    contigs = bam.references
    bam.close()

    with multiprocessing.Pool(num_cpus) as pool:
        func = partial(extract_region_splice_junctions, bam_file, gene_id_map=gene_id_map)
        results = pool.map(func, contigs)

    # flatten into df 
    flat_results = [item for sublist in results for item in sublist]
    df = pd.DataFrame(flat_results, columns=['Chromosome', 'Start', 'End', 'Strand', 'Gene_ID'])
    return df

def summarize_junctions(df, min_reads):
    """
    Summarize and filter splice junctions based on a minimum number of supporting reads.

    Parameters:
        df (pd.DataFrame): DataFrame containing splice junctions data with columns 'Chromosome', 'Start', 'End', 'Strand'.
        min_reads (int): Minimum number of reads supporting a junction to be retained.

    Returns:
        pd.DataFrame: DataFrame of unique junctions with counts of reads supporting each junction.
    """

    logging.info("Collapsing splice junctions")
    print(f"{len(df)} junctions")
    print(f"{df.drop_duplicates().shape[0]} unique junctions")

    # group junctions 
    grouped = df.groupby(['Chromosome', 'Start', 'End', 'Strand', 'Gene_ID']).size().reset_index(name='Counts')
    
    # filter junctions based on individual support 
    filtered = grouped[grouped['Counts'] >= min_reads]
    
    print(f"{filtered.drop_duplicates().shape[0]} unique junctions after filtering")
    
    return filtered


def calculate_overlaps(group, overlap_threshold):
    """Calculate overlapping introns within a DataFrame group."""
    overlaps = []
    rows = list(group.itertuples())  
    for i, row_i in enumerate(rows):
        for j in range(i + 1, len(rows)):
            row_j = rows[j]
            overlap_length = min(row_i.End, row_j.End) - max(row_i.Start, row_j.Start)
            if overlap_length > 0:
                ratio_i = overlap_length / row_i.Length
                ratio_j = overlap_length / row_j.Length
                if ratio_i >= overlap_threshold and ratio_j >= overlap_threshold:
                    overlaps.append((row_i.Index, row_j.Index))
    return overlaps

def process_group(group, overlap_threshold):
    """Process each group for parallel execution."""
    group['Length'] = group['End'] - group['Start']
    total_counts = group['Counts'].sum()
    group['Relative_Frequency'] = group['Counts'] / total_counts if total_counts > 0 else 0
    overlaps = calculate_overlaps(group, overlap_threshold)
    for i, j in overlaps:
        group.at[i, 'Overlapping_Index'] = j
        group.at[j, 'Overlapping_Index'] = i
    return group

def filter_junctions(df, overlap_threshold=0.60, frequency_threshold=0.002, num_cores=None):
    """
    Identify introns that reciprocally overlap more than a specified threshold and filter them by relative frequency using parallel processing.
    """
    logging.info("Filtering splice junctions")
    initial_count = len(df)
    print(f"Initial number of junctions: {initial_count}")

    if num_cores is None:
        num_cores = multiprocessing.cpu_count()

    # create groups to reduce overlap search space 
    grouped = df.groupby(['Chromosome', 'Strand', 'Gene_ID'])
    print("Number of groups:", grouped.ngroups)

    filtered_df_list = []

    for name, group in grouped:
        total_counts = group['Counts'].sum()
        group['Relative_Frequency'] = group['Counts'] / total_counts if total_counts > 0 else 0
        
        # filter by relative frequency 
        filtered_df_list.append(group[group['Relative_Frequency'] > frequency_threshold])

    # concatenate all filtered groups into a single DataFrame
    filtered_df = pd.concat(filtered_df_list, ignore_index=True)

    # report on intial filtering 
    filtered_count = len(filtered_df)
    filtered_out_count = initial_count - filtered_count
    print(f"Filtered out {filtered_out_count} introns, retaining {filtered_count} introns based on the gene frequency threshold of {frequency_threshold}.")

    re_grouped = filtered_df.groupby(['Chromosome', 'Strand', 'Gene_ID'])

    # multiprocessing pool
    pool = multiprocessing.Pool(num_cores)
    func = partial(process_group, overlap_threshold=overlap_threshold)
    results = pool.map(func, [group for _, group in re_grouped])
    pool.close()
    pool.join()

    concatenated = pd.concat(results)

    # write results to tmpfile 
    # concatenated.to_csv("/home/150/as7425/toLocal/concat_junctions.bed", sep="\t", index=False)

    # filter based on relative frequency of reciprocally-overlapping introns 
    to_remove = set()
    for idx, row in concatenated.iterrows():
        if pd.notna(row.get('Overlapping_Index')):
            
            overlap_idx = row['Overlapping_Index']
            if overlap_idx in concatenated.index:
            
                freq_main = row['Relative_Frequency']
                
                freq_overlap = concatenated.at[overlap_idx, 'Relative_Frequency']
                # print(f"freq_main is {freq_main} and freq_overlap is {freq_overlap} and index is {idx} and overlap index is {overlap_idx}")
                 
                if freq_main / freq_overlap > 5:
                    
                    # print(f"freq_main is {freq_main} and freq_overlap is {freq_overlap} and index is {idx} and overlap index is {overlap_idx} and ratio is {freq_main / freq_overlap}")
                    
                    to_remove.add(overlap_idx)

    # print(f"list to remove is {list(to_remove)}")
    filtered_final = concatenated.drop(list(to_remove))

    print(f"Filtered out {len(to_remove)} introns based on overlap threshold and frequency disparity.")
    print(f"Remaining introns: {len(filtered_final)}")

    filtered_final = filtered_final.reindex(columns=['Chromosome', 'Start', 'End', 'Gene_ID', 'Counts', 'Strand'])

    return filtered_final

# # BENCHMARKING
        
# # # Example usage
# bam_path = '/g/data/lf10/as7425/2023_mrna-biogenesis-maps/analysis/2024-05-21_ASR013-HeLa-POINT-RNA004-rep2/ASR013_HeLa-POINT-RNA004-rep2_primary_genome_alignments_modCalled.bam'
# output_path = '/scratch/lf10/as7425/ASR013_parallel_splice_junctions.tsv'
# num_cpus = multiprocessing.cpu_count()  # Use all available CPUs
# # # parallel_extract_splice_junctions(bam_path, output_path, num_cpus)

# import time

# def benchmark(bam_path, output_path, num_cpus):
#     start_time = time.time()
#     data = parallel_extract_splice_junctions(bam_path, output_path, num_cpus)
#     summarize_junctions(data, 2)
#     end_time = time.time()
#     print(f"Time taken with {num_cpus} CPU(s): {end_time - start_time} seconds")

# # # Benchmark for 1 CPU
# # benchmark(bam_path, 'output_1_cpu.tsv', 1)

# # # Benchmark for all CPUs
# benchmark(bam_path, 'output_all_cpus.tsv', multiprocessing.cpu_count())
