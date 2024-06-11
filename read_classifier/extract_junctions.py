import pandas as pd 
import pysam
import multiprocessing
from functools import partial

def extract_region_splice_junctions(bam_file, region):
    """
    Extract splice junctions for a specific region from a BAM file.
    """
    bam = pysam.AlignmentFile(bam_file, "rb")
    results = []
    for read in bam.fetch(region=region):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        current_pos = read.reference_start
        for length, operation in read.cigartuples:
            if operation in {0, 7, 8}:  # match
                current_pos += length
            elif operation == 3:  # splice
                start = current_pos
                end = current_pos + length
                strand = '-' if read.is_reverse else '+'
                results.append([read.reference_name, start, end, strand])
                current_pos += length
            elif operation == 2:  # deletion or softclip
                current_pos += length
            elif operation == 1:  # insertion 
                continue
    bam.close()
    return results

def parallel_extract_splice_junctions(bam_file, output_file, num_cpus):
    """
    Extract splice junctions from a BAM file in parallel using multiple CPUs.
    """
    bam = pysam.AlignmentFile(bam_file, "rb")
    # Get the list of contigs
    contigs = bam.references
    bam.close()

    # Create a pool of workers
    with multiprocessing.Pool(num_cpus) as pool:
        func = partial(extract_region_splice_junctions, bam_file)
        results = pool.map(func, contigs)

    # Flatten the list of lists and create a DataFrame
    flat_results = [item for sublist in results for item in sublist]
    df = pd.DataFrame(flat_results, columns=['Chromosome', 'Start', 'End', 'Strand'])
    return df

    # # Write results to file
    # with open(output_file, "w") as out:
    #     out.write("chromosome\tstart\tend\tstrand\n")
    #     for result in results:
    #         out.writelines(result)


def summarize_junctions(df, min_reads):
    """
    Summarize and filter splice junctions based on a minimum number of supporting reads.

    Parameters:
        df (pd.DataFrame): DataFrame containing splice junctions data with columns 'Chromosome', 'Start', 'End', 'Strand'.
        min_reads (int): Minimum number of reads supporting a junction to be retained.

    Returns:
        pd.DataFrame: Filtered DataFrame of unique junctions meeting the read count criterion.
    """
    print(f"{len(df)} junctions")
    print(f"{df.drop_duplicates().shape[0]} unique junctions")

    # group by junction characteristics and count occurrences
    grouped = df.groupby(['Chromosome', 'Start', 'End', 'Strand']).size().reset_index(name='Counts')
    
    # filter based on the count of reads supporting each junction
    filtered = grouped[grouped['Counts'] >= min_reads]
    
    print(f"{filtered.shape[0]} junctions after filtering")
    print(f"{filtered.drop_duplicates().shape[0]} unique junctions after filtering")

    # drop the 'Counts' column as it's no longer needed, returning only the unique junctions
    return filtered.drop(columns=['Counts'])


# ###########
        
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
