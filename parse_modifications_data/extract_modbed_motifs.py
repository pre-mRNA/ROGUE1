import sys
import argparse
from pyfaidx import Fasta
from multiprocessing import Pool, cpu_count
from collections import Counter, defaultdict
import logging
import time
import subprocess

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def process_chromosome_chunk(args):
    chrom, bed_lines, fasta_path = args
    results = []
    with Fasta(fasta_path) as fasta:
        for line in bed_lines:
            fields = line.strip().split()
            start, strand = int(fields[1]), fields[5]
            seq_region = str(fasta[chrom][start-2:start+3]).upper()  # Convert to uppercase here
            if strand == '-':
                seq_region = reverse_complement(seq_region)
            results.append(f"{line.strip()}\t{seq_region}\n")
    return results

def parallel_sort(input_file, output_file, num_cpus):
    logging.info("Performing parallel sort on the output file...")
    start_time = time.time()
    cmd = f"sort -k1,1 -k2,2n --parallel={num_cpus} -o {output_file} {input_file}"
    subprocess.run(cmd, shell=True, check=True)
    logging.info(f"Parallel sort completed in {time.time() - start_time:.2f} seconds")

def main(args):
    overall_start_time = time.time()

    # group intervals by chromosome
    logging.info("Grouping BED lines by chromosome...")
    chrom_groups = defaultdict(list)
    with open(args.ibed, 'r') as f:
        for line in f:
            chrom = line.split()[0]
            chrom_groups[chrom].append(line)

    logging.info(f"Processing {sum(len(lines) for lines in chrom_groups.values())} lines across {len(chrom_groups)} chromosomes...")

    # process chromosomes in parallel
    start_time = time.time()
    num_cpus = cpu_count()
    with Pool(processes=num_cpus) as pool:
        results = pool.map(process_chromosome_chunk, 
                           [(chrom, lines, args.fasta) for chrom, lines in chrom_groups.items()])
    logging.info(f"Processing completed in {time.time() - start_time:.2f} seconds")

    # write results to temp file
    temp_output = args.obed + ".temp"
    logging.info("Writing output to temporary file...")
    start_time = time.time()
    with open(temp_output, 'w') as out:
        for chunk in results:
            out.writelines(chunk)
    logging.info(f"Output written in {time.time() - start_time:.2f} seconds")

    # parallel sort on the temporary file
    parallel_sort(temp_output, args.obed, num_cpus)

    # rm temporary file
    subprocess.run(f"rm {temp_output}", shell=True, check=True)

    # calculate central nucleotide frequency distribution
    logging.info("Calculating frequency distribution...")
    start_time = time.time()
    central_nucleotides = []
    error_count = 0
    with open(args.obed, 'r') as f:
        for line_number, line in enumerate(f, 1):
            try:
                central_nt = line.strip().split('\t')[-1][2]
                central_nucleotides.append(central_nt)
            except IndexError:
                central_nucleotides.append('N')  # or any other placeholder
                error_count += 1
                logging.warning(f"Error processing line {line_number}: {line.strip()}")

    counter = Counter(central_nucleotides)
    total = sum(counter.values())
    frequencies = {nt: count / total * 100 for nt, count in counter.items()}
    logging.info("Frequency distribution of central nucleotide:")
    
    for nt, freq in sorted(frequencies.items()):
        logging.info(f"{nt}: {freq:.2f}%")
    logging.info(f"Frequency calculation completed in {time.time() - start_time:.2f} seconds")
    logging.info(f"Total lines with errors: {error_count}")
    logging.info(f"Total execution time: {time.time() - overall_start_time:.2f} seconds")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process BED and FASTA files to find and report 5-mer sequences.")
    parser.add_argument("-ibed", type=str, required=True, help="Input BED file path.")
    parser.add_argument("-fasta", type=str, required=True, help="Input FASTA file path.")
    parser.add_argument("-obed", type=str, required=True, help="Output BED file path with 5-mer sequences appended.")
    args = parser.parse_args()
    
    main(args)