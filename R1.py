import sys
import os
import warnings
import re 
import argparse
import subprocess
import pandas as pd
from pandas import isna 
import numpy as np

from tempfile import NamedTemporaryFile

from concurrent.futures import ThreadPoolExecutor, as_completed, ProcessPoolExecutor
import concurrent.futures
from scipy.cluster.hierarchy import fclusterdata
from typing import Tuple
from functools import partial
import csv 

from typing import Dict

### local imports 
sys.path.append("/home/150/as7425/R1/read_classifier/")
from bam_to_bed import bam_to_bed
from process_genome import generate_genome_file, read_chromosomes_from_genome_file, filter_bed_by_chromosome_inplace
from gtf_to_bed import gtf_to_bed, extend_gene_bed 
from parallel_bed_operations import split_bed_file, write_sorted_chunk
from intersect import run_bedtools
from parse_intersect import parse_output
from process_read_end_positions import calculate_distance_for_unique_positions, calculate_distance_to_annotated_ends
from process_gtf import parse_attributes, parse_gtf_line, get_gene_exon_table, get_biotypes 
from parse_classifications import parse_read_classification, summarize_error_cases, classify_read, print_classification_summary

### mod table 
sys.path.append("/home/150/as7425/R1/parse_modifications_data/")
from create_mod_table import extract_modifications

# set logging level 
import logging 
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

# main 
def main(bam_file, gtf_file, output_table):

    # fetch output directory 
    output_dir = os.path.dirname(output_table) 

    # generate a genome file for bedtools sort 
    genome_file = generate_genome_file(bam_file, output_dir) 

    # Extract modifications from the BAM file
    mod_output_file = os.path.join(output_dir, "modifications.csv")
    with open(mod_output_file, 'w') as mod_file:
        # Write header to the CSV file
        mod_file.write("ReadID,Reference,Strand,Modifications\n")
        
        # Iterate through the modifications yielded by extract_modifications
        for modifications in extract_modifications(bam_file):
            if modifications is not None:  # Check if modifications is not None
                mod_file.write(f"{modifications}\n")  # Write modifications to file
                print(modifications)  # Optionally, also print modifications
            else:
                # This else block can be used if you want to log no modification cases
                # For example, you could write a specific line, or simply not handle it here
                print("No modifications or missing tags for some reads.")

    return 

    # create a bed file of gene regions
    gene_bed_file = gtf_to_bed(gtf_file, "gene", output_dir) 

    # filter for relevant chromsoomes
    chromosomes = read_chromosomes_from_genome_file(genome_file)
    filter_bed_by_chromosome_inplace(gene_bed_file, chromosomes)

    # create a bed file of downstream of gene regions 
    dog_bed_file = extend_gene_bed(gene_bed_file, output_dir, genome_file) 

    # intersect reads against, exons, genes, and downstream of gene regions 
    intersect_files = run_bedtools(bam_file, gtf_file, dog_bed_file, genome_file, output_dir) 

    # parse the bedtools intersect files 
    df = parse_output(intersect_files[0], intersect_files[1], intersect_files[2], intersect_files[3], intersect_files[4])
    
    # Write df to a CSV file for testing
    df.to_csv(output_dir + "/test_df_end.csv", index=False)
    
    # fetch distances to nearest transcript end sites 
    end_position_df = calculate_distance_to_annotated_ends(df, gtf_file, output_dir)

    # save the output as CSV 
    end_position_df.to_csv(output_table, sep="\t", index=False)

    # fetch the gene biotypes 
    biotypes = get_biotypes(gtf_file)

    # classify each read 
    result_df = parse_read_classification(output_table, biotypes)

    # print some summary statistics 
    print_classification_summary(result_df)

    # Rename 'exon_count' to 'annotated_exon_count'
    result_df.rename(columns={'exon_count': 'annotated_exon_count'}, inplace=True)

    # Omit specific columns
    columns_to_omit = ['read_end_strand', 'strand_sign', 'chromosome', 'position', 'strand', 'biotypes_info']
    result_df.drop(columns=columns_to_omit, inplace=True)

    # Save the table to output
    result_df.to_csv(output_table, sep='\t', index=False)

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Classify chromatin-bound direct RNA-Seq reads for splicing and transcriptional status')
    parser.add_argument('-b', '--bam_file', type=str, help='Path to the input bam file.')
    parser.add_argument('-g', '--gtf_file', type=str, help='Path to the GTF file.')
    parser.add_argument('-o', '--output_table', type=str, default=None, help='Optional path to the output file.')

    args = parser.parse_args()

    bam_file = args.bam_file
    gtf_file = args.gtf_file 
    output_table = args.output_table

    main(bam_file, gtf_file, output_table)


