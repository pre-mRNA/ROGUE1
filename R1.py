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
from map_mm_tag_to_genome_position import extract_modifications
from extract_polyA_length import fetch_polyA_pt

# set logging level 
import logging 
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

# main 
def main(bam_file, gtf_file, output_table, calculate_modifications, calculate_polyA):

    # fetch output directory 
    output_dir = os.path.dirname(output_table) 

    # generate a genome file for bedtools sort 
    genome_file = generate_genome_file(bam_file, output_dir) 

    # calculate mod probs if flag is set 
    if calculate_modifications:
        modifications_dict = {}
        for modification in extract_modifications(bam_file):
            if modification is not None:
                read_id, reference, strand, mods = modification.split('\t')
                modifications_dict[read_id] = mods
            else:
                modifications_dict[read_id] = "No modifications or missing tags"
        modifications_df = pd.DataFrame(list(modifications_dict.items()), columns=['read_id', 'Modifications'])
        print(modifications_df.head())

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
    
    # merge the modifications if mods are calculated 
    if calculate_modifications:
        df = pd.merge(df, modifications_df, on='read_id', how='left')

    # merge the polyA tail length calls if requested 
    if calculate_polyA:
        polyA_lengths = fetch_polyA_pt(bam_file)
        polyA_df = pd.DataFrame(list(polyA_lengths.items()), columns=['read_id', 'PolyA_Length'])
        df = pd.merge(df, polyA_df, on='read_id', how='left')

    # write df to a CSV file for testing
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

    # omit specific columns from output 
    # columns_to_omit = ['read_end_strand', 'strand_sign', 'chromosome', 'position', 'strand', 'biotypes_info']
    # keep biotypes for now?
    columns_to_omit = ['read_end_strand', 'strand_sign', 'chromosome', 'position', 'strand']  
    result_df.drop(columns=columns_to_omit, inplace=True)

    # Save the table to output
    result_df.to_csv(output_table, sep='\t', index=False)

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Classify chromatin-bound direct RNA-Seq reads for splicing and transcriptional status')
    parser.add_argument('-b', '--bam_file', type=str, help='Path to the input bam file.')
    parser.add_argument('-g', '--gtf_file', type=str, help='Path to the GTF file.')
    parser.add_argument('-o', '--output_table', type=str, default=None, help='Optional path to the output file.')
    parser.add_argument('-m', '--modifications', action='store_true', help='Calculate modifications data if this flag is present.')
    parser.add_argument('-p', '--polyA', action='store_true', help='Calculate polyA tail lengths if this flag is present.')



    args = parser.parse_args()

    bam_file = args.bam_file
    gtf_file = args.gtf_file 
    output_table = args.output_table

    main(args.bam_file, args.gtf_file, args.output_table, args.modifications, args.polyA)

