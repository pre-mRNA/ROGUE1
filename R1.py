import sys
import os
import argparse
import logging 
import pandas as pd
import multiprocessing
import tempfile

# function imports 
sys.path.append("/home/150/as7425/R1/read_classifier/")
from process_genome import generate_genome_file, read_chromosomes_from_genome_file, filter_bed_by_chromosome_inplace
from gtf_to_bed import gtf_to_bed, extend_gene_bed, sort_bed_file
from intersect import run_bedtools
from parse_intersect import parse_output
from process_read_end_positions import calculate_distance_to_read_ends, get_transcript_ends
from process_gtf import get_biotypes 
from parse_classifications import parse_read_classification, print_classification_summary
from extract_junctions import parallel_extract_splice_junctions, summarize_junctions, filter_junctions 

sys.path.append("/home/150/as7425/R1/parse_modifications_data/")
from map_mm_tag_to_genome_position import extract_modifications
from extract_polyA_length import fetch_polyA_pt

# set logging level 
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

# get cpu count
cpu_count = multiprocessing.cpu_count()

def main(bam_file, gtf_file, output_table, calculate_modifications, calculate_polyA, junction_distance, index_path):

    output_dir = os.path.dirname(output_table) 

    # check if the output directory exists and is writeable and print to log
    if not os.path.exists(output_dir):
        logging.info(f"Output directory {output_dir} does not exist. Creating it now.")
        os.makedirs(output_dir)
    else:
        logging.info(f"Output directory for ROGUE1 table is {output_dir}.")

    # generate a genome file for bedtools sort 
    genome_file = generate_genome_file(bam_file, output_dir) 

    # create a bed file of gene regions
    gene_bed_file = gtf_to_bed(gtf_file, "gene", output_dir) 

    # filter for relevant chromsoomes
    chromosomes = read_chromosomes_from_genome_file(genome_file)
    filter_bed_by_chromosome_inplace(gene_bed_file, chromosomes)

    # create a bed file of downstream of gene regions 
    dog_bed_file = extend_gene_bed(gene_bed_file, output_dir, genome_file) 

    # if index_path is provided, read BED files from the index
    if index_path:
        pas_bed = os.path.join(index_path, "PAS.bed")
        dog_bed = os.path.join(index_path, "updated_downstream_of_gene.bed")
        exon_bed = os.path.join(index_path, "updated_exon.bed")
        intron_bed = os.path.join(index_path, "updated_intron.bed")
        gene_bed = os.path.join(index_path, "updated_gene.bed")
        
        # Check if all required files exist
        required_files = [pas_bed, dog_bed, exon_bed, intron_bed, gene_bed]
        if all(os.path.exists(f) for f in required_files):
            logging.info(f"Using index files from {index_path}")
            
            # Create temporary directory for sorted files
            with tempfile.TemporaryDirectory() as temp_dir:
                # Sort all index files
                sorted_pas_bed = sort_bed_file(pas_bed, os.path.join(temp_dir, "sorted_PAS.bed"))
                sorted_dog_bed = sort_bed_file(dog_bed, os.path.join(temp_dir, "sorted_updated_downstream_of_gene.bed"))
                sorted_exon_bed = sort_bed_file(exon_bed, os.path.join(temp_dir, "sorted_updated_exon.bed"))
                sorted_intron_bed = sort_bed_file(intron_bed, os.path.join(temp_dir, "sorted_updated_intron.bed"))
                sorted_gene_bed = sort_bed_file(gene_bed, os.path.join(temp_dir, "sorted_updated_gene.bed"))
                
                # Update the paths to use the sorted files
                pas_bed, dog_bed, exon_bed, intron_bed, gene_bed = sorted_pas_bed, sorted_dog_bed, sorted_exon_bed, sorted_intron_bed, sorted_gene_bed
                

                # intersect reads against, exons, genes, and downstream of gene regions 
                intersect_files = run_bedtools(bam_file, gtf_file, genome_file, output_dir, 
                                            pas_bed=pas_bed, dog_bed=dog_bed, exon_bed=exon_bed, 
                                            intron_bed=intron_bed, gene_bed=gene_bed)

                # parse the bedtools intersect files 
                df = parse_output(intersect_files[0], intersect_files[1], intersect_files[2], intersect_files[3], intersect_files[4])
                
                # store a dict of read_id and gene_id 
                gene_id_map = df.set_index('read_id')['gene_id'].to_dict()

                # extract and merge the modifications if mods are calculated 
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
                    df = pd.merge(df, modifications_df, on='read_id', how='left')

                # merge the polyA tail length calls if requested 
                if calculate_polyA:
                    polyA_lengths = fetch_polyA_pt(bam_file)
                    polyA_df = pd.DataFrame(list(polyA_lengths.items()), columns=['read_id', 'polya_length'])
                    df = pd.merge(df, polyA_df, on='read_id', how='left')

                # calculate distance to junctions if requested 
                if junction_distance: 

                    # extract all junctions 
                    all_junctions = parallel_extract_splice_junctions(bam_file, num_cpus = cpu_count, gene_id_map=gene_id_map)
                    all_junctions.to_csv(output_dir + "/all_junctions.bed", sep="\t", index=False)

                    # collapse to unique junctions with support level of at least 2 alignments 
                    collapsed_junctions = summarize_junctions(all_junctions, 2)

                    # save junctions to bed file for future inspection 
                    collapsed_junctions.to_csv(output_dir + "/collapsed_junctions.bed", sep="\t", index=False)

                    # filter junctions
                    final_junctions = filter_junctions(collapsed_junctions)

                    # save junctions 
                    with open(output_dir + "/final_junctions.bed", 'w') as file:
                        file.write("track name=junctions\n")
                        final_junctions.to_csv(file, sep='\t', index=False, header=False)

                    # convert junctions to donor final nucelotides 
                    splice_donors = final_junctions.assign(
                        Start=lambda df: df.apply(lambda x: x['Start'] - 1 if x['Strand'] == '+' else x['End'], axis=1),
                        End=lambda df: df.apply(lambda x: x['Start'] if x['Strand'] == '+' else x['End'] + 1, axis=1)
                    ).rename(columns={'Chromosome': 'chromosome', 'Start': 'position', 'Strand': 'strand'})

                    # calculate distances to nearest splice donors 
                    df = calculate_distance_to_read_ends(df, splice_donors, "splice_donors")

    # fetch the gene biotypes 
    biotypes = get_biotypes(gtf_file)
    biotypes_df = pd.DataFrame(list(biotypes.items()), columns=['gene_id', 'biotype_info'])
    # print(biotypes_df.head())
    df = pd.merge(df, biotypes_df, on='gene_id', how='left')
    df['biotype_info'] = df['biotype_info'].apply(lambda x: x if isinstance(x, dict) else {})
    df_biotype_expanded = pd.json_normalize(df['biotype_info'])
    df = pd.concat([df.drop('biotype_info', axis=1), df_biotype_expanded], axis=1)
    df.fillna({'biotype': 'unknown', 'gene_name': 'unknown', 'exon_count': 0}, inplace=True)

    # calculate distances between read ends and nearest transcript end sites 
    transcript_ends = get_transcript_ends(gtf_file, output_dir)
    end_position_df = calculate_distance_to_read_ends(df, transcript_ends, "polyA")

    # save the output as CSV 
    end_position_df.to_csv(output_table, sep="\t", index=False)

    # classify each read 
    # result_df = parse_read_classification(output_table, biotypes)

    # print some summary statistics 
    # print_classification_summary(result_df)

    # # Rename 'exon_count' to 'annotated_exon_count'
    # result_df.rename(columns={'exon_count': 'annotated_exon_count'}, inplace=True)

    # # omit specific columns from output 
    # # columns_to_omit = ['read_end_strand', 'strand_sign', 'chromosome', 'position', 'strand', 'biotypes_info']
    # # keep biotypes for now?
    # columns_to_omit = ['read_end_strand', 'strand_sign', 'chromosome', 'position', 'strand']  
    # result_df.drop(columns=columns_to_omit, inplace=True)

    # # Save the table to output
    # result_df.to_csv(output_table, sep='\t', index=False)

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Classify chromatin-bound direct RNA-Seq reads for splicing and transcriptional status')
    parser.add_argument('-b', '--bam_file', type=str, help='Path to the input bam file.')
    parser.add_argument('-g', '--gtf_file', type=str, help='Path to the GTF file.')
    parser.add_argument('-o', '--output_table', type=str, default=None, help='Optional path to the output file.')
    parser.add_argument('-m', '--modifications', action='store_true', help='Calculate modifications data if this flag is present.')
    parser.add_argument('-p', '--polyA', action='store_true', help='Calculate polyA tail lengths if this flag is present.')
    parser.add_argument('-j', '--junction_distance', action='store_true', help='Calculate distance between read 3 prime end and splice donors.')
    parser.add_argument('--index', type=str, help='Path to the index directory containing pre-computed BED files')

    args = parser.parse_args()

    bam_file = args.bam_file
    gtf_file = args.gtf_file 
    output_table = args.output_table
    
    
    main(args.bam_file, args.gtf_file, args.output_table, args.modifications, args.polyA, args.junction_distance, args.index)

