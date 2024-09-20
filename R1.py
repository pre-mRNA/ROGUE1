import sys
import os
import argparse
import logging 
import pandas as pd
import multiprocessing
import tempfile
import atexit
import shutil

# function imports 
sys.path.append("/home/150/as7425/R1/read_classifier/")
from process_genome import generate_genome_file, read_chromosomes_from_genome_file, filter_bed_by_chromosome_inplace
from gtf_to_bed import gtf_to_bed, extend_gene_bed, sort_bed_file, number_exons_and_introns_in_bed
from intersect import run_bedtools
from parse_intersect import parse_output
from process_read_end_positions import calculate_distance_to_read_ends, get_transcript_ends
from fetch_read_end_feature import get_end_feature
from process_gtf import get_biotypes 
from extract_junctions import parallel_extract_splice_junctions, summarize_junctions, filter_junctions, process_intron_to_junctions 
from polyA_z_score import classify_polya_reads

sys.path.append("/home/150/as7425/R1/parse_modifications_data/")
from map_mm_tag_to_genome_position import extract_modifications
from extract_polyA_length import fetch_polyA_pt

sys.path.append("/home/150/as7425/R1/classify_splicing")
from canonical_splicing_from_acceptor import classify_splicing

# set logging level 
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

# get cpu count
cpu_count = multiprocessing.cpu_count()

# prepare to clean up temp files 
temp_dir = tempfile.mkdtemp()

def cleanup_temp_dir():
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)

atexit.register(cleanup_temp_dir)

# load a pre-made genome index 
def load_index_files(index_path, temp_dir):
    pas_bed = os.path.join(index_path, "PAS.bed")
    dog_bed = os.path.join(index_path, "updated_downstream_of_gene.bed")
    exon_bed = os.path.join(index_path, "updated_exon.bed")
    intron_bed = os.path.join(index_path, "updated_intron.bed")
    gene_bed = os.path.join(index_path, "updated_gene.bed")

    required_files = [pas_bed, dog_bed, exon_bed, intron_bed, gene_bed]
    if all(os.path.exists(f) for f in required_files):
        logging.info(f"Using index files from {index_path}")
        sorted_pas_bed = sort_bed_file(pas_bed, os.path.join(temp_dir, "sorted_PAS.bed"))
        sorted_dog_bed = sort_bed_file(dog_bed, os.path.join(temp_dir, "sorted_updated_downstream_of_gene.bed"))
        sorted_exon_bed = sort_bed_file(exon_bed, os.path.join(temp_dir, "sorted_updated_exon.bed"))
        sorted_intron_bed = sort_bed_file(intron_bed, os.path.join(temp_dir, "sorted_updated_intron.bed"))
        sorted_gene_bed = sort_bed_file(gene_bed, os.path.join(temp_dir, "sorted_updated_gene.bed"))
        return sorted_pas_bed, sorted_dog_bed, sorted_exon_bed, sorted_intron_bed, sorted_gene_bed
    else:
        return None, None, None, None, None

def main(bam_file, gtf_file, output_table, calculate_modifications, calculate_polyA, junction_distance, index_path, record_exons, co_threshold, post_threshold):
    output_dir = os.path.dirname(output_table) 

    if not os.path.exists(output_dir):
        logging.info(f"Output directory {output_dir} does not exist. Creating it now.")
        os.makedirs(output_dir)
    else:
        logging.info(f"Output directory for ROGUE1 table is {output_dir}.")

    genome_file = generate_genome_file(bam_file, output_dir)

    pas_bed, dog_bed, exon_bed, intron_bed, gene_bed = None, None, None, None, None

    # import annotations from index, if provided 
    if index_path:
        pas_bed, dog_bed, exon_bed, intron_bed, gene_bed = load_index_files(index_path, temp_dir)
        if not all([pas_bed, dog_bed, exon_bed, intron_bed, gene_bed]):
            logging.warning("Some required files are missing in the index path. Creating new index files.")
            index_path = None

    if not index_path:

        # create the gene bed file 
        gene_bed = gtf_to_bed(gtf_file, "gene", output_dir)
        
        # sort and number the exon bed file 
        unsorted_exon_bed = gtf_to_bed(gtf_file, "exon", output_dir)
        exon_bed = os.path.join(output_dir, 'numbered_exon.bed')
        number_exons_and_introns_in_bed(unsorted_exon_bed, exon_bed, 'exon') # add exon numbers to bed col5

        # sort and number the intron bed file
        unsorted_intron_bed = gtf_to_bed(gtf_file, "intron", output_dir)
        intron_bed = os.path.join(output_dir, 'numbered_intron.bed')
        number_exons_and_introns_in_bed(unsorted_intron_bed, intron_bed, 'intron') # add exon numbers to bed col5

        # create the DOG bed file 
        dog_bed = extend_gene_bed(gene_bed, output_dir, genome_file)

        # filter annotations for relevant chromosomes 
        chromosomes = read_chromosomes_from_genome_file(genome_file)
        for bed_file in [gene_bed, exon_bed, intron_bed, dog_bed]:
            filter_bed_by_chromosome_inplace(bed_file, chromosomes)

    # intersect aligned regions against annotations 
    intersect_files = run_bedtools(bam_file, genome_file, output_dir, dog_bed, exon_bed, intron_bed, gene_bed, num_files=104)

    # parse intersection output 
    df = parse_output(intersect_files[0], intersect_files[1], intersect_files[2], intersect_files[3], 
                      intersect_files[4], intersect_files[5], record_exons, genome_file, output_dir, 
                      num_files=104, original_exon_bed=exon_bed, original_intron_bed=intron_bed, 
                      original_dog_bed=dog_bed)
    
    logging.info(f"Processed overlaps for {len(df)} reads")
    
    # store a dict of read_id and gene_id 
    gene_id_map = df.set_index('read_id')['gene_id'].to_dict()

    # get read end feature 
    end_features = get_end_feature(intersect_files[6], genome_file, output_dir, dog_bed, exon_bed, intron_bed, gene_id_map)

    # ensure all reads are present in end_features
    missing_reads = set(df['read_id']) - set(end_features.index)
    if missing_reads:
        logging.warning(f"{len(missing_reads)} reads are missing from end_features. Adding them with 'unclassified' feature.")
        missing_df = pd.DataFrame([{'read_id': read_id, 'gene_id': df.loc[df['read_id'] == read_id, 'gene_id'].iloc[0], 'feature': 'unclassified'} for read_id in missing_reads])
        missing_df.set_index('read_id', inplace=True)
        end_features = pd.concat([end_features, missing_df])

    df = pd.merge(df, end_features, on=['read_id', 'gene_id'], how='left')
    df['feature'] = df['feature'].fillna('unclassified')
    df.rename(columns={'feature': 'three_prime_feature'}, inplace=True)

    # reorder columns
    cols = df.columns.tolist()
    three_prime_feature_index = cols.index('three_prime_feature')
    cols = cols[:df.columns.get_loc('read_end_strand')+1] + ['three_prime_feature'] + cols[df.columns.get_loc('read_end_strand')+1:three_prime_feature_index] + cols[three_prime_feature_index+1:]
    df = df[cols]

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

        # use the junctions in the annotation if an index is provided 
        if index_path and intron_bed:
            logging.info("Using pre-computed introns from index")
            final_junctions = process_intron_to_junctions(intron_bed)
        else:

        # extract all junctions 
            logging.info("Performing de novo junction discovery from alignments")
            all_junctions = parallel_extract_splice_junctions(bam_file, num_cpus = cpu_count, gene_id_map=gene_id_map)
            all_junctions.to_csv(output_dir + "/all_junctions.bed", sep="\t", index=False)

            # collapse to unique junctions with support level of at least 2 alignments 
            logging.info("Summarising de novo junctions")
            collapsed_junctions = summarize_junctions(all_junctions, 2)

            # save junctions to bed file for future inspection 
            collapsed_junctions.to_csv(output_dir + "/collapsed_junctions.bed", sep="\t", index=False)

            # filter junctions
            logging.info("Filtering de novo junctions")
            final_junctions = filter_junctions(collapsed_junctions)

        # save junctions 
        with open(output_dir + "/final_junctions.bed", 'w') as file:
            file.write("track name=junctions\n")
            final_junctions.to_csv(file, sep='\t', index=False, header=False)

        # convert junctions to donor final nucelotides 
        splice_donors = final_junctions.assign(
            Start=lambda df: df.apply(lambda x: x['Start'] if x['Strand'] == '+' else x['End'], axis=1),
            End=lambda df: df.apply(lambda x: x['Start'] if x['Strand'] == '+' else x['End'], axis=1)
        ).rename(columns={'Chromosome': 'chromosome', 'Start': 'position', 'Strand': 'strand'})

        # calculate distances to nearest splice donors 
        df = calculate_distance_to_read_ends(df, splice_donors, "splice_donors")

    # fetch the gene biotypes 
    biotypes = get_biotypes(gtf_file)
    biotypes_df = pd.DataFrame(list(biotypes.items()), columns=['gene_id', 'biotype_info'])
    print(df.head())
    df = pd.merge(df, biotypes_df, on='gene_id', how='left')
    df['biotype_info'] = df['biotype_info'].apply(lambda x: x if isinstance(x, dict) else {})
    df_biotype_expanded = pd.json_normalize(df['biotype_info'])
    df = pd.concat([df.drop('biotype_info', axis=1), df_biotype_expanded], axis=1)
    df.fillna({'biotype': 'unknown', 'gene_name': 'unknown', 'exon_count': 0}, inplace=True)

    # classify canonical splicing from acceptors 
    if record_exons:
        logging.info("Classifying splicing status for reads...")
        df = classify_splicing(df)
        logging.info("Splicing classification complete.")

    # calculate distances between read ends and nearest transcript end sites 
    transcript_ends = get_transcript_ends(gtf_file, output_dir)
    end_position_df = calculate_distance_to_read_ends(df, transcript_ends, "polyA")

    # classify reads based on Poly(A) lengths
    if calculate_polyA:
        end_position_df = classify_polya_reads(end_position_df, output_dir, co_threshold, post_threshold)

    # save the output as CSV 
    end_position_df.to_csv(output_table, sep="\t", index=False)

    # classify each read 
    # result_df = parse_read_classification(output_table, biotypes)

    # print some summary statistics 
    # print_classification_summary(result_df)

    # # rename 'exon_count' to 'annotated_exon_count'
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
    parser.add_argument('--record_exons', action='store_true', help='Record the IDs of exons and introns that each read overlaps')
    parser.add_argument('--co_threshold', type=float, default=1.6, help='Z-score threshold for co-transcriptional classification.')
    parser.add_argument('--post_threshold', type=float, default=1.9, help='Z-score threshold for post-transcriptional classification.')

    args = parser.parse_args()
    
    main(args.bam_file, args.gtf_file, args.output_table, args.modifications, args.polyA, args.junction_distance, args.index, args.record_exons, args.co_threshold, args.post_threshold)

