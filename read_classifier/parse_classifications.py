import os
import pandas as pd
import logging

# This script processes RNA sequencing data to classify reads based on gene biotype information.
# It extracts relevant information from a provided tab-delimited file, enriches the data with gene biotypes,
# and classifies the reads into categories based on predefined criteria. It also handles and reports any potential errors.

def parse_read_classification(file_path, biotypes):
    """
    Parse the initial read classification from a given file, enriching it with biotype information.

    Parameters:
    - file_path (str): Path to the input file containing read data.
    - biotypes (dict): A dictionary mapping gene_ids to their biotype information.

    Returns:
    - DataFrame: A DataFrame with the enriched and classified read data.
    """
    
    # get the sample name
    sample_name = os.path.basename(file_path).split('_')[0]

    # read the data
    df = pd.read_csv(file_path, sep="\t", low_memory=False)

    # Add biotype, gene_name, and exon_count information
    df['biotypes_info'] = df['gene_id'].map(biotypes)
    df['gene_biotype'] = df['biotypes_info'].apply(lambda x: x.get('biotype') if isinstance(x, dict) else None)
    df['gene_name'] = df['biotypes_info'].apply(lambda x: x.get('gene_name') if isinstance(x, dict) else None)
    df['exon_count'] = df['biotypes_info'].apply(lambda x: x.get('exon_count') if isinstance(x, dict) else None)

    # add a new column for sample name
    df['sample'] = sample_name

    # add classification
    result = df.apply(lambda row: classify_read(row), axis=1)
    df['classification'], df['transcription_classification'], df['splicing_classification'] = zip(*result)

    # report rows with errors in classifiyng the absolute distance to polyA site 
    # likely caused by positions at the ends of chromosomes?
    error_df = df[df['classification'].isin(['ERROR_MARKER', 'EXCEPTION_MARKER'])]
    print(f"Rows with exceptions or errors: {len(error_df)}")

    # Summarize the manually collected error cases
    summarize_error_cases()

    # Check for missing mappings in the df 
    missing_mappings = df['biotypes_info'].isnull().sum()
    total_records = len(df)
    logging.warning(f"{missing_mappings} out of {total_records} records have missing gene_id mappings.")


    return df


# global list to collect error cases for function classify_read
error_cases = []

# function to summarise error cases for function classify_read
def summarize_error_cases():
    num_error_cases = len(error_cases)
    if num_error_cases > 0:
        print(f"Found {num_error_cases} error cases.")
        print("Preview of error cases:")
        print(error_cases[:10])  # preview first 10 error cases


# classify reads for transcription and processing based on the thresholds 
def classify_read(row, dog_threshold=25, polya_dist_threshold=25, intronic_threshold=9):
    
    try:
        # initialize transcription, splicing status 
         
        transcription_status = ''
        splicing_status = ''

        # Updated intergenic classification
        if row['gene_id'] == '0' or isna(row['gene_id']):
            return ('intergenic', 'N/A', 'N/A')

        # Insert print to check for NaN
        if pd.isna(row['exon_count']):
            print(f"Row with NaN exon_count: {row}")
        
        # Explicitly cast thresholds to integers for subsequent comparisons 
        dog_threshold = int(dog_threshold)
        polya_dist_threshold = int(polya_dist_threshold)
        intronic_threshold = int(intronic_threshold)
        
        # calculate minimum absolute distance to annotated polyA sites, handling NA cases 
        if pd.notna(row['downstream_distance']) and pd.notna(row['upstream_distance']):
            min_abs_distance = int(min(abs(int(row['downstream_distance'])), abs(int(row['upstream_distance']))))
        elif pd.notna(row['downstream_distance']):
            min_abs_distance = int(abs(int(row['downstream_distance'])))
        elif pd.notna(row['upstream_distance']):
            min_abs_distance = int(abs(int(row['upstream_distance'])))
        else:
            # Handle the case when both distance values are NaN
            error_cases.append(row)  # Append the row to error_cases
            return ('unclassified', 'unclassified', 'unclassified')

        # Define biotype sets
        single_exon_biotypes = {'rRNA', 'snRNA', 'snoRNA', 'tRNA', 'miRNA'}
        multi_exon_biotypes = {'protein_coding', 'ncRNA', 'pre_miRNA', 'pseudogene', 'transposable_element'}
        
        # Determine prefix only for single-exon genes or biotypes
        splice_count = int(row['splice_count'])
        exon_count = int(row['exon_count'])
        prefix = 'noncanonical_splicing_' if (splice_count >= 1 and 
                                            (exon_count == 1 or row['gene_biotype'] in single_exon_biotypes)) else ''
        
        # Check if the gene is single exon
        is_single_exon = exon_count == 1 or row['gene_biotype'] in single_exon_biotypes
        
        # Single exon biotypes or other single exon genes classification
        dog_overlap = int(row['DOG_overlap'])
        if is_single_exon:
            splicing_status = 'noncanonical_splicing' if prefix else 'N/A'
            
            if dog_overlap > dog_threshold and min_abs_distance > polya_dist_threshold:
                transcription_status = 'readthrough'
            elif dog_overlap > dog_threshold and min_abs_distance <= polya_dist_threshold:
                transcription_status = 'readthrough_cleaved'
            elif dog_overlap < dog_threshold and min_abs_distance < polya_dist_threshold:
                transcription_status = 'processed'
            elif dog_overlap == 0 and min_abs_distance > polya_dist_threshold:
                transcription_status = 'partially_transcribed'
            else:
                transcription_status = 'single-exon_unclassified'
                
            return (f"{prefix}{transcription_status}", transcription_status, splicing_status)

        # Multi-exon biotypes classification
        elif row['gene_biotype'] in multi_exon_biotypes:
            # Determine transcription status
            if dog_overlap > dog_threshold and min_abs_distance > polya_dist_threshold:
                transcription_status = 'readthrough'
            elif dog_overlap > dog_threshold and min_abs_distance <= polya_dist_threshold:
                transcription_status = 'readthrough_cleaved'
            elif dog_overlap <= dog_threshold and min_abs_distance <= polya_dist_threshold:
                transcription_status = 'cleaved'
            elif dog_overlap == 0 and min_abs_distance >= polya_dist_threshold:
                transcription_status = 'partially_transcribed'
            else:
                transcription_status = 'multiexon_unclassified'
            
            # Determine splicing status
            intronic_alignment = int(row['intronic_alignment'])
            if splice_count > 0 and intronic_alignment < intronic_threshold:
                splicing_status = 'spliced'
            elif splice_count > 0 and intronic_alignment >= intronic_threshold:
                splicing_status = 'partially_spliced'
            elif splice_count == 0 and intronic_alignment < intronic_threshold:
                splicing_status = 'ambiguous'
            elif splice_count == 0 and intronic_alignment > intronic_threshold:
                splicing_status = 'unspliced'
            else:
                splicing_status = 'splicing_unclassified'
            
            return (f"{transcription_status}_{splicing_status}", transcription_status, splicing_status)
        
        else:
            return 'unclassified'
        
    except Exception as e:
        logging.error(f"Exception in row {row}: {e}")
        return ('unclassified', 'unclassified', 'unclassified')
    


# print a summary of the classification 
def print_classification_summary(result_df):
    total_count = result_df.shape[0]
    if total_count == 0:
        print("No data to summarize.")
        return

    print("\nRead Classification Summary Stats:")
    print(result_df['classification'].value_counts())

    # updated classification categories
    single_exon_classifications = ['partially_transcribed', 'processed', 'readthrough', 
                                   'readthrough_cleaved',
                                   'noncanonical_splicing_partially_transcribed',
                                   'noncanonical_splicing_processed',
                                   'noncanonical_splicing_readthrough', 
                                   'single-exon_unclassified']

    multi_exon_classifications = ['readthrough_spliced', 'readthrough_partially_spliced',
                                'readthrough_ambiguous', 'readthrough_unspliced', 'readthrough_cleaved',
                                'cleaved_spliced', 'cleaved_partially_spliced',
                                'cleaved_ambiguous', 'cleaved_unspliced',
                                'partially_transcribed_spliced', 'partially_transcribed_partially_spliced',
                                'partially_transcribed_ambiguous', 'partially_transcribed_unspliced',
                                'multiexon_unclassified']  


    def filter_and_count(reads, classifications):
        filtered = reads[reads['classification'].isin(classifications)]
        return filtered.shape[0], filtered['classification'].value_counts()

    single_exon_count, single_exon_counts = filter_and_count(result_df, single_exon_classifications)
    multi_exon_count, multi_exon_counts = filter_and_count(result_df, multi_exon_classifications)
    unclassified_count = total_count - (single_exon_count + multi_exon_count)

    def calc_percentage(numerator, denominator):
        return round(numerator / denominator * 100, 1) if denominator != 0 else 0.0

    def print_stats(title, count, total, classifications, counts):
        print(f"\n{title} Summary Stats:")
        print(f"percentage of {title.lower()} reads: {calc_percentage(count, total)}%")
        for classification in classifications:
            c = counts.get(classification, 0)
            print(f"{classification}: {c} ({calc_percentage(c, count)}%)")

    print_stats("Single Exon", single_exon_count, total_count, single_exon_classifications, single_exon_counts)
    print_stats("Multi Exon", multi_exon_count, total_count, multi_exon_classifications, multi_exon_counts)

    print("\nUnclassified Reads Summary Stats:")
    print(f"percentage of unclassified reads: {calc_percentage(unclassified_count, total_count)}%")

