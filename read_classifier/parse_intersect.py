import os
import warnings
import pandas as pd 
import numpy as np 


# parse the intersect ouput
def parse_output(exon_overlap_file, gene_overlap_file, dog_overlap_file, bed_file, end_coordinates):


    for file in [exon_overlap_file, gene_overlap_file, dog_overlap_file]:
        if not os.path.isfile(file):
            raise Exception(f"Essential file {file} does not exist. Unable to proceed.")
        elif os.stat(file).st_size == 0:
            warnings.warn(f"Warning: Essential file {file} is empty. Some data may not be processed.")
            
    # define column names for exon_overlaps and gene_overlaps
    exon_cols = ['read_chrom', 'read_fragment_start', 'read_fragment_end', 'read_id', 'read-alignment_length', 'read_strand', 'exon_chrom', 'exon_start', 'exon_end', 'exon_gene_id', 'exon_gene_name', 'exon_strand', 'exon_base_overlap']
    gene_cols = ['read_chrom', 'read_fragment_start', 'read_fragment_end', 'read_id', 'read-alignment_length', 'read_strand', 'gene_chrom', 'gene_start', 'gene_end', 'gene_id', 'gene_biotype', 'gene_strand', 'gene_base_overlap']
    bam_to_bed_cols = ['chrom', 'start', 'end', 'name', 's-a_length', 'strand']

    # read in the overlap data
    exon_df = pd.read_csv(exon_overlap_file, sep="\t", header=None, names = exon_cols, low_memory=False)
    gene_df = pd.read_csv(gene_overlap_file, sep="\t", header=None, names = gene_cols, low_memory=False)
    bam_df = pd.read_csv(bed_file, sep="\t", header=None, names = bam_to_bed_cols, low_memory=False)

    # select the most 5' overlap coordinate for each overlap
    exon_df['overlap_coord'] = exon_df.apply(lambda row: row['read_fragment_start'] if row['read_strand'] == '+' else row['read_fragment_end'], axis=1)

    # split exon overlaps into + strand and - strand reads
    positive_strand_exon_df = exon_df[exon_df['read_strand'] == '+']
    negative_strand_exon_df = exon_df[exon_df['read_strand'] == '-']

    # sort them differently, getting the most 5' overlap of the read against the exon 
    positive_strand_exon_df = positive_strand_exon_df.sort_values('overlap_coord', ascending=True)
    negative_strand_exon_df = negative_strand_exon_df.sort_values('overlap_coord', ascending=False)

    # combine them back together and select the most 5' overlap
    sorted_exon_df = pd.concat([positive_strand_exon_df, negative_strand_exon_df])
    top_exon_overlap = sorted_exon_df.drop_duplicates(subset='read_id', keep='first')

    # get the selected gene for each read based on exon overlap and store in top_exon_gene 
    top_exon_gene = {row['read_id']: row['exon_gene_id'] for _, row in top_exon_overlap.iterrows()}

    # create a DataFrame that only includes overlaps with the selected gene for each read
    selected_exon_df = exon_df[exon_df.apply(lambda row: top_exon_gene[row['read_id']] == row['exon_gene_id'], axis=1)]

    # calculate the total overlap of each read with its corresponding gene
    exon_overlap_group = selected_exon_df.groupby(['read_id', 'exon_gene_id'])['exon_base_overlap'].sum().reset_index()

    # rename the column to make it clear what it represents
    # exon_overlap_group.rename(columns={'exon_base_overlap': 'total_exon_base_overlap'}, inplace=True)

    # preview the DataFrame
    print("Preview of exon_overlap_group:")
    print(exon_overlap_group.head())

    # identify the gene to which each read has the most alignment
    gene_overlap_group = gene_df.groupby(['read_id', 'gene_id'])['gene_base_overlap'].sum().reset_index()
    gene_overlap_group.sort_values('gene_base_overlap', ascending=False, inplace=True)
    top_gene_overlap = gene_overlap_group.drop_duplicates(subset='read_id', keep='first')
    
    # create a dict to map from read_id to gene_id with highest gene overlap
    top_gene = {row['read_id']: row['gene_id'] for _, row in top_gene_overlap.iterrows()}

    # update the dict with exon information where it's available
    for read_id, exon_gene_id in top_exon_gene.items():
        top_gene[read_id] = exon_gene_id
    
    # use top_gene for filtering the gene_df
    gene_df_filtered = gene_df[gene_df.apply(lambda row: top_gene.get(row['read_id']) == row['gene_id'], axis=1)]

    # preview of gene_df_filtered
    print("Preview of gene_df_filtered:")
    print(gene_df_filtered.head())

    # now, sum the base overlap against the gene
    gene_overlap_sum = gene_df_filtered.groupby(['read_id', 'gene_id', 'read-alignment_length'])['gene_base_overlap'].sum().reset_index()

    # merge the overlap with gene and the overlap with exons, based on on the read_id and gene_id in the two datasets 
    gene_overlap_sum = pd.merge(gene_overlap_sum, exon_overlap_group, left_on=['read_id', 'gene_id'], right_on=['read_id', 'exon_gene_id'], how='left')
    gene_overlap_sum['exon_base_overlap'].fillna(0, inplace=True)

    # preview the gene_overlap sum
    print("Preview of gene_overlap_sum:")
    print(gene_overlap_sum.head())

    # find the overlap with DOG regions 
    print("Calculating DOG coverage")

    # first, read in the downstream of gene coverage data
    dog_cols = ['read_chrom', 'read_fragment_start', 'read_fragment_end', 'read_id', 'read-alignment_length', 'read_strand', 'gene_chrom', 'gene_start', 'gene_end', 'dog_gene_id', 'gene_biotype', 'gene_strand', 'dog_base_overlap']
    dog_df = pd.read_csv(dog_overlap_file, sep="\t", header=None, names = dog_cols, low_memory=False)

    # Read and handle the potentially empty DOG file
    if os.path.isfile(dog_overlap_file) and os.stat(dog_overlap_file).st_size > 0:
        dog_df = pd.read_csv(dog_overlap_file, sep="\t", header=None, names=dog_cols, low_memory=False)
    else:
        # Create an empty DataFrame with the necessary columns to ensure subsequent operations do not fail
        dog_df = pd.DataFrame(columns=dog_cols)
        print("DOG data is empty, defaulting to empty DataFrame.")

    # Initialize an empty DataFrame for dog_overlap_group if dog_df is empty
    if dog_df.empty:
        dog_overlap_group = pd.DataFrame(columns=['read_id', 'dog_gene_id', 'dog_base_overlap'])
    else:
        dog_overlap_group = dog_df.groupby(['read_id', 'dog_gene_id'])['dog_base_overlap'].sum().reset_index()
        dog_overlap_group.set_index(['read_id', 'dog_gene_id'], inplace=True)

    print("Preview of dog_overlap_group:")
    print(dog_overlap_group.head())
    
    # preview the DOG overlap 
    print("Preview of dog_df:")
    print(dog_df.head())

    # Proceed with the rest of the function
    # If dog_df is empty, the following operations would still proceed without error
    if not dog_df.empty:
        dog_overlap_group = dog_df.groupby(['read_id', 'dog_gene_id'])['dog_base_overlap'].sum().reset_index()
        dog_overlap_group.set_index(['read_id', 'dog_gene_id'], inplace=True)
    else:
        # Ensure the dataframe structure is consistent even if empty
        dog_overlap_group = pd.DataFrame(columns=['read_id', 'dog_gene_id', 'dog_base_overlap']).set_index(['read_id', 'dog_gene_id'])

    print("Preview of dog_overlap_group:")
    print(dog_overlap_group.head())

    # preview the sums - delete this later
    print("Preview of dog_overlap_group:")
    print(dog_overlap_group.head())

    # select the overlap between the alignment and the dog region of the gene to which we have the best exonic alignment 
    # First, create a new dictionary top_dog_gene which maps from read_id to dog_gene_id with highest DOG overlap
    dog_overlap_group.sort_values('dog_base_overlap', ascending=False, inplace=True)
    top_dog_overlap = dog_overlap_group.reset_index().drop_duplicates(subset='read_id', keep='first')
    top_dog_gene = {row['read_id']: row['dog_gene_id'] for _, row in top_dog_overlap.iterrows()}

    # Update the top_exon_gene dict to include genes from top_gene and top_dog_gene when there is no exonic alignment
    top_exon_gene = {row['read_id']: row['exon_gene_id'] if row['exon_base_overlap'] > 0 else top_gene.get(row['read_id'], top_dog_gene.get(row['read_id'], None)) for _, row in top_exon_overlap.iterrows()}

    # Next steps remain the same
    dog_overlap_group.reset_index(inplace=True)
    top_exon_gene_df = pd.DataFrame(list(top_exon_gene.items()), columns=['read_id', 'gene_id'])
    dog_overlap_filtered = pd.merge(dog_overlap_group, top_exon_gene_df, left_on=['read_id', 'dog_gene_id'], right_on=['read_id', 'gene_id'], how='inner')

    # Fill in NaN values with 0 in 'dog_base_overlap'
    dog_overlap_filtered['dog_base_overlap'].fillna(0, inplace=True)
    dog_overlap_filtered.reset_index(inplace=True)

    # preview of dog_overlap_sum
    print("Preview of dog_overlap_filtered:")
    print(dog_overlap_filtered.head())

    # merge dog sum into our gene_overlap_sum dataframe
    gene_overlap_sum = gene_overlap_sum.merge(dog_overlap_filtered[['read_id', 'dog_base_overlap']], how='left', left_on='read_id', right_on='read_id')
    gene_overlap_sum.rename(columns={'dog_base_overlap': 'DOG_overlap'}, inplace=True)
       
    # merge with bam_df to include all reads, regardless of overlap status
    bam_df.rename(columns={'name': 'read_id', 's-a_length': 'read-alignment_length'}, inplace=True)
        
    # take unique observations of read_id to merge into the gene_overlap_sum dataframe 
    bam_df = bam_df.drop_duplicates(subset='read_id')

    # merge in all the reads to ensure we don't lose cases with 0 exonic or genic overlap 
    # gene_overlap_sum = pd.merge(bam_df[['read_id', 'read-alignment_length']], gene_overlap_sum, how='left', on='read_id')
    gene_overlap_sum = pd.merge(bam_df[['read_id', 'read-alignment_length']], gene_overlap_sum[gene_overlap_sum.columns[~gene_overlap_sum.columns.isin(['read-alignment_length'])]], how='left', on='read_id')

    print("Preview of gene_overlap_sum after merge")
    print(gene_overlap_sum.head())

    # handle the original read-alignment-length column and split into numerical read-length and alingment-length columns 
    gene_overlap_sum[['read_length', 'alignment_length', 'splice_count']] = gene_overlap_sum['read-alignment_length'].str.split(',', expand=True)
    gene_overlap_sum['read_length'] = gene_overlap_sum['read_length'].astype(int)
    gene_overlap_sum['alignment_length'] = gene_overlap_sum['alignment_length'].astype(int)
    gene_overlap_sum['splice_count'] = gene_overlap_sum['splice_count'].astype(int)
    gene_overlap_sum.drop(columns='read-alignment_length', inplace=True)

    # add intronic alignment, unclassified length and unaligned length, and reorder the cols 
    gene_overlap_sum['intronic_alignment'] = gene_overlap_sum['gene_base_overlap'] - gene_overlap_sum['exon_base_overlap']
    gene_overlap_sum['unclassified_length'] = gene_overlap_sum['alignment_length'] - gene_overlap_sum['gene_base_overlap'] - gene_overlap_sum['DOG_overlap']
    gene_overlap_sum['unaligned_length'] = gene_overlap_sum['read_length'] - gene_overlap_sum['alignment_length']
    column_order = ['read_id', 'gene_id', 'read_length', 'alignment_length', 'splice_count', 'gene_base_overlap', 'exon_base_overlap', 'intronic_alignment', 'DOG_overlap', 'unclassified_length', 'unaligned_length']
    gene_overlap_sum = gene_overlap_sum.reindex(columns=column_order)

    # convert end coordinate data to df 
    end_coord_df = pd.DataFrame(end_coordinates.items(), columns=['read_id', 'end_coordinates'])

    # merge into dataframe 
    gene_overlap_sum = gene_overlap_sum.merge(end_coord_df, on='read_id', how='left')
    gene_overlap_sum['end_coordinates'].fillna('NA', inplace=True)
    print("Number of 'NA' values in end coordinates: ", (gene_overlap_sum['end_coordinates'] == 'NA').sum())

    # replace missing values in gene_overlap_sum
    gene_overlap_sum = gene_overlap_sum.fillna({'gene_id': np.nan, 'read_length': 0, 'alignment_length': 0, 'splice_count': 0,'gene_base_overlap': 0, 'exon_base_overlap': 0, 'intronic_alignment': 0, 'DOG_overlap': 0, 'unclassified_length': 0, 'unaligned_length': 0})

    # read in the parsed bedtools intersect output 
    gene_overlap_sum[['read_end_chromosome', 'read_end_position', 'read_end_strand']] = gene_overlap_sum['end_coordinates'].str.split(':', expand=True)
    gene_overlap_sum['read_end_position'] = gene_overlap_sum['read_end_position'].astype(int)
    gene_overlap_sum.drop(columns=['end_coordinates'], inplace=True)
    gene_overlap_sum['strand_sign'] = gene_overlap_sum['read_end_strand'].map({'+': 1, '-': -1})

    if 'read_end_chromosome' not in gene_overlap_sum.columns:
        raise ValueError("The expected 'read_end_chromosome' column is missing after split and processing.")

    return gene_overlap_sum
