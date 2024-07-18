import os
import warnings
import pandas as pd 
import numpy as np 
from collections import defaultdict

def parse_output(exon_overlap_file, intron_overlap_file, gene_overlap_file, dog_overlap_file, bed_file, end_coordinates, record_exons):
    for file in [exon_overlap_file, intron_overlap_file, gene_overlap_file, dog_overlap_file]:
        if not os.path.isfile(file):
            raise Exception(f"Essential file {file} does not exist. Unable to proceed.")
        elif os.stat(file).st_size == 0:
            warnings.warn(f"Warning: Essential file {file} is empty. Some data may not be processed.")
            
    exon_cols = ['read_chrom', 'read_fragment_start', 'read_fragment_end', 'read_id', 'read-alignment_length', 'read_strand', 'exon_chrom', 'exon_start', 'exon_end', 'exon_gene_id', 'exon_id', 'exon_strand', 'exon_base_overlap']
    intron_cols = ['read_chrom', 'read_fragment_start', 'read_fragment_end', 'read_id', 'read-alignment_length', 'read_strand', 'intron_chrom', 'intron_start', 'intron_end', 'intron_gene_id', 'intron_id', 'intron_strand', 'intron_base_overlap']
    gene_cols = ['read_chrom', 'read_fragment_start', 'read_fragment_end', 'read_id', 'read-alignment_length', 'read_strand', 'gene_chrom', 'gene_start', 'gene_end', 'gene_id', 'gene_biotype', 'gene_strand', 'gene_base_overlap']
    bam_to_bed_cols = ['chrom', 'start', 'end', 'name', 's-a_length', 'strand']

    exon_df = pd.read_csv(exon_overlap_file, sep="\t", header=None, names=exon_cols, low_memory=False)
    intron_df = pd.read_csv(intron_overlap_file, sep="\t", header=None, names=intron_cols, low_memory=False)
    gene_df = pd.read_csv(gene_overlap_file, sep="\t", header=None, names=gene_cols, low_memory=False)
    bam_df = pd.read_csv(bed_file, sep="\t", header=None, names=bam_to_bed_cols, low_memory=False)

    # select the 5' most coordinate for each overlap
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

    # filter for the selected exon overlaps
    selected_exon_df = exon_df[exon_df.apply(lambda row: top_exon_gene[row['read_id']] == row['exon_gene_id'], axis=1)]

    # Cclculate total overlap between each read with its corresponding gene
    exon_overlap_group = selected_exon_df.groupby(['read_id', 'exon_gene_id'])['exon_base_overlap'].sum().reset_index()

    # calculate total overlap between each read with its corresponding intron
    intron_overlap_group = intron_df.groupby(['read_id', 'intron_gene_id'])['intron_base_overlap'].sum().reset_index()

    # identify the gene to which each read has the most alignment
    gene_overlap_group = gene_df.groupby(['read_id', 'gene_id'])['gene_base_overlap'].sum().reset_index()
    gene_overlap_group.sort_values('gene_base_overlap', ascending=False, inplace=True)
    top_gene_overlap = gene_overlap_group.drop_duplicates(subset='read_id', keep='first')
    
    # create a dict to map from read_id to gene_id with highest gene overlap
    top_gene = {row['read_id']: row['gene_id'] for _, row in top_gene_overlap.iterrows()}

    # update dict with exon information 
    for read_id, exon_gene_id in top_exon_gene.items():
        top_gene[read_id] = exon_gene_id

    gene_df_filtered = gene_df[gene_df.apply(lambda row: top_gene.get(row['read_id']) == row['gene_id'], axis=1)]

    # get total overlap with target gene 
    gene_overlap_sum = gene_df_filtered.groupby(['read_id', 'gene_id', 'read-alignment_length'])['gene_base_overlap'].sum().reset_index()

    # merge
    gene_overlap_sum = pd.merge(gene_overlap_sum, exon_overlap_group, left_on=['read_id', 'gene_id'], right_on=['read_id', 'exon_gene_id'], how='left')
    gene_overlap_sum = pd.merge(gene_overlap_sum, intron_overlap_group, left_on=['read_id', 'gene_id'], right_on=['read_id', 'intron_gene_id'], how='left')
    
    gene_overlap_sum['exon_base_overlap'].fillna(0, inplace=True)
    gene_overlap_sum['intron_base_overlap'].fillna(0, inplace=True)

    # Rename the column for gene - exon overlap
    gene_overlap_sum['gene_exon_bases'] = gene_overlap_sum['gene_base_overlap'] - gene_overlap_sum['exon_base_overlap']

    # handle DOGs
    dog_cols = ['read_chrom', 'read_fragment_start', 'read_fragment_end', 'read_id', 'read-alignment_length', 'read_strand', 'gene_chrom', 'gene_start', 'gene_end', 'dog_gene_id', 'gene_biotype', 'gene_strand', 'dog_base_overlap']
    dog_df = pd.read_csv(dog_overlap_file, sep="\t", header=None, names=dog_cols, low_memory=False)

    if dog_df.empty:
        dog_overlap_group = pd.DataFrame(columns=['read_id', 'dog_gene_id', 'dog_base_overlap'])
    else:
        dog_overlap_group = dog_df.groupby(['read_id', 'dog_gene_id'])['dog_base_overlap'].sum().reset_index()

    # select the overlap between the alignment and the dog region of the gene to which we have the best exonic alignment 
    dog_overlap_group.sort_values('dog_base_overlap', ascending=False, inplace=True)
    top_dog_overlap = dog_overlap_group.drop_duplicates(subset='read_id', keep='first')
    top_dog_gene = {row['read_id']: row['dog_gene_id'] for _, row in top_dog_overlap.iterrows()}

    # update the top_exon_gene dict to include genes from top_gene and top_dog_gene when there is no exonic alignment
    top_exon_gene = {row['read_id']: row['exon_gene_id'] if row['exon_base_overlap'] > 0 else top_gene.get(row['read_id'], top_dog_gene.get(row['read_id'], None)) for _, row in top_exon_overlap.iterrows()}

    top_exon_gene_df = pd.DataFrame(list(top_exon_gene.items()), columns=['read_id', 'gene_id'])
    dog_overlap_filtered = pd.merge(dog_overlap_group, top_exon_gene_df, left_on=['read_id', 'dog_gene_id'], right_on=['read_id', 'gene_id'], how='inner')

    dog_overlap_filtered['dog_base_overlap'].fillna(0, inplace=True)

    # merge dog sum 
    gene_overlap_sum = gene_overlap_sum.merge(dog_overlap_filtered[['read_id', 'dog_base_overlap']], how='left', left_on='read_id', right_on='read_id')
    gene_overlap_sum.rename(columns={'dog_base_overlap': 'DOG_overlap'}, inplace=True)

    # merge with bam_df to include all reads, regardless of overlap status
    bam_df.rename(columns={'name': 'read_id', 's-a_length': 'read-alignment_length'}, inplace=True)
    bam_df = bam_df.drop_duplicates(subset='read_id')
    gene_overlap_sum = pd.merge(bam_df[['read_id', 'read-alignment_length']], gene_overlap_sum[gene_overlap_sum.columns[~gene_overlap_sum.columns.isin(['read-alignment_length'])]], how='left', on='read_id')

    # handle the original read-alignment-length column; split into read-length and alignment-length cols
    gene_overlap_sum[['read_length', 'alignment_length', 'splice_count']] = gene_overlap_sum['read-alignment_length'].str.split(',', expand=True)
    gene_overlap_sum['read_length'] = gene_overlap_sum['read_length'].astype(int)
    gene_overlap_sum['alignment_length'] = gene_overlap_sum['alignment_length'].astype(int)
    gene_overlap_sum['splice_count'] = gene_overlap_sum['splice_count'].astype(int)
    gene_overlap_sum.drop(columns='read-alignment_length', inplace=True)

    # consider unclassified and unalgined length 
    gene_overlap_sum['unclassified_length'] = gene_overlap_sum['alignment_length'] - gene_overlap_sum['gene_base_overlap'] - gene_overlap_sum['DOG_overlap']
    gene_overlap_sum['unaligned_length'] = gene_overlap_sum['read_length'] - gene_overlap_sum['alignment_length']

    # handle read end coordinate
    end_coord_df = pd.DataFrame(end_coordinates.items(), columns=['read_id', 'end_coordinates'])
    gene_overlap_sum = gene_overlap_sum.merge(end_coord_df, on='read_id', how='left')
    gene_overlap_sum['end_coordinates'].fillna('NA', inplace=True)

    gene_overlap_sum = gene_overlap_sum.fillna({
        'gene_id': np.nan, 'read_length': 0, 'alignment_length': 0, 'splice_count': 0,
        'gene_base_overlap': 0, 'exon_base_overlap': 0, 'intron_base_overlap': 0,
        'gene_exon_bases': 0, 'DOG_overlap': 0, 'unclassified_length': 0, 'unaligned_length': 0
    })

    # read in each parse bedtools intersect output 
    gene_overlap_sum[['read_end_chromosome', 'read_end_position', 'read_end_strand']] = gene_overlap_sum['end_coordinates'].str.split(':', expand=True)
    gene_overlap_sum['read_end_position'] = gene_overlap_sum['read_end_position'].astype(int)
    gene_overlap_sum.drop(columns=['end_coordinates'], inplace=True)
    gene_overlap_sum['strand_sign'] = gene_overlap_sum['read_end_strand'].map({'+': 1, '-': -1})

    # if record_exons, record exon and intron IDs with overlaps in the target gene
    # NOTE: we could consider adding a threshold to introns, since noisy reads 
    # can leak into the intron and cause spurious classification of intron spanning reads 
    if record_exons:
        read_exon_ids = defaultdict(lambda: defaultdict(set))
        read_intron_ids = defaultdict(lambda: defaultdict(set))
        
        for _, row in exon_df.iterrows():
            read_exon_ids[row['read_id']][row['exon_gene_id']].add(row['exon_id'])
        
        for _, row in intron_df.iterrows():
            read_intron_ids[row['read_id']][row['intron_gene_id']].add(row['intron_id'])

        def collapse_ids(id_set):
            return ','.join(sorted(set(id.split('_')[-1] for id in id_set)))

        gene_overlap_sum['exon_ids'] = gene_overlap_sum.apply(lambda row: collapse_ids(read_exon_ids.get(row['read_id'], {}).get(row['gene_id'], [])), axis=1)
        gene_overlap_sum['intron_ids'] = gene_overlap_sum.apply(lambda row: collapse_ids(read_intron_ids.get(row['read_id'], {}).get(row['gene_id'], [])), axis=1)

        # check for gene_id consistency in the exon and intron overlaps 
        def check_gene_id_consistency(row):
            gene_id = row['gene_id']
            exon_gene_id = row['exon_gene_id'] if pd.notna(row['exon_gene_id']) else gene_id
            intron_gene_id = row['intron_gene_id'] if pd.notna(row['intron_gene_id']) else gene_id
            if not (gene_id == exon_gene_id == intron_gene_id):
                raise ValueError(f"Inconsistent gene_ids for read {row['read_id']}: gene={gene_id}, exon={exon_gene_id}, intron={intron_gene_id}")
            return True

        gene_overlap_sum.apply(check_gene_id_consistency, axis=1)

    # update col order 
    column_order = ['read_id', 'gene_id', 'read_length', 'alignment_length', 'splice_count', 
                    'gene_base_overlap', 'exon_base_overlap', 'intron_base_overlap', 'gene_exon_bases', 'DOG_overlap', 
                    'unclassified_length', 'unaligned_length', 'read_end_chromosome', 
                    'read_end_position', 'read_end_strand', 'strand_sign']
    if record_exons:
        column_order.extend(['exon_ids', 'intron_ids'])

    gene_overlap_sum = gene_overlap_sum.reindex(columns=column_order)

    # final check 
    if 'read_end_chromosome' not in gene_overlap_sum.columns:
        raise ValueError("The expected 'read_end_chromosome' column is missing after split and processing.")

    # summarise 
    print("Number of 'NA' values in end coordinates: ", (gene_overlap_sum['read_end_chromosome'] == 'NA').sum())

    # check columns 
    missing_columns = set(column_order) - set(gene_overlap_sum.columns)
    if missing_columns:
        raise ValueError(f"The following expected columns are missing: {', '.join(missing_columns)}")

    return gene_overlap_sum