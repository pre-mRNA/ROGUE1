import pandas as pd
import numpy as np

def get_highest_intron(df):
    intron_numbers = df['intron_ids'].fillna('').str.split(',').explode()
    intron_numbers = pd.to_numeric(intron_numbers, errors='coerce')
    return int(intron_numbers.max()) if not intron_numbers.empty and not pd.isna(intron_numbers.max()) else 0

def reclassify_dog_regions(group):
    highest_exon = get_highest_intron(group)
    
    def reclassify(row):
        if pd.isna(row['three_prime_feature']):
            return 'unclassified'
        if row['three_prime_feature'] == 'region_DOG' and row['dog_base_overlap'] < 20 and row['exon_base_overlap'] > 20:
            return f"exon_{highest_exon}"
        return row['three_prime_feature']
    
    group['three_prime_feature'] = group.apply(reclassify, axis=1)
    return group

def classify_read(row):
    exon_ids = [int(x) for x in str(row['exon_ids']).split(',') if x.strip()] if pd.notna(row['exon_ids']) else []
    intron_ids = [int(x) for x in str(row['intron_ids']).split(',') if x.strip()] if pd.notna(row['intron_ids']) else []
    
    if row['splice_count'] > 0:
        if len(exon_ids) != row['splice_count'] + 1:
            return 'cryptic_splicing'
    
    if len(exon_ids) + len(intron_ids) < 2:
        return 'insufficient_features'
    
    if len(intron_ids) == 1 and (not exon_ids or (len(exon_ids) == 1 and exon_ids[0] == intron_ids[0])):
        return 'ambiguous'
    
    if row['splice_count'] > 0 and intron_ids:
        if len(exon_ids) > 1 and max(exon_ids) > min(intron_ids):
            return 'partially_spliced'
        
    if row['splice_count'] == 0:
        if len(intron_ids) >= 1 and len(exon_ids) >= 1:
            return 'fully_unspliced'
    
    if row['intron_base_overlap'] > 75 and row['splice_count'] > 0:
        return 'partially_spliced'
    elif row['intron_base_overlap'] < 25 and row['splice_count'] >= 1:
        return 'spliced'
    else:
        if row['intron_base_overlap'] <= 75 and row['intron_base_overlap'] >= 25:
            return 'ambiguous'
        elif row['splice_count'] == 0 and row['intron_base_overlap'] < 75:
            return 'no_splice_low_intron'
        else:
            return 'other'

def classify_splicing(df):
    
    # handle empty 'three_prime_feature' values
    df['three_prime_feature'] = df['three_prime_feature'].fillna('unclassified')
    
    # identify the highest intron for each gene
    highest_introns = df.groupby('gene_id').apply(get_highest_intron).reset_index()
    highest_introns.columns = ['gene_id', 'highest_intron']
    df = pd.merge(df, highest_introns, on='gene_id', how='left')
    df['highest_intron'] = df['highest_intron'].fillna(0).astype(int)
    
    # reclassify DOG regions
    df = df.groupby('gene_id', group_keys=False).apply(reclassify_dog_regions)
    
    # classify splicing status 
    df['canonical_acceptor_splicing_status'] = df.apply(classify_read, axis=1)
    
    return df