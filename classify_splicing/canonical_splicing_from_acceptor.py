import pandas as pd
import numpy as np

def get_highest_intron(intron_ids):
    """
    Computes the highest intron number from a comma-separated string of intron IDs.
    Returns 0 if no valid intron numbers are found.
    """
    if pd.isna(intron_ids) or intron_ids.strip() == '':
        return 0
    
    # split the intron_ids by comma and convert to numeric, ignoring errors
    intron_numbers = pd.to_numeric(intron_ids.split(','), errors='coerce')
    highest = intron_numbers.max()
    return int(highest) if not pd.isna(highest) else 0

def reclassify_dog_regions(df):
    """
    Reclassifies the 'three_prime_feature' for DOG regions based on specific overlap criteria.
    Utilizes the 'highest_intron' value for each row.
    """
    
    # classification condition
    condition = (
        (df['three_prime_feature'] == 'region_DOG') &
        (df['dog_base_overlap'] < 20) &
        (df['exon_base_overlap'] > 20)
    )
    
    # classify 
    df.loc[condition, 'three_prime_feature'] = 'exon_' + (df.loc[condition, 'highest_intron'] + 1).astype(str)

    df['three_prime_feature'] = df['three_prime_feature'].fillna('unclassified')
    
    return df

def classify_read(row):
    """
    Classifies the splicing status of a read based on exon and intron overlaps and splice counts.
    """
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
        if 25 <= row['intron_base_overlap'] <= 75:
            return 'ambiguous'
        elif row['splice_count'] == 0 and row['intron_base_overlap'] < 75:
            return 'no_splice_low_intron'
        else:
            return 'other'

def classify_splicing(df):
    """
    Classifies splicing status for each read in the DataFrame without losing any rows.
    """

    df['three_prime_feature'] = df['three_prime_feature'].fillna('unclassified')
    
    # get the highest observed intron for each read 
    df['highest_intron'] = df['intron_ids'].apply(get_highest_intron)
    
    # calculate the longest observed intron for each gene 
    df['highest_intron'] = df.groupby('gene_id')['highest_intron'].transform('max').fillna(0).astype(int)
    
    # reclassify DOG regions based on the computed 'highest_intron'
    df = reclassify_dog_regions(df)
    
    # classify splicing status for each read
    df['canonical_acceptor_splicing_status'] = df.apply(classify_read, axis=1)
    
    return df
