import numpy as np
from scipy import stats
import logging

def calculate_polya_zscores(df, rRNA_cutoff = 1000):
    """
    Calculate polyA tail length  z-scores for all reads, against a reference biotype (rRNA). 
    :param df: DataFrame containing read information including 'biotype' and 'polya_length'
    :return: DataFrame with additional 'polya_zscore' column, or the original DataFrame if calculation is skipped
    """

    # filter rRNA reads and drop NaN values
    rRNA_lengths = df[(df['biotype'] == 'rRNA') & (df['polya_length'].notna())]['polya_length']
    
    # check if we have at least rRNA_cutoff rRNA reads
    if len(rRNA_lengths) < rRNA_cutoff:
        logging.info(f"Skipping polyA z-score calculation: Less than {rRNA_cutoff} rRNA reads available.")
        return df
    
    # calculate mean and standard deviation of rRNA polyA lengths
    rRNA_mean = rRNA_lengths.mean()
    rRNA_std = rRNA_lengths.std()
    
    logging.info(f"rRNA polyA length stats: mean = {rRNA_mean:.2f}, std = {rRNA_std:.2f}")
    
    def calculate_zscore(x):
        return (x - rRNA_mean) / rRNA_std
    
    # calculate z-scores for all reads
    df['polya_zscore'] = df['polya_length'].apply(calculate_zscore)
    
    # z-score logging 
    for biotype in df['biotype'].unique():
        biotype_zscores = df[df['biotype'] == biotype]['polya_zscore'].dropna()
        if len(biotype_zscores) > 0:
            mean_zscore = biotype_zscores.mean()
            median_zscore = biotype_zscores.median()
            prop_exceeding = (biotype_zscores > 2.33).mean()
            logging.info(f"Biotype: {biotype}")
            logging.info(f"  Number of reads: {len(biotype_zscores)}")
            logging.info(f"  Mean z-score: {mean_zscore:.2f}")
            logging.info(f"  Median z-score: {median_zscore:.2f}")
            logging.info(f"  Proportion of reads exceeding z-score of 2.33: {prop_exceeding:.2%}")
    
    return df

