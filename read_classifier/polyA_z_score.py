import logging
import numpy as np
import pandas as pd
from collections import defaultdict
import os

def classify_polya_reads(df, output_dir, co_threshold=1.6, post_threshold=1.9):
    
    # check for required columns
    required_columns = ['polya_length', 'biotype', 'alignment_length', 'gene_id', 'read_id']
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        logging.warning(f"missing required columns in dataframe: {', '.join(missing_columns)}. cannot proceed.")
        return df  # return the original DataFrame unchanged
    
    # define background biotypes and target biotype
    background_biotypes = ['snRNA', 'rRNA', 'Mt_rRNA']
    target_biotype = 'protein_coding'
    
    classification_summary = []
    gene_post_counts = defaultdict(int)
    gene_total_counts = defaultdict(int)
    df['polya_classification'] = np.nan
    df['polya_z_score'] = np.nan
    
    
    # ensure 'polya_length' is numeric
    df['polya_length'] = pd.to_numeric(df['polya_length'], errors='coerce')
    
    # define mask for valid 'polya_length' and 'alignment_length' > 200
    valid_mask = (df['polya_length'].notna()) & (df['alignment_length'] > 200)
    
    # filter for protein_coding genes
    protein_coding_mask = (df['biotype'] == target_biotype) & valid_mask
    protein_coding_df = df[protein_coding_mask].copy()
    
    logging.info("extracting background and target data")
    background_data = df[df['biotype'].isin(background_biotypes) & (df['polya_length'].notna())]['polya_length'].values
    
    # error checks
    if len(background_data) < 20:
        logging.warning(f"not enough background reads ({len(background_data)} reads). skipping classification.")
        return df
    if len(protein_coding_df) < 20:
        logging.warning(f"not enough protein_coding reads ({len(protein_coding_df)} reads). skipping classification.")
        return df
    
    # calculate mean and std of background_data
    background_mean = np.mean(background_data)
    background_std = np.std(background_data)
    if background_std == 0:
        logging.warning("standard deviation of background data is zero. cannot compute z-scores. skipping classification.")
        return df
    
    logging.info("calculating z-scores for protein_coding reads")
    # calculate Z-scores
    protein_coding_df['polya_z_score'] = (protein_coding_df['polya_length'] - background_mean) / background_std
    
    # classify from thresholds
    def classify_z_score(z):
        if z < co_threshold:
            return 'co_transcriptional'
        elif z > post_threshold:
            return 'post_transcriptional'
        else:
            return 'uncertain'
    
    protein_coding_df['polya_classification'] = protein_coding_df['polya_z_score'].apply(classify_z_score)
    
    # update the main dataframe with classification results
    df.loc[protein_coding_df.index, 'polya_z_score'] = protein_coding_df['polya_z_score']
    df.loc[protein_coding_df.index, 'polya_classification'] = protein_coding_df['polya_classification']
    
    # summary stats 
    classification_counts = protein_coding_df['polya_classification'].value_counts()
    co_transcriptional_count = classification_counts.get('co_transcriptional', 0)
    post_transcriptional_count = classification_counts.get('post_transcriptional', 0)
    uncertain_count = classification_counts.get('uncertain', 0)
    total_classified = co_transcriptional_count + post_transcriptional_count
    co_percent = (co_transcriptional_count / total_classified) * 100 if total_classified > 0 else 0
    post_percent = (post_transcriptional_count / total_classified) * 100 if total_classified > 0 else 0
    uncertain_percent = (uncertain_count / len(protein_coding_df)) * 100 if len(protein_coding_df) > 0 else 0
    
    classification_summary.append({
        'total_protein_coding_reads': len(protein_coding_df),
        'co_transcriptional_reads': co_transcriptional_count,
        'post_transcriptional_reads': post_transcriptional_count,
        'uncertain_reads': uncertain_count,
        'co_percent': round(co_percent, 2),
        'post_percent': round(post_percent, 2),
        'uncertain_percent': round(uncertain_percent, 2)
    })
    
    # save the summary data 
    classification_summary_df = pd.DataFrame(classification_summary)
    summary_csv_path = os.path.join(output_dir, 'classification_summary.csv')
    classification_summary_df.to_csv(summary_csv_path, index=False)
    logging.info(f"classification summary saved at {summary_csv_path}")
    
    total_post_trans = classification_summary_df['post_transcriptional_reads'].iloc[0]
    total_protein_coding = classification_summary_df['total_protein_coding_reads'].iloc[0]
    overall_post_percent = (total_post_trans / total_protein_coding) * 100 if total_protein_coding > 0 else 0
    logging.info(f"overall percentage of post-transcriptional reads: {overall_post_percent:.2f}%")
    
    gene_post_counts = df[df['polya_classification'] == 'post_transcriptional']['gene_id'].value_counts().to_dict()
    gene_total_counts = df[df['biotype'] == 'protein_coding']['gene_id'].value_counts().to_dict()
    
    gene_stats = pd.DataFrame({
        'gene_id': list(gene_total_counts.keys()),
        'total_protein_coding_reads': list(gene_total_counts.values()),
        'post_transcriptional_reads': [gene_post_counts.get(gene_id, 0) for gene_id in gene_total_counts.keys()]
    })
    
    # calculate percentage of post_transcriptional reads per gene_id
    gene_stats['post_percent'] = (gene_stats['post_transcriptional_reads'] / gene_stats['total_protein_coding_reads']) * 100
    gene_stats['post_percent'] = gene_stats['post_percent'].round(2)
    
    # calculate the median percentage
    median_post_percent = gene_stats['post_percent'].median()
    logging.info(f"median percentage of post-transcriptional reads per gene_id: {median_post_percent:.2f}%")
    
    # save gene_stats 
    gene_stats_csv_path = os.path.join(output_dir, 'gene_id_post_transcriptional_stats.csv')
    gene_stats.to_csv(gene_stats_csv_path, index=False)
    logging.info(f"gene-wise post-transcriptional statistics saved at {gene_stats_csv_path}")
    
    return df
