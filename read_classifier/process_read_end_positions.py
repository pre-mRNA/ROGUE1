from concurrent.futures import ThreadPoolExecutor, as_completed, ProcessPoolExecutor
from functools import partial
from typing import Tuple
import pandas as pd
import numpy as np
import logging 
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

from gtf_to_bed import gtf_to_bed

def get_transcript_ends(gtf_file, output_dir): 

    # extract the annotated trascript end sites from the gtf file 
    bed_file = gtf_to_bed(gtf_file, 'transcript', output_dir)
    
    transcript_end_coordinates = pd.read_csv(
        bed_file, sep='\t', header=None,
        names=['chromosome', 'position', 'end_position', 'transcript_id', 'gene_id', 'strand'],
        dtype={'chromosome': str, 'position': 'int64', 'end_position': 'int64', 'transcript_id': str, 'gene_id': str, 'strand': str}
    )

    logging.info(f"Extracted {len(transcript_end_coordinates)} features from BED file.")

    return transcript_end_coordinates


def calculate_distance_for_unique_positions(unique_positions: pd.DataFrame, transcript_end_coordinates: pd.DataFrame, feature_name) -> pd.DataFrame:
    logging.info("ran calculate_distance_for_unique_positions")
    
    try:
        if unique_positions.empty or transcript_end_coordinates.empty:
            print("Received an empty DataFrame.")
            return pd.DataFrame()
        
        strand_sign = unique_positions['strand_sign'].iloc[0]
        
        # using .loc to avoid potential warnings
        transcript_end_coordinates_chunk = transcript_end_coordinates.loc[
            transcript_end_coordinates['strand'] == ('+' if strand_sign == 1 else '-'), :
        ].copy()
        transcript_end_coordinates_chunk.loc[:, 'read_end_chromosome'] = transcript_end_coordinates_chunk['chromosome'].astype(str)
        # transcript_end_coordinates_chunk.loc[:, 'position'] = transcript_end_coordinates_chunk['position'].astype('int64')

        # logging.info("TES chunk head")
        # print(transcript_end_coordinates_chunk.head())

        # logging.info("Data types in unique_positions:\n%s", unique_positions.dtypes)
        # logging.info("Data types in transcript_end_coordinates_chunk:\n%s", transcript_end_coordinates_chunk.dtypes)

        merged_df_downstream = pd.merge_asof(
            unique_positions.sort_values('read_end_position'), 
            transcript_end_coordinates_chunk.sort_values('position'), 
            left_on='read_end_position', 
            right_on='position', 
            by='read_end_chromosome',
            direction='forward' if strand_sign == 1 else 'backward'
        )

        # logging.info("downstream")
        # print(merged_df_downstream.head())

        merged_df_upstream = pd.merge_asof(
            unique_positions.sort_values('read_end_position'), 
            transcript_end_coordinates_chunk.sort_values('position'), 
            left_on='read_end_position', 
            right_on='position', 
            by='read_end_chromosome',
            direction='backward' if strand_sign == 1 else 'forward'
        )

        # logging.info("upstream")
        # print(merged_df_upstream.head())

        # append the feature name
        merged_df_downstream[f'{feature_name}_downstream_distance'] = (merged_df_downstream['position'] - merged_df_downstream['read_end_position']) * strand_sign
        merged_df_upstream[f'{feature_name}_upstream_distance'] = (merged_df_upstream['position'] - merged_df_upstream['read_end_position']) * strand_sign

        merged_df = pd.merge(
            merged_df_downstream[['read_end_chromosome', 'read_end_position', 'read_end_strand', f'{feature_name}_downstream_distance']],
            merged_df_upstream[['read_end_chromosome', 'read_end_position', 'read_end_strand', f'{feature_name}_upstream_distance']],
            on=['read_end_chromosome', 'read_end_position', 'read_end_strand'],
            validate="one_to_one"
        )
        return merged_df
    
    
    except Exception as e:
        print(f"An error occurred: {e}")
        return pd.DataFrame()

def calculate_distance_to_read_ends(df: pd.DataFrame, target_feature, feature_name) -> pd.DataFrame:
   
    logging.info(f"Calculating distance between read ends and {feature_name}")
    
    try:
        # create an index of unique read end positions 
        # we only calculate distances for the unique read end postitions to save some computational load 
        unique_positions = df[['read_end_chromosome', 'read_end_position', 'read_end_strand', 'strand_sign']].drop_duplicates()
        print(f"Unique positions extracted: {len(unique_positions)} records.")

        # print(unique_positions.head())
        
        results = []
        
        # in parallel, calculate the minimum distance between read end positions and annotated transcripts end sites 
        with ProcessPoolExecutor() as executor:
            future_to_position = {executor.submit(calculate_distance_for_unique_positions, group, target_feature, feature_name): group for _, group in unique_positions.groupby(['read_end_chromosome', 'strand_sign'])}
            for future in as_completed(future_to_position):
                result = future.result()
                if not result.empty:
                    results.append(result)

        merged_unique_positions = pd.concat(results, ignore_index=True)
        if merged_unique_positions.empty:
            print("No results were returned from the calculation function.")
            return pd.DataFrame()

        df_final = pd.merge(df, merged_unique_positions, on=['read_end_chromosome', 'read_end_position', 'read_end_strand'], how='left', validate="many_to_one")
        return df_final
    
    except Exception as e:
        print(f"An error occurred while calculating distances: {e}")
        return pd.DataFrame()

