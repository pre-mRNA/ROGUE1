import os
import warnings
import pandas as pd
import numpy as np
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import multiprocessing
from functools import reduce
import logging 

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

def process_overlap_group(df, group_cols, value_col):
    return df.groupby(group_cols)[value_col].sum().reset_index()

def filter_and_sum_overlaps(overlap_df, type_name, overlap_threshold):
    filtered_df = overlap_df[overlap_df[f'{type_name}_base_overlap'] >= overlap_threshold]
    return filtered_df.groupby(['read_id', f'{type_name}_id', f'{type_name}_gene_id'])[f'{type_name}_base_overlap'].sum().reset_index()

def process_chunk(chunk_data, overlap_threshold):
    exon_df, intron_df, gene_df, dog_df = chunk_data

    def get_most_upstream_feature(read_data, strand):
        if strand == '+':
            return min(read_data, key=lambda x: x['start'])
        else:
            return max(read_data, key=lambda x: x['end'])

    def assign_gene(read_data):
        if not read_data:
            return None

        strand = read_data[0]['strand']
        most_upstream = get_most_upstream_feature(read_data, strand)

        if most_upstream['type'] == 'exon':
            return most_upstream['gene_id']

        exon_intron_features = [f for f in read_data if f['type'] in ['exon', 'intron']]
        if exon_intron_features:
            gene_total_overlap = defaultdict(int)
            for f in exon_intron_features:
                gene_total_overlap[f['gene_id']] += f['overlap']
            return max(gene_total_overlap, key=gene_total_overlap.get)

        dog_features = [f for f in read_data if f['type'] == 'dog']
        if dog_features:
            return max(dog_features, key=lambda x: x['overlap'])['gene_id']

        return None

    read_data = defaultdict(list)
    for df, feature_type in [(exon_df, 'exon'), (intron_df, 'intron'), (gene_df, 'gene'), (dog_df, 'dog')]:
        for _, row in df.iterrows():
            gene_id_col = 'gene_id' if feature_type == 'gene' else f'{feature_type}_gene_id'
            read_data[row['read_id']].append({
                'type': feature_type,
                'gene_id': row[gene_id_col],
                'start': row[f'{feature_type}_start'],
                'end': row[f'{feature_type}_end'],
                'strand': row[f'{feature_type}_strand'],
                'overlap': row[f'{feature_type}_base_overlap']
            })

    top_gene = {read_id: assign_gene(data) for read_id, data in read_data.items()}
    top_gene = {k: v for k, v in top_gene.items() if v is not None}

    # initialise overlap groups 
    exon_overlap_group = pd.DataFrame(columns=['read_id', 'exon_gene_id', 'exon_base_overlap'])
    intron_overlap_group = pd.DataFrame(columns=['read_id', 'intron_gene_id', 'intron_base_overlap'])
    gene_overlap_group = pd.DataFrame(columns=['read_id', 'gene_id', 'gene_base_overlap'])
    dog_overlap_group = pd.DataFrame(columns=['read_id', 'dog_gene_id', 'dog_base_overlap'])

    if not exon_df.empty:
        exon_overlap_group = process_overlap_group(exon_df[exon_df.apply(lambda row: top_gene.get(row['read_id']) == row['exon_gene_id'], axis=1)], ['read_id', 'exon_gene_id'], 'exon_base_overlap')

    if not intron_df.empty:
        intron_overlap_group = process_overlap_group(intron_df[intron_df.apply(lambda row: top_gene.get(row['read_id']) == row['intron_gene_id'], axis=1)], ['read_id', 'intron_gene_id'], 'intron_base_overlap')
    
    if not gene_df.empty:
        gene_overlap_group = process_overlap_group(gene_df[gene_df.apply(lambda row: top_gene.get(row['read_id']) == row['gene_id'], axis=1)], ['read_id', 'gene_id'], 'gene_base_overlap')
    
    if not dog_df.empty:
        dog_overlap_group = process_overlap_group(dog_df[dog_df.apply(lambda row: top_gene.get(row['read_id']) == row['dog_gene_id'], axis=1)], ['read_id', 'dog_gene_id'], 'dog_base_overlap')

    exon_filtered_summed = filter_and_sum_overlaps(exon_df, 'exon', overlap_threshold)
    intron_filtered_summed = filter_and_sum_overlaps(intron_df, 'intron', overlap_threshold)

    all_overlap_groups = [exon_overlap_group, intron_overlap_group, gene_overlap_group, dog_overlap_group]
    gene_overlap_sum = reduce(lambda left, right: pd.merge(left, right, on='read_id', how='outer'), all_overlap_groups)

    overlap_columns = ['exon_base_overlap', 'intron_base_overlap', 'gene_base_overlap', 'dog_base_overlap']
    gene_overlap_sum.loc[:, overlap_columns] = gene_overlap_sum.loc[:, overlap_columns].infer_objects(copy = False)

    gene_overlap_sum['gene_id'] = gene_overlap_sum.apply(lambda row: top_gene.get(row['read_id'], np.nan), axis=1)

    return gene_overlap_sum, exon_filtered_summed, intron_filtered_summed

def parse_output(exon_overlap_file, intron_overlap_file, gene_overlap_file, dog_overlap_file, bed_file, end_coordinates, record_exons, genome_file, output_dir, num_files, original_exon_bed, original_intron_bed, original_dog_bed, overlap_threshold=7):
    
    missing_files = []
    empty_files = []
    for file in [exon_overlap_file, intron_overlap_file, gene_overlap_file, dog_overlap_file]:
        if not os.path.isfile(file):
            missing_files.append(file)
        elif os.stat(file).st_size == 0:
            empty_files.append(file)

    if missing_files:
        raise Exception(f"Essential file(s) do not exist: {', '.join(missing_files)}. Unable to proceed.")
    if empty_files:
        logging.warning(f"The following file(s) are empty and will be skipped: {', '.join(empty_files)}")

    exon_cols = ['read_chrom', 'read_fragment_start', 'read_fragment_end', 'read_id', 'read-alignment_length', 'read_strand', 'exon_chrom', 'exon_start', 'exon_end', 'exon_gene_id', 'exon_id', 'exon_strand', 'exon_base_overlap']
    intron_cols = ['read_chrom', 'read_fragment_start', 'read_fragment_end', 'read_id', 'read-alignment_length', 'read_strand', 'intron_chrom', 'intron_start', 'intron_end', 'intron_gene_id', 'intron_id', 'intron_strand', 'intron_base_overlap']
    gene_cols = ['read_chrom', 'read_fragment_start', 'read_fragment_end', 'read_id', 'read-alignment_length', 'read_strand', 'gene_chrom', 'gene_start', 'gene_end', 'gene_id', 'gene_biotype', 'gene_strand', 'gene_base_overlap']
    dog_cols = ['read_chrom', 'read_fragment_start', 'read_fragment_end', 'read_id', 'read-alignment_length', 'read_strand', 'dog_chrom', 'dog_start', 'dog_end', 'dog_gene_id', 'dog_id', 'dog_strand', 'dog_base_overlap']
    bam_to_bed_cols = ['chrom', 'start', 'end', 'name', 's-a_length', 'strand']


    # read in all alignments 
    bam_df = pd.read_csv(bed_file, sep="\t", header=None, names=bam_to_bed_cols, low_memory=False)
    bam_df.rename(columns={'name': 'read_id', 's-a_length': 'read-alignment_length'}, inplace=True)
    
    # read in the overlaps 
    exon_df = pd.read_csv(exon_overlap_file, sep="\t", header=None, names=exon_cols, low_memory=False)
    intron_df = pd.read_csv(intron_overlap_file, sep="\t", header=None, names=intron_cols, low_memory=False)
    gene_df = pd.read_csv(gene_overlap_file, sep="\t", header=None, names=gene_cols, low_memory=False)
    
    if os.stat(dog_overlap_file).st_size == 0:
        print("Warning: DOG overlap file is empty. Proceeding without DOG information.")
        dog_df = pd.DataFrame(columns=dog_cols)
    else:
        dog_df = pd.read_csv(dog_overlap_file, sep="\t", header=None, names=dog_cols, low_memory=False)

    bam_df = pd.read_csv(bed_file, sep="\t", header=None, names=bam_to_bed_cols, low_memory=False)

    # unqiue read IDs across all overlap files
    all_read_ids = set(exon_df['read_id']) | set(intron_df['read_id']) | set(gene_df['read_id']) | set(dog_df['read_id'])

    # distribute read IDs across CPUs 
    num_cpus = multiprocessing.cpu_count()
    read_id_chunks = np.array_split(list(all_read_ids), num_cpus)

    def create_chunk_data(read_ids):
        return (
            exon_df[exon_df['read_id'].isin(read_ids)],
            intron_df[intron_df['read_id'].isin(read_ids)],
            gene_df[gene_df['read_id'].isin(read_ids)],
            dog_df[dog_df['read_id'].isin(read_ids)]
        )

    # parallel chunking stratetgy
    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        chunk_results = list(executor.map(process_chunk, 
                                          [create_chunk_data(chunk) for chunk in read_id_chunks],
                                          [overlap_threshold] * len(read_id_chunks)))

    gene_overlap_sum = pd.concat([result[0] for result in chunk_results])
    exon_filtered_summed = pd.concat([result[1] for result in chunk_results])
    intron_filtered_summed = pd.concat([result[2] for result in chunk_results])

    bam_df.rename(columns={'name': 'read_id', 's-a_length': 'read-alignment_length'}, inplace=True)
    bam_df = bam_df.drop_duplicates(subset='read_id')
    gene_overlap_sum = pd.merge(bam_df[['read_id', 'read-alignment_length']], gene_overlap_sum, how='left', on='read_id')


    # read lengths, alignment lengths, splice counts
    gene_overlap_sum[['read_length', 'alignment_length', 'splice_count']] = gene_overlap_sum['read-alignment_length'].str.split(',', expand=True)
    gene_overlap_sum['read_length'] = pd.to_numeric(gene_overlap_sum['read_length'], errors='coerce').fillna(0).astype(int)
    gene_overlap_sum['alignment_length'] = pd.to_numeric(gene_overlap_sum['alignment_length'], errors='coerce').fillna(0).astype(int)
    gene_overlap_sum['splice_count'] = pd.to_numeric(gene_overlap_sum['splice_count'], errors='coerce').fillna(0).astype(int)
    gene_overlap_sum.drop(columns='read-alignment_length', inplace=True)

    # unclassified and unaligned lengths
    gene_overlap_sum['unclassified_length'] = gene_overlap_sum['alignment_length'] - gene_overlap_sum['gene_base_overlap'] - gene_overlap_sum['dog_base_overlap']
    gene_overlap_sum['unaligned_length'] = gene_overlap_sum['read_length'] - gene_overlap_sum['alignment_length']

    end_coord_df = pd.DataFrame(end_coordinates.items(), columns=['read_id', 'end_coordinates'])
    gene_overlap_sum = gene_overlap_sum.merge(end_coord_df, on='read_id', how='left')
    gene_overlap_sum['end_coordinates'].fillna('NA', inplace=True)

    gene_overlap_sum = gene_overlap_sum.fillna({
        'gene_id': np.nan, 'read_length': 0, 'alignment_length': 0, 'splice_count': 0,
        'gene_base_overlap': 0, 'exon_base_overlap': 0, 'intron_base_overlap': 0,
        'dog_base_overlap': 0, 'unclassified_length': 0, 'unaligned_length': 0
    })

    gene_overlap_sum[['read_end_chromosome', 'read_end_position', 'read_end_strand']] = gene_overlap_sum['end_coordinates'].str.split(':', expand=True)
    gene_overlap_sum['read_end_position'] = pd.to_numeric(gene_overlap_sum['read_end_position'], errors='coerce').fillna(0).astype(int)
    gene_overlap_sum.drop(columns=['end_coordinates'], inplace=True)
    gene_overlap_sum['strand_sign'] = gene_overlap_sum['read_end_strand'].map({'+': 1, '-': -1})

    if record_exons:
        read_exon_ids = defaultdict(lambda: defaultdict(set))
        read_intron_ids = defaultdict(lambda: defaultdict(set))

        for _, row in exon_filtered_summed.iterrows():
            read_exon_ids[row['read_id']][row['exon_gene_id']].add(row['exon_id'])

        for _, row in intron_filtered_summed.iterrows():
            read_intron_ids[row['read_id']][row['intron_gene_id']].add(row['intron_id'])

        def collapse_ids(id_set):
            return ','.join(sorted(set(id.split('_')[-1] for id in id_set)))

        gene_overlap_sum['exon_ids'] = gene_overlap_sum.apply(lambda row: collapse_ids(read_exon_ids.get(row['read_id'], {}).get(row['gene_id'], [])), axis=1)
        gene_overlap_sum['intron_ids'] = gene_overlap_sum.apply(lambda row: collapse_ids(read_intron_ids.get(row['read_id'], {}).get(row['gene_id'], [])), axis=1)

    column_order = ['read_id', 'gene_id', 'read_length', 'alignment_length', 'splice_count',
                    'gene_base_overlap', 'exon_base_overlap', 'intron_base_overlap', 'dog_base_overlap',
                    'unclassified_length', 'unaligned_length', 
                    'read_end_chromosome', 'read_end_position', 'read_end_strand', 'strand_sign']
    
    if record_exons:
            column_order.extend(['exon_ids', 'intron_ids'])

    gene_overlap_sum = gene_overlap_sum.reindex(columns=column_order)

    if 'read_end_chromosome' not in gene_overlap_sum.columns:
        raise ValueError("The expected 'read_end_chromosome' column is missing after split and processing.")

    print("Number of 'NA' values in end coordinates: ", (gene_overlap_sum['read_end_chromosome'] == 'NA').sum())

    missing_columns = set(column_order) - set(gene_overlap_sum.columns)
    if missing_columns:
        raise ValueError(f"The following expected columns are missing: {', '.join(missing_columns)}")

    # summarise
    print(f"Total number of reads processed: {len(gene_overlap_sum)}")
    print(f"Number of unique genes: {gene_overlap_sum['gene_id'].nunique()}")

    return gene_overlap_sum