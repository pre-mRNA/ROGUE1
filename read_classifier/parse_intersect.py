import os
import warnings
import pandas as pd 
import numpy as np 
from collections import defaultdict
import tempfile
import subprocess
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
import multiprocessing
from functools import reduce

def process_overlap_group(df, group_cols, value_col):
    return df.groupby(group_cols)[value_col].sum().reset_index()

def parse_output(exon_overlap_file, intron_overlap_file, gene_overlap_file, dog_overlap_file, bed_file, end_coordinates, record_exons, genome_file, output_dir, num_files, original_exon_bed, original_intron_bed, original_dog_bed):
 
    # check input files 
    
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
        warnings.warn(f"Warning: The following file(s) are empty. Some data may not be processed: {', '.join(empty_files)}")
            
    exon_cols = ['read_chrom', 'read_fragment_start', 'read_fragment_end', 'read_id', 'read-alignment_length', 'read_strand', 'exon_chrom', 'exon_start', 'exon_end', 'exon_gene_id', 'exon_id', 'exon_strand', 'exon_base_overlap']
    intron_cols = ['read_chrom', 'read_fragment_start', 'read_fragment_end', 'read_id', 'read-alignment_length', 'read_strand', 'intron_chrom', 'intron_start', 'intron_end', 'intron_gene_id', 'intron_id', 'intron_strand', 'intron_base_overlap']
    gene_cols = ['read_chrom', 'read_fragment_start', 'read_fragment_end', 'read_id', 'read-alignment_length', 'read_strand', 'gene_chrom', 'gene_start', 'gene_end', 'gene_id', 'gene_biotype', 'gene_strand', 'gene_base_overlap']
    dog_cols = ['read_chrom', 'read_fragment_start', 'read_fragment_end', 'read_id', 'read-alignment_length', 'read_strand', 'dog_chrom', 'dog_start', 'dog_end', 'dog_gene_id', 'dog_id', 'dog_strand', 'dog_base_overlap']
    bam_to_bed_cols = ['chrom', 'start', 'end', 'name', 's-a_length', 'strand']

    # read overlap files 

    try:
        exon_df = pd.read_csv(exon_overlap_file, sep="\t", header=None, names=exon_cols, low_memory=False)
        intron_df = pd.read_csv(intron_overlap_file, sep="\t", header=None, names=intron_cols, low_memory=False)
        gene_df = pd.read_csv(gene_overlap_file, sep="\t", header=None, names=gene_cols, low_memory=False)
        bam_df = pd.read_csv(bed_file, sep="\t", header=None, names=bam_to_bed_cols, low_memory=False)

        # prepare for empty DOG df, e.g. when handling total RNA 
        if os.stat(dog_overlap_file).st_size == 0:
            print("Warning: DOG overlap file is empty. Proceeding without DOG information.")
            dog_df = pd.DataFrame(columns=dog_cols)
        else:
            dog_df = pd.read_csv(dog_overlap_file, sep="\t", header=None, names=dog_cols, low_memory=False)
    except IOError as e:
        raise Exception(f"Error reading input files: {str(e)}")

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

        # if the most upstream feature is an exon, we're done
        if most_upstream['type'] == 'exon':
            return most_upstream['gene_id']

        # check if there are any exon or intron alignments
        exon_intron_features = [f for f in read_data if f['type'] in ['exon', 'intron']]
        if exon_intron_features:
            # if there are multiple genes with equal exon/intron alignment, use total overlap as tiebreaker
            gene_total_overlap = defaultdict(int)
            for f in exon_intron_features:
                gene_total_overlap[f['gene_id']] += f['overlap']
            return max(gene_total_overlap, key=gene_total_overlap.get)

        # if no exon or intron alignments, check for DOG alignments
        dog_features = [f for f in read_data if f['type'] == 'dog']
        if dog_features:
            # if there are DOG alignments, use the one with the most overlap
            return max(dog_features, key=lambda x: x['overlap'])['gene_id']

        # if no alignments at all, return none
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

    # assign genes
    top_gene = {read_id: assign_gene(data) for read_id, data in read_data.items()}
    top_gene = {k: v for k, v in top_gene.items() if v is not None}

    # process overlaps in parallel 
    num_cores = multiprocessing.cpu_count()
    with ProcessPoolExecutor(max_workers=num_cores) as executor:
        exon_overlap_future = executor.submit(process_overlap_group, exon_df[exon_df.apply(lambda row: top_gene.get(row['read_id']) == row['exon_gene_id'], axis=1)], ['read_id', 'exon_gene_id'], 'exon_base_overlap')
        intron_overlap_future = executor.submit(process_overlap_group, intron_df[intron_df.apply(lambda row: top_gene.get(row['read_id']) == row['intron_gene_id'], axis=1)], ['read_id', 'intron_gene_id'], 'intron_base_overlap')
        gene_overlap_future = executor.submit(process_overlap_group, gene_df[gene_df.apply(lambda row: top_gene.get(row['read_id']) == row['gene_id'], axis=1)], ['read_id', 'gene_id'], 'gene_base_overlap')
        dog_overlap_future = executor.submit(process_overlap_group, dog_df[dog_df.apply(lambda row: top_gene.get(row['read_id']) == row['dog_gene_id'], axis=1)], ['read_id', 'dog_gene_id'], 'dog_base_overlap')

    exon_overlap_group = exon_overlap_future.result()
    intron_overlap_group = intron_overlap_future.result()
    gene_overlap_group = gene_overlap_future.result()
    dog_overlap_group = dog_overlap_future.result()

    # merge overlaps 
    all_overlap_groups = [exon_overlap_group, intron_overlap_group, gene_overlap_group, dog_overlap_group]
    gene_overlap_sum = reduce(lambda left, right: pd.merge(left, right, on='read_id', how='outer'), all_overlap_groups)

    # handle NaN overlaps 
    overlap_columns = ['exon_base_overlap', 'intron_base_overlap', 'gene_base_overlap', 'dog_base_overlap']
    gene_overlap_sum[overlap_columns] = gene_overlap_sum[overlap_columns].fillna(0)

    # ensure gene_id is consistent 
    gene_overlap_sum['gene_id'] = gene_overlap_sum.apply(lambda row: top_gene.get(row['read_id'], np.nan), axis=1)

    # merge with bam_df to include all reads
    bam_df.rename(columns={'name': 'read_id', 's-a_length': 'read-alignment_length'}, inplace=True)
    bam_df = bam_df.drop_duplicates(subset='read_id')
    gene_overlap_sum = pd.merge(bam_df[['read_id', 'read-alignment_length']], gene_overlap_sum, how='left', on='read_id')

    gene_overlap_sum[['read_length', 'alignment_length', 'splice_count']] = gene_overlap_sum['read-alignment_length'].str.split(',', expand=True)
    gene_overlap_sum['read_length'] = pd.to_numeric(gene_overlap_sum['read_length'], errors='coerce').fillna(0).astype(int)
    gene_overlap_sum['alignment_length'] = pd.to_numeric(gene_overlap_sum['alignment_length'], errors='coerce').fillna(0).astype(int)
    gene_overlap_sum['splice_count'] = pd.to_numeric(gene_overlap_sum['splice_count'], errors='coerce').fillna(0).astype(int)
    gene_overlap_sum.drop(columns='read-alignment_length', inplace=True)

    # unclassified and unaligned lengths
    gene_overlap_sum['unclassified_length'] = gene_overlap_sum['alignment_length'] - gene_overlap_sum['gene_base_overlap'] - gene_overlap_sum['dog_base_overlap']
    gene_overlap_sum['unaligned_length'] = gene_overlap_sum['read_length'] - gene_overlap_sum['alignment_length']

    # end coordinates
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

    def create_gene_chunks(df, num_files):

        gene_groups = df.groupby('gene_id')
        total_reads = len(df)
        target_chunk_size = total_reads // num_files
        
        chunks = defaultdict(list)
        current_chunk = 0
        current_chunk_size = 0
        
        for gene_id, group in gene_groups:
            gene_size = len(group)
            
            if current_chunk_size + gene_size > target_chunk_size and current_chunk < num_files - 1:
                current_chunk += 1
                current_chunk_size = 0
            
            chunks[current_chunk].append(gene_id)
            current_chunk_size += gene_size
        
        return chunks

    def create_and_sort_3prime_bed(df, output_dir, gene_chunks):
        temp_bed_files = []
        
        def process_chunk(chunk_genes, i):
            chunk = df[df['gene_id'].isin(chunk_genes)].sort_values(by=['gene_id', 'read_end_chromosome', 'read_end_position'])
            temp_file = tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=f'_3prime_{i}.sorted.bed', dir=output_dir)
            for _, row in chunk.iterrows():
                temp_file.write(f"{row['gene_id']}_{row['read_end_chromosome']}\t{row['read_end_position']}\t{row['read_end_position']+1}\t{row['read_id']}\t.\t{row['read_end_strand']}\n")
            temp_file.close()
            return temp_file.name

        with ThreadPoolExecutor() as executor:
            futures = []
            for i, chunk_genes in gene_chunks.items():
                futures.append(executor.submit(process_chunk, chunk_genes, i))
            
            temp_bed_files = [future.result() for future in futures]
        
        return temp_bed_files

    def create_and_sort_feature_bed(original_exon_bed, original_intron_bed, original_dog_bed, output_dir, gene_chunks):
        temp_bed_files = []
        
        def process_chunk(chunk_genes, i):
            temp_file = tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=f'_features_{i}.bed', dir=output_dir)
            for bed_file, feature_type in [(original_exon_bed, 'exon'), (original_intron_bed, 'intron'), (original_dog_bed, 'DOG')]:
                with open(bed_file, 'r') as f:
                    for line in f:
                        fields = line.strip().split('\t')
                        gene_id = fields[3]
                        if gene_id in chunk_genes:
                            chrom, start, end, feature_id, score, strand = fields
                            temp_file.write(f"{gene_id}_{chrom}\t{start}\t{end}\t{feature_type}_{feature_id}\t{score}\t{strand}\n")
            temp_file.close()
            
            sorted_file = f"{temp_file.name}.sorted"
            subprocess.run(f"sort -k1,1 -k2,2n {temp_file.name} > {sorted_file}", shell=True, check=True)
            os.unlink(temp_file.name)  # remove unsorted file 
            return sorted_file

        with ThreadPoolExecutor() as executor:
            futures = []
            for i, chunk_genes in gene_chunks.items():
                futures.append(executor.submit(process_chunk, chunk_genes, i))
            
            temp_bed_files = [future.result() for future in futures]
        
        return temp_bed_files

    # create gene chunks 
    gene_chunks = create_gene_chunks(gene_overlap_sum, num_files)

    # create and sort BED files
    three_prime_bed_files = create_and_sort_3prime_bed(gene_overlap_sum, output_dir, gene_chunks)
    feature_bed_files = create_and_sort_feature_bed(original_exon_bed, original_intron_bed, original_dog_bed, output_dir, gene_chunks)

    # intersect read 3' coordinates to features from their assigned genes 
    with ThreadPoolExecutor(max_workers=num_files) as executor:
        tasks = []
        for i in range(len(gene_chunks)):
            intersect_output = f"{output_dir}/three_prime_feature_intersect_{i}.bed"
            intersect_cmd = f"bedtools intersect -a {three_prime_bed_files[i]} -b {feature_bed_files[i]} -sorted -wo > {intersect_output}"
            tasks.append(executor.submit(subprocess.run, intersect_cmd, shell=True, check=True))

        for task in as_completed(tasks):
            task.result()

    # process results 
    feature_map = {}
    error_cases = []

    for i in range(len(gene_chunks)):
        intersect_output = f"{output_dir}/three_prime_feature_intersect_{i}.bed"
        with open(intersect_output, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                read_id = fields[3]
                feature_type = fields[9].split('_')[0]
                if read_id in feature_map:
                    error_cases.append((read_id, fields[0].split('_')[0]))
                else:
                    parts = fields[9].split('_')
                    feature_map[read_id] = '_'.join(fields[10].split('_')[1:])

    # add feature information 
    gene_overlap_sum['three_prime_feature'] = gene_overlap_sum['read_id'].map(lambda x: feature_map.get(x, "unclassified"))

    # log error cases for 3' coordinate features 
    with open(f"{output_dir}/three_prime_feature_errors.log", "w") as f:
        f.write(f"Number of reads with multiple 3' end features: {len(error_cases)}\n")
        for read_id, gene_id in error_cases:
            f.write(f"Read ID: {read_id}, Gene ID: {gene_id}\n")

    # if record_exons, record exon and intron IDs with overlaps in the target gene
    # NOTE: we could consider adding a threshold to introns, since noisy reads 
    # can leak into the intron and cause spurious classification of intron spanning reads 
    if record_exons:
        read_exon_ids = defaultdict(lambda: defaultdict(set))
        read_intron_ids = defaultdict(lambda: defaultdict(set))

        for _, row in exon_df.iterrows():
            
            # check for valid gene_id assignment
            if top_gene.get(row['read_id']):
                read_exon_ids[row['read_id']][row['exon_gene_id']].add(row['exon_id'])

        for _, row in intron_df.iterrows():
            
            # check for valid gene_id assignment
            if top_gene.get(row['read_id']):
                read_intron_ids[row['read_id']][row['intron_gene_id']].add(row['intron_id'])

        def collapse_ids(id_set):
            return ','.join(sorted(set(id.split('_')[-1] for id in id_set)))

        gene_overlap_sum['exon_ids'] = gene_overlap_sum.apply(
            lambda row: collapse_ids(read_exon_ids.get(row['read_id'], {}).get(row['gene_id'], [])), axis=1
        )
        gene_overlap_sum['intron_ids'] = gene_overlap_sum.apply(
            lambda row: collapse_ids(read_intron_ids.get(row['read_id'], {}).get(row['gene_id'], [])), axis=1
        )

    # update col order 
    column_order = ['read_id', 'gene_id', 'read_length', 'alignment_length', 'splice_count', 
                    'gene_base_overlap', 'exon_base_overlap', 'intron_base_overlap', 'gene_exon_bases', 'DOG_overlap', 
                    'unclassified_length', 'unaligned_length', 'read_end_chromosome', 
                    'read_end_position', 'read_end_strand', 'strand_sign', 'three_prime_feature']
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