import pandas as pd
import numpy as np
import argparse
import os
import logging
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
from collections import defaultdict

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def read_bed(file_path):
    df = pd.read_csv(file_path, sep='\t', header=None, 
                     names=['chrom', 'start', 'end', 'gene_id', 'score', 'strand'],
                     dtype={'chrom': 'category', 'start': np.int32, 'end': np.int32, 'gene_id': 'category', 'score': 'category', 'strand': 'category'})
    return df

def read_clusters(file_path):
    clusters_df = pd.read_csv(file_path, sep='\t', header=None, 
                              names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'distance'],
                              dtype={'chrom': 'category', 'start': np.int32, 'end': np.int32, 'name': str, 'score': np.int32, 'strand': 'category', 'distance': np.float32})
    clusters_df[['gene_id', 'category']] = clusters_df['name'].str.split('_', n=1, expand=True)
    clusters_df['gene_id'] = clusters_df['gene_id'].astype('category')
    clusters_df['category'] = clusters_df['category'].astype('category')
    return clusters_df

def read_fasta_index(fai_path):
    fai_df = pd.read_csv(fai_path, sep='\t', header=None, names=['chrom', 'length', 'offset', 'linebases', 'linewidth'])
    return fai_df.set_index('chrom')['length'].to_dict()

def find_downstream_positions(clusters_df, gene_df):
    # start with post-transcriptional cluster 
    post_trans_clusters = clusters_df[clusters_df['category'] == 'post_transcriptional']
    gene_info = gene_df.set_index('gene_id')[['chrom', 'start', 'end', 'strand']]
    
    downstream_positions = {}
    
    # don't return genes that don't have polyA clusters 
    for gene_id, group in post_trans_clusters.groupby('gene_id', observed=True):
        if gene_id not in gene_info.index:
            continue
        
        gene_chrom, gene_start, gene_end, gene_strand = gene_info.loc[gene_id]
        
        valid_clusters = group[group['strand'] == gene_strand]
        
        # only considers clusters downstream of the TSS 
        if gene_strand == '+':
            valid_clusters = valid_clusters[valid_clusters['start'] >= gene_start]
        else:  
            valid_clusters = valid_clusters[valid_clusters['end'] <= gene_end]
        
        if not valid_clusters.empty:
            downstream_positions[gene_id] = valid_clusters.to_dict('records')
    
    return downstream_positions

def number_exons_and_introns(exons, introns, gene_id, strand):
    if exons.empty and introns.empty:
        return pd.DataFrame(), pd.DataFrame()

    if not exons.empty and not introns.empty:
        regions = pd.concat([exons, introns], ignore_index=True)
        region_types = ['exon'] * len(exons) + ['intron'] * len(introns)
    elif not exons.empty:
        regions = exons.copy().reset_index(drop=True)
        region_types = ['exon'] * len(exons)
    else:
        regions = introns.copy().reset_index(drop=True)
        region_types = ['intron'] * len(introns)
    
    regions['region_type'] = region_types
    
    if strand == '+':
        regions = regions.sort_values('start')
    else:
        regions = regions.sort_values('start', ascending=False)
    
    exon_count = 1
    intron_count = 1
    numbers = []
    for _, row in regions.iterrows():
        if row['region_type'] == 'exon':
            numbers.append(f"exon_{exon_count}")
            exon_count += 1
        else:
            numbers.append(f"intron_{intron_count}")
            intron_count += 1
    
    regions['number'] = numbers
    regions['combined'] = gene_id + '_' + regions['number']
    
    numbered_exons = regions[regions['region_type'] == 'exon'].drop('region_type', axis=1)
    numbered_introns = regions[regions['region_type'] == 'intron'].drop('region_type', axis=1)
    
    return numbered_exons, numbered_introns

def create_initial_dog(gene_entry, pas_position, dog_length=2500):
    if gene_entry['strand'] == '+':
        dog_start = pas_position
        dog_end = pas_position + dog_length
    else:
        dog_start = pas_position - dog_length
        dog_end = pas_position

    return pd.DataFrame({
        'chrom': [gene_entry['chrom']],
        'start': [dog_start],
        'end': [dog_end],
        'gene_id': [gene_entry['gene_id']],
        'score': [gene_entry['gene_id'] + "_region_DOG"],
        'strand': [gene_entry['strand']]
    })

def adjust_for_genome_boundaries(dog, chrom_lengths):
    if dog.empty:
        return dog
    
    chrom_length = chrom_lengths.get(dog['chrom'].iloc[0], float('inf'))
    dog['start'] = dog['start'].clip(lower=0)
    dog['end'] = dog['end'].clip(upper=chrom_length)
    
    return dog[dog['start'] < dog['end']]  

def adjust_regions_for_gene(gene_id, exon_df, intron_df, gene_df, cluster, chrom_lengths):
    gene_exons = exon_df[exon_df['gene_id'] == gene_id].sort_values('start')
    gene_introns = intron_df[intron_df['gene_id'] == gene_id].sort_values('start')
    gene_entry = gene_df[gene_df['gene_id'] == gene_id].iloc[0]
    
    # record the PAS adjustment distance
    pas_adjustment = 0

    # if no cluster is found, keep the original gene structure
    if cluster is None:
        original_dog = create_initial_dog(gene_entry, gene_entry['end'] if gene_entry['strand'] == '+' else gene_entry['start'])
        original_dog = adjust_for_genome_boundaries(original_dog, chrom_lengths)
        numbered_exons, numbered_introns = number_exons_and_introns(gene_exons, gene_introns, gene_id, gene_entry['strand'])
        return numbered_exons, numbered_introns, original_dog, pd.DataFrame([gene_entry]), pd.DataFrame(), "No cluster provided", pas_adjustment

    # select a polyA cluster 
    pas_position = cluster['end'] if cluster['strand'] == '+' else cluster['start']
    
    # discard genes where the most downstream cluster is before/after gene boundary
    if (cluster['strand'] == '+' and pas_position < gene_entry['start']) or \
       (cluster['strand'] == '-' and pas_position > gene_entry['end']):
        error_message = (f"PAS position {pas_position} is out of range for gene {gene_id}")
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame([gene_entry]), pd.DataFrame([cluster]), error_message, pas_adjustment

    # trim or extend the gene to the new PAS
    new_gene_entry = gene_entry.copy()
    if cluster['strand'] == '+':
        # record adjustment distance
        pas_adjustment = abs(gene_entry['end'] - pas_position)
        new_gene_entry['end'] = pas_position
    else:
        pas_adjustment = abs(gene_entry['start'] - pas_position)
        new_gene_entry['start'] = pas_position

    # delete all exons and introns downstream of the PAS
    if cluster['strand'] == '+':
        exons_to_keep = gene_exons[gene_exons['start'] < pas_position]
        introns_to_keep = gene_introns[gene_introns['start'] < pas_position]
    else:
        exons_to_keep = gene_exons[gene_exons['end'] > pas_position]
        introns_to_keep = gene_introns[gene_introns['end'] > pas_position]

    # handle exons/introns that span the PAS
    spanning_exon = exons_to_keep[(exons_to_keep['start'] < pas_position) & (exons_to_keep['end'] > pas_position)]
    spanning_intron = introns_to_keep[(introns_to_keep['start'] < pas_position) & (introns_to_keep['end'] > pas_position)]

    if not spanning_intron.empty:
        introns_to_keep = introns_to_keep.drop(spanning_intron.index)
        # extend/trim the relevant exon
        if cluster['strand'] == '+':
            last_exon = exons_to_keep.iloc[-1]
            if last_exon['end'] < pas_position:
                exons_to_keep.loc[last_exon.name, 'end'] = pas_position
        else:
            first_exon = exons_to_keep.iloc[0]
            if first_exon['start'] > pas_position:
                exons_to_keep.loc[first_exon.name, 'start'] = pas_position
    elif not spanning_exon.empty:
        # trim the spanning exon
        exon_index = spanning_exon.index[0]
        if cluster['strand'] == '+':
            exons_to_keep.loc[exon_index, 'end'] = pas_position 
        else:
            exons_to_keep.loc[exon_index, 'start'] = pas_position
    else:
        # no spanning element, just ensure we extend the exon if needed
        if cluster['strand'] == '+':
            last_exon = exons_to_keep.iloc[-1]
            if last_exon['end'] < pas_position:
                exons_to_keep.loc[last_exon.name, 'end'] = pas_position
        else:
            first_exon = exons_to_keep.iloc[0]
            if first_exon['start'] > pas_position:
                exons_to_keep.loc[first_exon.name, 'start'] = pas_position

    # continuity check
    total_exon_intron_length = (exons_to_keep['end'] - exons_to_keep['start']).sum() + \
                               (introns_to_keep['end'] - introns_to_keep['start']).sum()
    adjusted_gene_length = new_gene_entry['end'] - new_gene_entry['start']
    if total_exon_intron_length != adjusted_gene_length:
        logging.debug(f"Original gene: start={gene_entry['start']}, end={gene_entry['end']}, length={gene_entry['end'] - gene_entry['start']}")
        logging.debug(f"PAS position: {pas_position}")
        logging.debug(f"Strand: {cluster['strand']}")
        logging.debug(f"Exons to keep after adjustment:\n{exons_to_keep}")
        logging.debug(f"Introns to keep after adjustment:\n{introns_to_keep}")
        logging.debug(f"New gene entry: start={new_gene_entry['start']}, end={new_gene_entry['end']}, length={adjusted_gene_length}")
        
        raise ValueError(f"Adjusted gene {gene_id} is not continuous. "
                         f"Exon+Intron length ({total_exon_intron_length}) != Gene length ({adjusted_gene_length})")

    new_dog = create_initial_dog(new_gene_entry, pas_position)
    new_dog = adjust_for_genome_boundaries(new_dog, chrom_lengths)

    # final checks
    error_message = None
    if new_dog.empty:
        error_message = f"DOG region is empty for gene {gene_id}"
    else:
        if len(new_dog) > 1:
            error_message = f"More than one DOG region found for gene {gene_id}"
        elif new_dog['strand'].iloc[0] != gene_entry['strand']:
            error_message = f"DOG strand does not match gene strand for gene {gene_id}"
        elif (cluster['strand'] == '+' and new_dog['start'].iloc[0] != pas_position) or \
             (cluster['strand'] == '-' and new_dog['end'].iloc[0] != pas_position):
            error_message = f"DOG is not adjacent to PAS for gene {gene_id}"

    numbered_exons, numbered_introns = number_exons_and_introns(exons_to_keep, introns_to_keep, gene_id, new_gene_entry['strand'])

    return numbered_exons, numbered_introns, new_dog, pd.DataFrame([new_gene_entry]), pd.DataFrame([cluster]), error_message, pas_adjustment

def process_gene_chunk(chunk_data):
    exon_df, intron_df, gene_df, downstream_positions_dict, chrom_lengths, gene_ids_chunk = chunk_data
    chunk_exons, chunk_introns, chunk_dogs, chunk_genes, chunk_clusters = [], [], [], [], []
    discordant_pas = []
    cross_chromosome_clusters = defaultdict(int)
    genes_without_valid_clusters = []
    pas_adjustments = []

    for gene_id in gene_ids_chunk:
        try:
            gene_info = gene_df[gene_df['gene_id'] == gene_id].iloc[0]
            gene_chrom = gene_info['chrom']
            gene_strand = gene_info['strand']

            gene_clusters = downstream_positions_dict.get(gene_id, [])
            cross_chrom_count = sum(1 for c in gene_clusters if c['chrom'] != gene_chrom)
            if cross_chrom_count > 0:
                cross_chromosome_clusters[gene_id] = cross_chrom_count

            valid_clusters = [c for c in gene_clusters if c['chrom'] == gene_chrom and c['strand'] == gene_strand]
            cluster = None
            if valid_clusters:
                for c in valid_clusters:
                    c['score'] = int(c['score'])
                cluster = max(valid_clusters, key=lambda x: x['score'])

            exon, intron, dog, gene, cluster_data, error_message, pas_adjustment = adjust_regions_for_gene(
                gene_id, exon_df, intron_df, gene_df, cluster, chrom_lengths)

            chunk_exons.append(exon)
            chunk_introns.append(intron)
            chunk_dogs.append(dog)
            chunk_genes.append(gene)
            if cluster is not None:
                chunk_clusters.append(cluster_data)

            if error_message:
                discordant_pas.append(error_message)

            if cluster is None:
                genes_without_valid_clusters.append(gene_id)

            if pas_adjustment > 0:
                pas_adjustments.append(pas_adjustment)

        except Exception as e:
            # if error, write debug info and re-raise to stop
            error_log_path = os.path.join(os.environ.get('OUTPUT_DIR', '.'), f'error_case_{gene_id}.log')
            with open(error_log_path, 'w') as f:
                f.write(f"Error for gene {gene_id}:\n")
                f.write(str(e) + "\n")
                f.write("Debug info:\n")
                f.write(f"gene_info:\n{gene_info}\n")
                f.write(f"gene_clusters:\n{gene_clusters}\n")
            raise e

    return (
        pd.concat(chunk_exons, ignore_index=True) if chunk_exons else pd.DataFrame(),
        pd.concat(chunk_introns, ignore_index=True) if chunk_introns else pd.DataFrame(),
        pd.concat(chunk_dogs, ignore_index=True) if chunk_dogs else pd.DataFrame(),
        pd.concat(chunk_genes, ignore_index=True) if chunk_genes else pd.DataFrame(),
        pd.concat(chunk_clusters, ignore_index=True) if chunk_clusters else pd.DataFrame(),
        discordant_pas,
        dict(cross_chromosome_clusters),
        genes_without_valid_clusters,
        pas_adjustments
    )

def main():
    parser = argparse.ArgumentParser(description="Process genomic intervals and create updated BED files.")
    parser.add_argument("--gene", required=True, help="Path to gene BED file")
    parser.add_argument("--exon", required=True, help="Path to exon BED file")
    parser.add_argument("--intron", required=True, help="Path to intron BED file")
    parser.add_argument("--clusters", required=True, help="Path to polyA clusters BED file")
    parser.add_argument("--fai", required=True, help="Path to genome FASTA index file")
    parser.add_argument("--output_dir", required=True, help="Directory to save output files")
    parser.add_argument("--cores", type=int, default=os.cpu_count(), help="Number of CPU cores to use")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    os.environ['OUTPUT_DIR'] = args.output_dir

    exon_df = read_bed(args.exon)
    intron_df = read_bed(args.intron)
    gene_df = read_bed(args.gene)
    clusters_df = read_clusters(args.clusters)
    chrom_lengths = read_fasta_index(args.fai)

    initial_gene_count = len(gene_df['gene_id'].unique())
    logging.info(f"Initial number of genes: {initial_gene_count}")

    downstream_positions_dict = find_downstream_positions(clusters_df, gene_df)

    if not downstream_positions_dict:
        logging.warning("No valid downstream positions found. Outputting original annotations.")
        exon_df.to_csv(os.path.join(args.output_dir, 'updated_exon.bed'), sep='\t', header=False, index=False)
        intron_df.to_csv(os.path.join(args.output_dir, 'updated_intron.bed'), sep='\t', header=False, index=False)
        gene_df.to_csv(os.path.join(args.output_dir, 'updated_gene.bed'), sep='\t', header=False, index=False)
        pd.DataFrame(columns=['chrom', 'start', 'end', 'gene_id', 'score', 'strand']).to_csv(
            os.path.join(args.output_dir, 'updated_downstream_of_gene.bed'), sep='\t', header=False, index=False)
        pd.DataFrame(columns=['chrom', 'start', 'end', 'name', 'score', 'strand', 'distance']).to_csv(
            os.path.join(args.output_dir, 'PAS.bed'), sep='\t', header=False, index=False)
        return

    gene_ids = gene_df['gene_id'].unique()
    gene_chunks = np.array_split(gene_ids, args.cores)
    chunk_data = [(exon_df, intron_df, gene_df, downstream_positions_dict, chrom_lengths, gene_chunk) for gene_chunk in gene_chunks]

    try:
        with ProcessPoolExecutor(max_workers=args.cores) as executor:
            results = list(tqdm(executor.map(process_gene_chunk, chunk_data), total=len(gene_chunks), desc="Processing genes"))
    except Exception as e:
        logging.error("An error occurred during processing. Check the error_case_* logs for details.")
        raise SystemExit(1)

    # unpack results
    final_exons, final_introns, final_dogs, final_genes, final_clusters, all_discordant_pas, all_cross_chromosome, all_genes_without_valid_clusters, all_pas_adjustments = map(list, zip(*results))
  
    final_exons = [df for df in final_exons if not df.empty]
    final_introns = [df for df in final_introns if not df.empty]
    final_dogs = [df for df in final_dogs if not df.empty]
    final_genes = [df for df in final_genes if not df.empty]
    final_clusters = [df for df in final_clusters if not df.empty]

    final_exons = pd.concat(final_exons) if final_exons else pd.DataFrame()
    final_introns = pd.concat(final_introns) if final_introns else pd.DataFrame()
    final_dogs = pd.concat(final_dogs) if final_dogs else pd.DataFrame()
    final_genes = pd.concat(final_genes) if final_genes else pd.DataFrame()
    final_clusters = pd.concat(final_clusters) if final_clusters else pd.DataFrame()

    all_discordant_pas = [item for sublist in all_discordant_pas for item in sublist]
    all_genes_without_valid_clusters = [item for sublist in all_genes_without_valid_clusters for item in sublist]
    all_pas_adjustments = [item for sublist in all_pas_adjustments for item in sublist]

    merged_cross_chromosome = {}
    for d in all_cross_chromosome:
        for gene, count in d.items():
            merged_cross_chromosome[gene] = merged_cross_chromosome.get(gene, 0) + count

    final_gene_count = len(final_genes['gene_id'].unique())
    logging.info(f"Final number of genes: {final_gene_count}")

    if final_gene_count != initial_gene_count:
        logging.error(f"Gene count mismatch: Started with {initial_gene_count}, ended with {final_gene_count}")
        raise ValueError("Gene count mismatch")

    final_exons.to_csv(os.path.join(args.output_dir, 'updated_exon.bed'), sep='\t', header=False, index=False, 
                       columns=['chrom', 'start', 'end', 'gene_id', 'combined', 'strand'])
    final_introns.to_csv(os.path.join(args.output_dir, 'updated_intron.bed'), sep='\t', header=False, index=False, 
                         columns=['chrom', 'start', 'end', 'gene_id', 'combined', 'strand'])
    final_dogs.to_csv(os.path.join(args.output_dir, 'updated_downstream_of_gene.bed'), sep='\t', header=False, index=False)
    final_genes.to_csv(os.path.join(args.output_dir, 'updated_gene.bed'), sep='\t', header=False, index=False)
    final_clusters.to_csv(os.path.join(args.output_dir, 'PAS.bed'), sep='\t', header=False, index=False)

    # Compute statistics on PAS adjustments
    if all_pas_adjustments:
        pas_series = pd.Series(all_pas_adjustments)
        mean_adj = pas_series.mean()
        median_adj = pas_series.median()
        q25 = pas_series.quantile(0.25)
        q75 = pas_series.quantile(0.75)
        min_adj = pas_series.min()
        max_adj = pas_series.max()
        total_adjusted_genes = len(all_pas_adjustments)
    else:
        mean_adj = median_adj = q25 = q75 = min_adj = max_adj = np.nan
        total_adjusted_genes = 0

    log_file = os.path.join(args.output_dir, 'processing_log.txt')
    with open(log_file, 'w') as f:
        f.write(f"Initial number of genes: {initial_gene_count}\n")
        f.write(f"Final number of genes: {final_gene_count}\n\n")

        f.write(f"Total discordant PAS: {len(all_discordant_pas)}\n\n")
        for pas in all_discordant_pas:
            f.write(f"{pas}\n")
        
        f.write(f"\nTotal genes with cross-chromosome clusters: {len(merged_cross_chromosome)}\n")
        for gene, count in merged_cross_chromosome.items():
            f.write(f"Gene {gene}: {count} cluster(s) on different chromosome(s)\n")
        
        f.write(f"\nTotal genes without valid clusters: {len(all_genes_without_valid_clusters)}\n")
        for gene in all_genes_without_valid_clusters:
            f.write(f"Gene without valid clusters: {gene}\n")

        f.write("\nPAS Adjustment Statistics:\n")
        f.write(f"Total adjusted genes: {total_adjusted_genes}\n")
        f.write(f"Mean adjustment: {mean_adj}\n")
        f.write(f"Median adjustment: {median_adj}\n")
        f.write(f"25% quantile: {q25}\n")
        f.write(f"75% quantile: {q75}\n")
        f.write(f"Min adjustment: {min_adj}\n")
        f.write(f"Max adjustment: {max_adj}\n")

    logging.info("All updated BED files have been created.")
    logging.info(f"Found {len(all_discordant_pas)} discordant PAS events.")
    logging.info(f"Found {len(merged_cross_chromosome)} genes with cross-chromosome clusters.")
    logging.info(f"Found {len(all_genes_without_valid_clusters)} genes without valid clusters.")
    logging.info(f"Total adjusted genes: {total_adjusted_genes}, mean={mean_adj}, median={median_adj}, q25={q25}, q75={q75}, min={min_adj}, max={max_adj}")
    logging.info(f"Details logged in {log_file}")

if __name__ == "__main__":
    main()
