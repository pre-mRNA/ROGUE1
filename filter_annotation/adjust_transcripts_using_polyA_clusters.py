import pandas as pd
import numpy as np
import argparse
import os
import logging
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
from collections import defaultdict

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

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
    post_trans_clusters = clusters_df[clusters_df['category'] == 'post_transcriptional']
    gene_info = gene_df.set_index('gene_id')[['chrom', 'start', 'end', 'strand']]
    
    downstream_positions = {}
    
    for gene_id, group in post_trans_clusters.groupby('gene_id'):
        if gene_id not in gene_info.index:
            continue
        
        gene_chrom, gene_start, gene_end, gene_strand = gene_info.loc[gene_id]
        
        valid_clusters = group[group['strand'] == gene_strand]
        
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

    regions = pd.concat([exons, introns]).reset_index(drop=True)
    regions['region_type'] = ['exon'] * len(exons) + ['intron'] * len(introns)
    
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

def subtract_regions(dog, regions):
    if dog.empty or regions.empty:
        return dog

    dog_intervals = pd.IntervalIndex.from_arrays(dog['start'], dog['end'], closed='both')
    region_intervals = pd.IntervalIndex.from_arrays(regions['start'], regions['end'], closed='both')

    remaining_intervals = dog_intervals.difference(region_intervals)
    
    if len(remaining_intervals) == 0:
        return pd.DataFrame(columns=dog.columns)
    
    new_dog = pd.DataFrame({
        'chrom': dog['chrom'].iloc[0],
        'start': [interval.left for interval in remaining_intervals],
        'end': [interval.right for interval in remaining_intervals],
        'gene_id': dog['gene_id'].iloc[0],
        'score': dog['score'].iloc[0],
        'strand': dog['strand'].iloc[0]
    })
    
    return new_dog

def trim_split_dog(dog):
    if len(dog) <= 1:
        return dog
    return dog.iloc[[0]]  # if dog is split, e.g by another exon, keep only the most upstream dog region 

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
    
    if cluster is None:
        exons_to_keep = gene_exons
        introns_to_keep = gene_introns
        new_dog = pd.DataFrame()
        new_gene_entry = gene_entry
        error_message = "No cluster provided"
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame([gene_entry]), pd.DataFrame(), error_message

    pas_position = cluster['end'] if cluster['strand'] == '+' else cluster['start']
    
    if (cluster['strand'] == '+' and pas_position < gene_entry['start']) or \
       (cluster['strand'] == '-' and pas_position > gene_entry['end']):
        error_message = f"PAS position {pas_position} is {'before' if cluster['strand'] == '+' else 'after'} gene {'start' if cluster['strand'] == '+' else 'end'} {gene_entry['start' if cluster['strand'] == '+' else 'end']} for gene {gene_id}"
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame([gene_entry]), pd.DataFrame([cluster]), error_message

    # adjust exons and introns based on PAS position
    affected_exons = gene_exons[(gene_exons['start'] <= pas_position) & (gene_exons['end'] >= pas_position)]
    affected_exons_indices = affected_exons.index
    if cluster['strand'] == '+':
        gene_exons.loc[affected_exons_indices, 'end'] = pas_position
        exons_to_keep = gene_exons[gene_exons['end'] <= pas_position]
        introns_to_keep = gene_introns[gene_introns['end'] < pas_position]
    else:
        gene_exons.loc[affected_exons_indices, 'start'] = pas_position
        exons_to_keep = gene_exons[gene_exons['start'] >= pas_position]
        introns_to_keep = gene_introns[gene_introns['start'] > pas_position]

    new_gene_entry = gene_entry.copy()
    new_gene_entry['end'] = pas_position if cluster['strand'] == '+' else gene_entry['end']
    new_gene_entry['start'] = gene_entry['start'] if cluster['strand'] == '+' else pas_position

    # create DOG regions from scratch and adjust for genome boundaries
    new_dog = create_initial_dog(new_gene_entry, pas_position)
    new_dog = adjust_for_genome_boundaries(new_dog, chrom_lengths)

    numbered_exons, numbered_introns = number_exons_and_introns(exons_to_keep, introns_to_keep, gene_id, new_gene_entry['strand'])

    error_message = None
    if not new_dog.empty:
        if len(new_dog) > 1:
            error_message = "More than one DOG region found for gene " + gene_id
        elif new_dog['strand'].iloc[0] != gene_entry['strand']:
            error_message = "DOG strand does not match gene strand for gene " + gene_id
        elif (cluster['strand'] == '+' and new_dog['start'].iloc[0] != pas_position) or \
             (cluster['strand'] == '-' and new_dog['end'].iloc[0] != pas_position):
            error_message = "DOG is not adjacent to PAS for gene " + gene_id + " (" + ('positive' if cluster['strand'] == '+' else 'negative') + " strand)"

    return numbered_exons, numbered_introns, new_dog, pd.DataFrame([new_gene_entry]), pd.DataFrame([cluster]), error_message


def process_gene_chunk(chunk_data):
    exon_df, intron_df, gene_df, downstream_positions_dict, chrom_lengths, gene_ids_chunk = chunk_data
    chunk_exons, chunk_introns, chunk_dogs, chunk_genes, chunk_clusters = [], [], [], [], []
    discordant_pas = []
    cross_chromosome_clusters = defaultdict(int)
    genes_without_valid_clusters = []
    
    for gene_id in gene_ids_chunk:
        gene_info = gene_df[gene_df['gene_id'] == gene_id].iloc[0]
        gene_chrom = gene_info['chrom']
        gene_strand = gene_info['strand']
        
        gene_clusters = downstream_positions_dict.get(gene_id, [])
        
        # check for incompatible chromosomes
        cross_chrom_count = sum(1 for c in gene_clusters if c['chrom'] != gene_chrom)
        if cross_chrom_count > 0:
            cross_chromosome_clusters[gene_id] = cross_chrom_count
        
        # filter for correct chromosome and strand
        valid_clusters = [c for c in gene_clusters if c['chrom'] == gene_chrom and c['strand'] == gene_strand]
        
        # select most downstream valid cluster 
        cluster = None
        if valid_clusters:
            if gene_strand == '+':
                cluster = max(valid_clusters, key=lambda x: x['end'])
            else:  
                cluster = min(valid_clusters, key=lambda x: x['start'])
        
        exon, intron, dog, gene, cluster_data, error_message = adjust_regions_for_gene(gene_id, exon_df, intron_df, gene_df, cluster, chrom_lengths)
        chunk_exons.append(exon)
        chunk_introns.append(intron)
        chunk_dogs.append(dog)
        chunk_genes.append(gene)
        chunk_clusters.append(cluster_data)
        if error_message:
            discordant_pas.append(error_message)
        
        if cluster is None:
            genes_without_valid_clusters.append(gene_id)
    
    return (
        pd.concat(chunk_exons, ignore_index=True),
        pd.concat(chunk_introns, ignore_index=True),
        pd.concat(chunk_dogs, ignore_index=True),
        pd.concat(chunk_genes, ignore_index=True),
        pd.concat(chunk_clusters, ignore_index=True),
        discordant_pas,
        dict(cross_chromosome_clusters),
        genes_without_valid_clusters
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

    exon_df = read_bed(args.exon)
    intron_df = read_bed(args.intron)
    gene_df = read_bed(args.gene)
    clusters_df = read_clusters(args.clusters)
    chrom_lengths = read_fasta_index(args.fai)

    downstream_positions_dict = find_downstream_positions(clusters_df, gene_df)
    
    if not downstream_positions_dict:
        logging.warning("No valid downstream positions found. Outputting original annotations.")
        exon_df.to_csv(os.path.join(args.output_dir, 'updated_exon.bed'), sep='\t', header=False, index=False)
        intron_df.to_csv(os.path.join(args.output_dir, 'updated_intron.bed'), sep='\t', header=False, index=False)
        gene_df.to_csv(os.path.join(args.output_dir, 'updated_gene.bed'), sep='\t', header=False, index=False)
        pd.DataFrame(columns=['chrom', 'start', 'end', 'gene_id', 'score', 'strand']).to_csv(os.path.join(args.output_dir, 'updated_downstream_of_gene.bed'), sep='\t', header=False, index=False)
        pd.DataFrame(columns=['chrom', 'start', 'end', 'name', 'score', 'strand', 'distance']).to_csv(os.path.join(args.output_dir, 'PAS.bed'), sep='\t', header=False, index=False)
        return

    gene_ids = gene_df['gene_id'].unique()
    gene_chunks = np.array_split(gene_ids, args.cores)
    chunk_data = [(exon_df, intron_df, gene_df, downstream_positions_dict, chrom_lengths, gene_chunk) for gene_chunk in gene_chunks]

    with ProcessPoolExecutor(max_workers=args.cores) as executor:
        results = list(tqdm(executor.map(process_gene_chunk, chunk_data), total=len(gene_chunks), desc="Processing genes"))

    final_exons, final_introns, final_dogs, final_genes, final_clusters, all_discordant_pas, all_cross_chromosome, all_genes_without_valid_clusters = map(list, zip(*results))
  
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

    merged_cross_chromosome = {}
    for d in all_cross_chromosome:
        for gene, count in d.items():
            if gene in merged_cross_chromosome:
                merged_cross_chromosome[gene] += count
            else:
                merged_cross_chromosome[gene] = count

    os.makedirs(args.output_dir, exist_ok=True)

    final_exons.to_csv(os.path.join(args.output_dir, 'updated_exon.bed'), sep='\t', header=False, index=False, 
                       columns=['chrom', 'start', 'end', 'gene_id', 'combined', 'strand'])
    final_introns.to_csv(os.path.join(args.output_dir, 'updated_intron.bed'), sep='\t', header=False, index=False, 
                         columns=['chrom', 'start', 'end', 'gene_id', 'combined', 'strand'])
    final_dogs.to_csv(os.path.join(args.output_dir, 'updated_downstream_of_gene.bed'), sep='\t', header=False, index=False)
    final_genes.to_csv(os.path.join(args.output_dir, 'updated_gene.bed'), sep='\t', header=False, index=False)
    final_clusters.to_csv(os.path.join(args.output_dir, 'PAS.bed'), sep='\t', header=False, index=False)

    log_file = os.path.join(args.output_dir, 'processing_log.txt')
    with open(log_file, 'w') as f:
        f.write(f"Total discordant PAS: {len(all_discordant_pas)}\n\n")
        for pas in all_discordant_pas:
            f.write(f"{pas}\n")
        
        f.write(f"\nTotal genes with cross-chromosome clusters: {len(merged_cross_chromosome)}\n")
        for gene, count in merged_cross_chromosome.items():
            f.write(f"Gene {gene}: {count} cluster(s) on different chromosome(s)\n")
        
        f.write(f"\nTotal genes without valid clusters: {len(all_genes_without_valid_clusters)}\n")
        for gene in all_genes_without_valid_clusters:
            f.write(f"Gene without valid clusters: {gene}\n")

    logging.info(f"All updated BED files have been created in the specified output directory.")
    logging.info(f"Found {len(all_discordant_pas)} genes with discordant PAS.")
    logging.info(f"Found {len(merged_cross_chromosome)} genes with cross-chromosome clusters.")
    logging.info(f"Found {len(all_genes_without_valid_clusters)} genes without valid clusters.")
    logging.info(f"Details logged in {log_file}")

if __name__ == "__main__":
    main()