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
        else:  # minus strand
            valid_clusters = valid_clusters[valid_clusters['end'] <= gene_end]
        
        if not valid_clusters.empty:
            downstream_positions[gene_id] = valid_clusters.to_dict('records')
    
    return downstream_positions

def number_exons_and_introns(exons, introns, gene_id, strand):
    regions = pd.concat([exons, introns])
    regions['region_type'] = ['exon'] * len(exons) + ['intron'] * len(introns)
    
    if strand == '+':
        regions = regions.sort_values('start')
    else:
        regions = regions.sort_values('start', ascending=False)
    
    exon_count = 1
    intron_count = 1
    for idx, row in regions.iterrows():
        if row['region_type'] == 'exon':
            regions.at[idx, 'name'] = f"{gene_id}_exon_{exon_count}"
            exon_count += 1
        else:
            regions.at[idx, 'name'] = f"{gene_id}_intron_{intron_count}"
            intron_count += 1
    
    numbered_exons = regions[regions['region_type'] == 'exon'].drop('region_type', axis=1)
    numbered_introns = regions[regions['region_type'] == 'intron'].drop('region_type', axis=1)
    
    return numbered_exons, numbered_introns

def adjust_regions_for_gene(gene_id, exon_df, intron_df, dog_df, gene_df, cluster, window=5):
    if cluster is None:
        return (
            exon_df[exon_df['gene_id'] == gene_id],
            intron_df[intron_df['gene_id'] == gene_id],
            dog_df[dog_df['gene_id'] == gene_id],
            gene_df[gene_df['gene_id'] == gene_id],
            pd.DataFrame(),
            None
        )
    
    gene_exons = exon_df[exon_df['gene_id'] == gene_id].sort_values('start')
    gene_introns = intron_df[intron_df['gene_id'] == gene_id].sort_values('start')
    gene_dog = dog_df[dog_df['gene_id'] == gene_id]
    gene_entry = gene_df[gene_df['gene_id'] == gene_id].iloc[0]
    
    pas_position = cluster['end'] if cluster['strand'] == '+' else cluster['start']
    
    if (cluster['strand'] == '+' and pas_position < gene_entry['start']) or \
       (cluster['strand'] == '-' and pas_position > gene_entry['end']):
        return (
            gene_exons,
            gene_introns,
            gene_dog,
            pd.DataFrame([gene_entry]),
            pd.DataFrame(),
            f"PAS position {pas_position} is {'before' if cluster['strand'] == '+' else 'after'} gene {'start' if cluster['strand'] == '+' else 'end'} {gene_entry['start' if cluster['strand'] == '+' else 'end']} for gene {gene_id}"
        )

    if cluster['strand'] == '+':
        containing_exon = gene_exons[(gene_exons['start'] <= pas_position) & (gene_exons['end'] > pas_position)]
        containing_intron = gene_introns[(gene_introns['start'] <= pas_position) & (gene_introns['end'] > pas_position)]
        
        if not containing_exon.empty:
            last_exon = containing_exon.iloc[0].copy()
            last_exon['end'] = pas_position
            exons_to_keep = gene_exons[gene_exons['end'] <= pas_position]
            exons_to_keep = pd.concat([exons_to_keep, pd.DataFrame([last_exon])])
        elif not containing_intron.empty:
            last_exon = gene_exons[gene_exons['end'] <= pas_position].iloc[-1].copy()
            last_exon['end'] = pas_position
            exons_to_keep = gene_exons[gene_exons['end'] < pas_position]
            exons_to_keep = pd.concat([exons_to_keep, pd.DataFrame([last_exon])])
        else:
            exons_to_keep = gene_exons
        
        introns_to_keep = gene_introns[gene_introns['end'] < pas_position]
        
        new_dog_start = pas_position + window
        new_dog_end = new_dog_start + (gene_dog['end'].iloc[0] - gene_dog['start'].iloc[0]) if not gene_dog.empty else new_dog_start + window
        
        new_gene_entry = gene_entry.copy()
        new_gene_entry['end'] = pas_position
        
    else:  # minus strand
        containing_exon = gene_exons[(gene_exons['start'] < pas_position) & (gene_exons['end'] >= pas_position)]
        containing_intron = gene_introns[(gene_introns['start'] < pas_position) & (gene_introns['end'] >= pas_position)]
        
        if not containing_exon.empty:
            first_exon = containing_exon.iloc[0].copy()
            first_exon['start'] = pas_position
            exons_to_keep = gene_exons[gene_exons['start'] >= pas_position]
            exons_to_keep = pd.concat([pd.DataFrame([first_exon]), exons_to_keep])
        elif not containing_intron.empty:
            first_exon = gene_exons[gene_exons['start'] >= pas_position].iloc[0].copy()
            first_exon['start'] = pas_position
            exons_to_keep = gene_exons[gene_exons['start'] > pas_position]
            exons_to_keep = pd.concat([pd.DataFrame([first_exon]), exons_to_keep])
        else:
            exons_to_keep = gene_exons
        
        introns_to_keep = gene_introns[gene_introns['start'] > pas_position]
        
        new_dog_end = pas_position - window
        new_dog_start = new_dog_end - (gene_dog['end'].iloc[0] - gene_dog['start'].iloc[0]) if not gene_dog.empty else new_dog_end - window
        
        new_gene_entry = gene_entry.copy()
        new_gene_entry['start'] = pas_position

    new_dog = pd.DataFrame({
        'chrom': cluster['chrom'],
        'start': new_dog_start,
        'end': new_dog_end,
        'gene_id': gene_id,
        'score': gene_dog['score'].iloc[0] if not gene_dog.empty else '.',
        'strand': cluster['strand']
    }, index=[0])
    
    exons_to_keep, introns_to_keep = number_exons_and_introns(exons_to_keep, introns_to_keep, gene_id, cluster['strand'] if cluster else gene_entry['strand'])

    return exons_to_keep, introns_to_keep, new_dog, pd.DataFrame([new_gene_entry]), pd.DataFrame([cluster]) if cluster else pd.DataFrame(), None

def process_gene_chunk(chunk_data):
    exon_df, intron_df, dog_df, gene_df, downstream_positions_dict, window, gene_ids_chunk = chunk_data
    chunk_exons, chunk_introns, chunk_dogs, chunk_genes, chunk_clusters = [], [], [], [], []
    discordant_pas = []
    cross_chromosome_clusters = defaultdict(int)
    genes_without_valid_clusters = []
    
    for gene_id in gene_ids_chunk:
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
            if gene_strand == '+':
                cluster = max(valid_clusters, key=lambda x: x['end'])
            else:  # minus strand
                cluster = min(valid_clusters, key=lambda x: x['start'])
        
        if cluster is None:
            
            # return original gene structure if no valid clusters
            genes_without_valid_clusters.append(gene_id)
            chunk_exons.append(exon_df[exon_df['gene_id'] == gene_id])
            chunk_introns.append(intron_df[intron_df['gene_id'] == gene_id])
            chunk_dogs.append(dog_df[dog_df['gene_id'] == gene_id])
            chunk_genes.append(pd.DataFrame([gene_info]))
            chunk_clusters.append(pd.DataFrame())
        else:
            exon, intron, dog, gene, cluster_data, error_message = adjust_regions_for_gene(gene_id, exon_df, intron_df, dog_df, gene_df, cluster, window)
            chunk_exons.append(exon)
            chunk_introns.append(intron)
            chunk_dogs.append(dog)
            chunk_genes.append(gene)
            chunk_clusters.append(cluster_data)
            if error_message:
                discordant_pas.append(error_message)
    
    return (
        pd.concat(chunk_exons),
        pd.concat(chunk_introns),
        pd.concat(chunk_dogs),
        pd.concat(chunk_genes),
        pd.concat(chunk_clusters),
        discordant_pas,
        dict(cross_chromosome_clusters),
        genes_without_valid_clusters
    )

def main():
    parser = argparse.ArgumentParser(description="Process genomic intervals and create updated BED files.")
    parser.add_argument("--gene", required=True, help="Path to gene BED file")
    parser.add_argument("--exon", required=True, help="Path to exon BED file")
    parser.add_argument("--intron", required=True, help="Path to intron BED file")
    parser.add_argument("--dog", required=True, help="Path to downstream of gene BED file")
    parser.add_argument("--clusters", required=True, help="Path to polyA clusters BED file")
    parser.add_argument("--output_dir", required=True, help="Directory to save output files")
    parser.add_argument("--window", type=int, default=5, help="Window size for adjusting regions")
    parser.add_argument("--cores", type=int, default=os.cpu_count(), help="Number of CPU cores to use")
    args = parser.parse_args()

    exon_df = read_bed(args.exon)
    intron_df = read_bed(args.intron)
    gene_df = read_bed(args.gene)
    dog_df = read_bed(args.dog)
    clusters_df = read_clusters(args.clusters)

    downstream_positions_dict = find_downstream_positions(clusters_df, gene_df)
    
    if not downstream_positions_dict:
        logging.warning("No valid downstream positions found. Outputting original annotations.")
        exon_df.to_csv(os.path.join(args.output_dir, 'updated_exon.bed'), sep='\t', header=False, index=False)
        intron_df.to_csv(os.path.join(args.output_dir, 'updated_intron.bed'), sep='\t', header=False, index=False)
        dog_df.to_csv(os.path.join(args.output_dir, 'updated_downstream_of_gene.bed'), sep='\t', header=False, index=False)
        gene_df.to_csv(os.path.join(args.output_dir, 'updated_gene.bed'), sep='\t', header=False, index=False)
        pd.DataFrame(columns=['chrom', 'start', 'end', 'name', 'score', 'strand', 'distance']).to_csv(os.path.join(args.output_dir, 'PAS.bed'), sep='\t', header=False, index=False)
        return

    gene_ids = gene_df['gene_id'].unique()
    gene_chunks = np.array_split(gene_ids, args.cores)
    chunk_data = [(exon_df, intron_df, dog_df, gene_df, downstream_positions_dict, args.window, gene_chunk) for gene_chunk in gene_chunks]

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

    # merge cross-chromosome dictionaries
    merged_cross_chromosome = {}
    for d in all_cross_chromosome:
        for gene, count in d.items():
            if gene in merged_cross_chromosome:
                merged_cross_chromosome[gene] += count
            else:
                merged_cross_chromosome[gene] = count


    os.makedirs(args.output_dir, exist_ok=True)

    final_exons.to_csv(os.path.join(args.output_dir, 'updated_exon.bed'), sep='\t', header=False, index=False)
    final_introns.to_csv(os.path.join(args.output_dir, 'updated_intron.bed'), sep='\t', header=False, index=False)
    final_dogs.to_csv(os.path.join(args.output_dir, 'updated_downstream_of_gene.bed'), sep='\t', header=False, index=False)
    final_genes.to_csv(os.path.join(args.output_dir, 'updated_gene.bed'), sep='\t', header=False, index=False)
    final_clusters.to_csv(os.path.join(args.output_dir, 'PAS.bed'), sep='\t', header=False, index=False)

    # log discordant PAS, cross-chromosome clusters, and genes without valid clusters
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