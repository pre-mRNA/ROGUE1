import logging
from typing import Dict
import pandas as pd 
import re 

# subfunction for parsing attributes from a gtf file 
def parse_attributes(attribute_field):

    """Parse the GTF attribute field into a dictionary of key-value pairs."""
    attribute_dict = {}
    attribute_list = attribute_field.split(';')

    for attribute in attribute_list:
        if attribute:
            key_value = attribute.strip().split(' ')
            if len(key_value) == 2:
                attribute_dict[key_value[0]] = re.sub('"', '', key_value[1])
    return attribute_dict

# function to parse GTF line by line - needed for get_gene_exon_table 
def parse_gtf_line(line: str) -> Dict[str, str]:
    """Parse a single line in GTF file."""
    try:
        fields = line.strip().split('\t')
        attributes_field = fields[8]
        attributes = {k: v.strip('"') for k, v in [attr.strip().split(' ') for attr in attributes_field.split(';') if attr.strip()]}
        return attributes
    except IndexError as e:
        logging.error(f"Error parsing line: {line}. Error: {e}")
        return {}
    


# function to return the minimum number of annotated exons per transcript from gtf - requires parse_gtf_line
def get_gene_exon_table(gtf_file):
    logging.info("Starting to parse GTF file.")
    # Read GTF into DataFrame
    df = pd.read_csv(gtf_file, comment='#', sep='\t', header=None, 
                     names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])

    # Filter only exons
    df = df[df['feature'] == 'exon']
    
    # If df is empty, early exit
    if df.empty:
        logging.warning("No exon entries found in GTF file.")
        return df

    # Extract gene information
    df['gene_id'] = df['attribute'].str.extract('gene_id "([^"]+)"')
    df['gene_name'] = df['attribute'].str.extract('gene_name "([^"]+)"')
    
    # Check for NaN values in gene_id and type mismatch
    if df['gene_id'].isna().any():
        logging.warning("Found NaN values in gene_id. Investigate the GTF file.")
    
    df['gene_id'] = df['gene_id'].astype(str)
    
    logging.info("GTF file parsed. Processing data.")

    # Count the number of exons for each gene
    exon_counts = df.groupby('gene_id').size().reset_index(name='exon_count')
    # Find minimum exon_count per gene
    min_exon_counts = exon_counts.groupby('gene_id')['exon_count'].min().reset_index()
    
    # Merge with original DataFrame to get gene names
    final_df = pd.merge(df.drop_duplicates(['gene_id']), min_exon_counts, on='gene_id', how='left')
    final_df = final_df[['gene_name', 'gene_id', 'exon_count']]

    logging.info("Data processed successfully. Returning final DataFrame.")
    

    return final_df

# function to get gene biotypes 
def get_biotypes(gtf_file):

    logging.info("Getting biotypes from GTF file")

    # specify the column names for the GTF file
    column_names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

    # load the GTF file
    gtf_df = pd.read_csv(gtf_file, sep="\t", comment='#', names=column_names)

    # filter rows where feature is 'gene'
    gtf_df = gtf_df[gtf_df['feature'] == 'gene']

    # parse the attribute column into a dictionary
    gtf_df['attributes'] = gtf_df['attribute'].apply(parse_attributes)

    # create new columns for gene_id and gene_biotype
    gtf_df['gene_id'] = gtf_df['attributes'].apply(lambda x: x.get('gene_id'))
    gtf_df['gene_biotype'] = gtf_df['attributes'].apply(lambda x: x.get('gene_biotype'))

    # print the head of the DataFrame
    print("head of the gtf DataFrame:\n", gtf_df.head())

    # print summary stats
    print("\nsummary stats:")
    print("number of unique genes: ", gtf_df['gene_id'].nunique())
    print("number of unique gene biotypes: ", gtf_df['gene_biotype'].nunique())

    # create a dictionary where keys are gene_ids and values are gene_biotypes
    gene_dict = dict(zip(gtf_df['gene_id'], gtf_df['gene_biotype']))

    # Get the gene-exon table for additional info
    gene_exon_table = get_gene_exon_table(gtf_file)
    gene_exon_dict = gene_exon_table.set_index('gene_id').T.to_dict('list')
    
    # Combine biotype, gene_name, and exon_count into a single dictionary
    for gene_id, biotype in gene_dict.items():
        gene_name, exon_count = gene_exon_dict.get(gene_id, [None, None])
        gene_dict[gene_id] = {'biotype': biotype, 'gene_name': gene_name, 'exon_count': exon_count}
    
    return gene_dict


