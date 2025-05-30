import os
import tempfile 
from tempfile import NamedTemporaryFile
import subprocess
import pandas as pd
import logging 
from run_command import run_command

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

# convert genes or exons to bed tracks 
def gtf_to_bed(gtf_file, feature_type, output_dir):

    bed_file = NamedTemporaryFile(dir=output_dir, delete=False, suffix=f"_{feature_type}_gtf_to_bed.bed").name

    if feature_type == "gene":

        # infer gene boundaries from exons
        logging.info("Extracting exons to infer gene boundaries")
        exon_bed = gtf_to_bed(gtf_file, "exon", output_dir)

        try:
            exon_df = pd.read_csv(exon_bed, sep='\t', header=None, names=['chromosome', 'start', 'end', 'gene_id', 'gene_name', 'strand'])
        except Exception as e:
            logging.error(f"Failed to read exon BED file: {exon_bed}. Error: {e}")
            raise e

        logging.info(f"Number of exons extracted: {len(exon_df)}")

        # group exons by gene and calculate outer boundaries 
        try:
            gene_df = exon_df.groupby('gene_id').agg({
                'chromosome': 'first',
                'start': 'min',
                'end': 'max',
                'gene_name': 'first',
                'strand': 'first'
            }).reset_index()

            gene_df = gene_df[['chromosome', 'start', 'end', 'gene_id', 'gene_name', 'strand']]

            # Save gene_df to a temporary file
            with tempfile.NamedTemporaryFile(mode='w+t', delete=False, suffix='_gene_unsorted.bed') as temp_gene_file:
                gene_df.to_csv(temp_gene_file.name, sep='\t', header=False, index=False)
                temp_gene_path = temp_gene_file.name
            logging.info(f"Unsorted Gene BED file saved to temporary file: {temp_gene_path}")

            # Sort the temporary gene file and save as bed_file
            sort_bed_file(temp_gene_path, bed_file)
            logging.info(f"Sorted Gene BED file created successfully: {bed_file}")

            # Clean up the temporary unsorted gene file
            os.unlink(temp_gene_path)
            logging.info(f"Temporary unsorted Gene BED file deleted: {temp_gene_path}")

        except Exception as e:
            logging.error(f"Failed to process gene boundaries. Error: {e}")
            raise e

        
    elif feature_type == "exon":
        # we use a complex set of commands to merge overlapping exon entries by gene, otherwise, the overlap lengths by gene(at exon level) are confounded by overlapping exons
        # inspired by https://www.biostars.org/p/469489/#469503
        cmd = f"cat {gtf_file} |  grep -v \"intron\" | grep -v \"protein_coding_CDS_not_defined\" | awk 'OFS=\"\\t\" {{if ($3==\"exon\") {{print $1,$4-1,$5,$10,$20,$7}}}}' | tr -d '\";' | awk 'OFS=\"\\t\" {{print $4, $5, $6, $2, $3, $1}}' | sed 's/\\t/;/' | sed 's/\\t/;/' | sort -k1,1 -k2,2n --parallel=104 --buffer-size=80G | bedtools merge -i - -c 4 -o first | sed 's/;/\\t/' | sed 's/;/\\t/' | awk 'OFS=\"\\t\" {{print $6, $4, $5, $1, $2, $3}}' | sort -k1,1 -k2,2n -k3,3n --parallel=48 --buffer-size=32G > {bed_file}"
        logging.info("Skipping intron-contianing and protein-coding CDS not defined biotypes while parsing exons")
        run_command(cmd, "extracting exons ")

    elif feature_type == "transcript":
        cmd = f"cat {gtf_file} | awk 'OFS=\"\\t\" {{if ($3==\"transcript\") {{print $1,$4-1,$5,$10,$18,$7}}}}' | tr -d '\";' | awk 'OFS=\"\\t\" {{if ($6==\"+\") {{print $1,$3,$3,$4,$5,$6}} else {{print $1,$2,$2,$4,$5,$6}}}}' | sort -k1,1 -k2,2n --parallel=48 --buffer-size=32G > {bed_file}"
        run_command(cmd, "Converting intron by subtracting exons from trimmed genes")

    elif feature_type == "intron":
        gene_bed = gtf_to_bed(gtf_file, "gene", output_dir)
        exon_bed = gtf_to_bed(gtf_file, "exon", output_dir)
        cmd = f"""
        bedtools subtract -s -a <(awk 'BEGIN{{FS=OFS="\\t"}} {{print $1 ";" $4, $2, $3, $4, $5, $6}}' {gene_bed} | sort -k1,1 -k2,2n) \
        -b <(awk 'BEGIN{{FS=OFS="\\t"}} {{print $1 ";" $4, $2, $3, $4, $5, $6}}' {exon_bed} | sort -k1,1 -k2,2n) | \
        awk 'BEGIN{{FS=OFS="\\t"}} {{split($1, a, ";"); $1=a[1]; print}}' > {bed_file}
        """
        run_command(cmd, "Converting intron by subtracting exons from trimmed genes")

    else:
            raise ValueError(f"Unknown feature type: {feature_type}. Valid options are 'gene', 'exon', and 'transcript'.")

    return bed_file

# number exons and introns by strand 
def number_exons_and_introns_in_bed(input_bed, output_bed, feature_type):
    df = pd.read_csv(input_bed, sep='\t', header=None, 
                     names=['chrom', 'start', 'end', 'gene_id', 'score', 'strand'],
                     low_memory=False)
    
    # number features grouped by gene and considering strand 
    def number_features(group):
        if group['strand'].iloc[0] == '+':
            group = group.sort_values('start')
        else:
            group = group.sort_values('start', ascending=False)
        
        group['number'] = range(1, len(group) + 1)
        return group

    df = df.groupby('gene_id').apply(number_features).reset_index(drop=True)
    
    # create exon/intron labels 
    df['combined'] = df['gene_id'] + '_' + feature_type + '_' + df['number'].astype(str)
    
    # write to a temporary file and sort 
    with tempfile.NamedTemporaryFile(mode='w+t', delete=False, suffix='.bed') as temp_file:
        df.to_csv(temp_file.name, sep='\t', header=False, index=False, 
                  columns=['chrom', 'start', 'end', 'gene_id', 'combined', 'strand'])
    
    sort_bed_file(temp_file.name, output_bed)
    os.unlink(temp_file.name)
    
# extends gene region ends to create downstream-of-gene regions 
def extend_gene_bed(gene_bed_file, output_dir, genome_file, extend_bases=500):
    extended_bed = NamedTemporaryFile(dir=output_dir, delete=False, suffix="_extendedGene.bed").name

    genome = pd.read_csv(genome_file, sep="\t", header=None, index_col=0).squeeze().to_dict() if pd.read_csv(genome_file, sep="\t", header=None, index_col=0).shape[1] == 1 else pd.read_csv(genome_file, sep="\t", header=None, index_col=0).to_dict()
    # use gene bed file 
    df = pd.read_csv(gene_bed_file, sep="\t", header=None)

    # extend depending on strand 
    df[[1, 2]] = df.apply(
        lambda x: [x[2], min(x[2] + extend_bases, genome.get(x[0], float('inf')))] if x[5] == "+" 
        else [max(x[1] - extend_bases, 0), x[1]], axis=1, result_type='expand')

    # check if there are chromosomes not found in the genome file
    # missing_chromosomes = set(df[0].unique()) - set(genome.keys())
    # if missing_chromosomes:
    #     print(f"Warning: The following chromosomes were not found in the genome file: {', '.join(missing_chromosomes)}")

    # debug 
    # print(df.head())  # check the first few entries to see how they are modified
    # print(genome)  # print the genome dictionary to ensure correct chromosome lengths

    # add region label 
    df[4] = df[3] + "_region_DOG"
    
    # sort equivalent to -k1,1 -k2,2n; type convert to enable this
    df[0] = df[0].astype(str)
    df[1] = df[1].astype(int)
    df = df.sort_values(by=[0, 1])

    # save
    df.to_csv(extended_bed, sep="\t", header=False, index=False)

    return extended_bed

def sort_bed_file(input_file, output_file):
    cmd = f"sort -k1,1 -k2,2n --parallel=104 --buffer-size=80G {input_file} > {output_file}"
    subprocess.run(cmd, shell=True, check=True)
    return output_file