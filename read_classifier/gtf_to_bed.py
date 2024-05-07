import os
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
        cmd = f"cat {gtf_file} | awk 'OFS=\"\\t\" {{if ($3==\"gene\") {{print $1,$4-1,$5,$10,$20,$7}}}}' | tr -d '\";' | sort -k1,1 -k2,2n --parallel=104 --buffer-size=80G > {bed_file}"
    
    elif feature_type == "exon":
        # we use a complex set of commands to merge overlapping exon entries by gene, otherwise, the overlap lengths by gene(at exon level) are confounded by overlapping exons
        # inspired by https://www.biostars.org/p/469489/#469503
        cmd = f"cat {gtf_file} |  grep -v \"intron\" | grep -v \"protein_coding_CDS_not_defined\" | awk 'OFS=\"\\t\" {{if ($3==\"exon\") {{print $1,$4-1,$5,$10,$20,$7}}}}' | tr -d '\";' | awk 'OFS=\"\\t\" {{print $4, $5, $6, $2, $3, $1}}' | sed 's/\\t/;/' | sed 's/\\t/;/' | sort -k1,1 -k2,2n --parallel=104 --buffer-size=80G | bedtools merge -i - -c 4 -o first | sed 's/;/\\t/' | sed 's/;/\\t/' | awk 'OFS=\"\\t\" {{print $6, $4, $5, $1, $2, $3}}' | sort -k1,1 -k2,2n -k3,3n --parallel=104 --buffer-size=80G > {bed_file}"
        logging.info("Skipping intron-contianing and protein-coding CDS not defined biotypes while parsing exons")

    elif feature_type == "transcript":
        cmd = f"cat {gtf_file} | awk 'OFS=\"\\t\" {{if ($3==\"transcript\") {{print $1,$4-1,$5,$10,$18,$7}}}}' | tr -d '\";' | awk 'OFS=\"\\t\" {{if ($6==\"+\") {{print $1,$3,$3,$4,$5,$6}} else {{print $1,$2,$2,$4,$5,$6}}}}' | sort -k1,1 -k2,2n --parallel=104 --buffer-size=80G > {bed_file}"

    else:
            raise ValueError(f"Unknown feature type: {feature_type}. Valid options are 'gene', 'exon', and 'transcript'.")

    run_command(cmd, f"Converting {feature_type} GTF to BED")
    return bed_file

    # takes as input the bedfile of genes and extends them to DOG regions by a fixed number of bases 
def extend_gene_bed(gene_bed_file, output_dir, genome_file, extend_bases=500):
    extended_bed = NamedTemporaryFile(dir=output_dir, delete=False, suffix="_extendedGene.bed").name

    # Read the genome file into a dictionary
    genome = pd.read_csv(genome_file, sep="\t", header=None, index_col=0, squeeze=True).to_dict()

    df = pd.read_csv(gene_bed_file, sep="\t", header=None)

    # Extend the end coordinates and cap them at the corresponding genome length
    # df[[1,2]] = df.apply(lambda x: [x[2], min(max(0, x[2] + extend_bases), genome.get(x[0], x[2]))] if x[5]=="+" else [max(0, min(x[1]-extend_bases, genome.get(x[0], x[1]))), x[1]], axis=1, result_type='expand')
    # Extend the end coordinates based on the strand
    # Extend the gene coordinates based on the strand
    # Extend the gene coordinates based on the strand
    df[[1, 2]] = df.apply(
        lambda x: [x[2], min(x[2] + extend_bases, genome.get(x[0], float('inf')))] if x[5] == "+" 
        else [max(x[1] - extend_bases, 0), x[1]], axis=1, result_type='expand')

    # Check if there are chromosomes not found in the genome file
    # missing_chromosomes = set(df[0].unique()) - set(genome.keys())
    # if missing_chromosomes:
    #     print(f"Warning: The following chromosomes were not found in the genome file: {', '.join(missing_chromosomes)}")

    print(df.head())  # Check the first few entries to see how they are modified
    print(genome)  # Print the genome dictionary to ensure correct chromosome lengths

    # Sort the DataFrame by the first and second column
    df.sort_values([0, 1], inplace=True)

    df.to_csv(extended_bed, sep="\t", header=False, index=False)

    return extended_bed