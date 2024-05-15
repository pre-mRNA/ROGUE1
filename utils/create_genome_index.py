import argparse
import os
import subprocess
import pandas as pd
import logging 
import sys
from tempfile import NamedTemporaryFile

# import gtf_to_bed
sys.path.append("/home/150/as7425/R1/read_classifier")
from gtf_to_bed import gtf_to_bed, extend_gene_bed

def setup_logging(level):
    logging.basicConfig(format='%(asctime)s - %(levelname)s: %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=level)

# run shell commands in python 
def run_command(command):
    try:
        subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {e.cmd}")
        logging.error(f"Error message: {e.stderr.decode()}")
        raise

# we don't use a bam when creating the genome index, so we can't use the bam header to create a genome file 
# instead, use the fasta index 
# see: https://www.biostars.org/p/70795/#211287
def create_genome_file_from_fasta_index(fasta_index_file, output_dir):
    """Creates a genome size file from a FASTA index."""
    genome_file_path = os.path.join(output_dir, "genome_sizes.txt")
    with open(fasta_index_file, 'r') as f_index, open(genome_file_path, 'w') as f_out:
        for line in f_index:
            parts = line.strip().split('\t')
            chrom = parts[0]
            size = parts[1]
            f_out.write(f"{chrom}\t{size}\n")
    return genome_file_path

# use bedtools to create the introns from collapsed exons and gene regions 
def create_introns(gene_bed, exon_bed, output_dir, verbose):
    intron_bed = os.path.join(output_dir, "introns.bed")
    command = f"bedtools subtract -s -a {gene_bed} -b {exon_bed} > {intron_bed}"
    run_command(command)
    if verbose:
        logging.info(f"Introns created at {intron_bed}, preview:")
        print(pd.read_csv(intron_bed, sep="\t", header=None).head())
    return intron_bed

# ensure that dog regions don't overlap with a downstream gene
def create_non_intersecting_dogs(gene_bed, extended_gene_bed, output_dir, verbose):
    dog_bed = os.path.join(output_dir, "dogs.bed")
    command = f"bedtools subtract -a {extended_gene_bed} -b {gene_bed} > {dog_bed}"
    run_command(command)
    if verbose:
        logging.info(f"DOG regions created at {dog_bed}, preview:")
        print(pd.read_csv(dog_bed, sep="\t", header=None).head())
    return dog_bed

# merge exon, intron and clipped dog features into a single bed file (our index)
def merge_features(output_dir, exon_ranges, intron_ranges, dog_ranges, verbose):


    exons_bed = exon_ranges
    introns_bed = intron_ranges
    dogs_bed = dog_ranges
    merged_bed = os.path.join(output_dir, "all_features.bed")

    # Prepare the command to adjust columns and set feature type
    cmd = f"awk 'BEGIN{{OFS=\"\\t\"}} {{$5=\"exon\"; $6=$6; $7=$4; $4=\".\"; print}}' {exons_bed} > {output_dir}/tmp_exons.bed; "
    cmd += f"awk 'BEGIN{{OFS=\"\\t\"}} {{$6=$5; $5=\"intron\"; $7=$4; $4=\".\"; print}}' {introns_bed} > {output_dir}/tmp_introns.bed; "
    cmd += f"awk 'BEGIN{{OFS=\"\\t\"}} {{$6=$5; $5=\"DOG\"; $7=$4; $4=\".\"; print}}' {dogs_bed} > {output_dir}/tmp_dogs.bed; "
    cmd += f"cat {output_dir}/tmp_exons.bed {output_dir}/tmp_introns.bed {output_dir}/tmp_dogs.bed | "
    # cmd += "sort -k1,1 -k2,2n | bedtools merge -i - -c 5,7 -o distinct -s > " + merged_bed
    cmd += "sort -k1,1 -k2,2n > " + merged_bed
    run_command(cmd)

    if verbose:
        logging.info(f"Merged features file created at {merged_bed}, preview:")
        print(pd.read_csv(merged_bed, sep="\t", header=None).head())

    return merged_bed

def merge_feature(bed_file, output_dir, feature):
    merged_file_path = os.path.join(output_dir, f"merged_{feature}.bed")
    
    # debug 
    # print(f"feature is {feature} and output dir is {output_dir} and path is {merged_file_path}")
    
    command = f"bedtools merge -i {bed_file} -s -c 4,5,6 -o 'distinct' > {merged_file_path}"
    subprocess.run(command, shell=True, check=True)
    
    # overwrite the original BED file with the merged results
    subprocess.run(f"mv {merged_file_path} {bed_file}", shell=True, check=True)
    logging.info(f"Merged exon BED file created and original overwritten: {bed_file}")



def color_features(input_bed, output_bed):

    df = pd.read_csv(input_bed, sep='\t', header=None, names=['chr', 'start', 'end', 'dot', 'feature', 'strand', 'gene_id'])

    colors = {
        'DOG': '255,0,0',   # red
        'exon': '0,255,0',  # green
        'intron': '0,0,255' # blue
    }
    
    df['color'] = df['feature'].map(colors)

    # blockcount for visualisation 
    df['blockCount'] = '.'

    # fit extended bed format 
    df = df[['chr', 'start', 'end', 'dot', 'feature', 'strand', 'gene_id', 'blockCount', 'color']]

    df.to_csv(output_bed, sep='\t', header=False, index=False)



def annotate_upstream_regions(gene_bed_file, output_dir, genome_file, extend_bases=500):
    upstream_bed = NamedTemporaryFile(dir=output_dir, delete=False, suffix="_upstreamRegions.bed").name

    # Load the genome sizes to ensure extensions do not exceed chromosome lengths
    genome = pd.read_csv(genome_file, sep="\t", header=None, index_col=0, squeeze=True).to_dict()

    # Read the gene BED file
    df = pd.read_csv(gene_bed_file, sep="\t", header=None)

    # Extend gene coordinates based on strand
    def extend_upstream(row):
        chrom, start, end, strand = row[0], row[1], row[2], row[5]
        if strand == "+":
            new_start = max(start - extend_bases, 0)  # Cap at 0 to avoid negative start coordinates
            return pd.Series([chrom, new_start, start, *row[3:5], strand, *row[6:]])
        else:
            new_end = min(end + extend_bases, genome.get(chrom, float('inf')))
            return pd.Series([chrom, end, new_end, *row[3:5], strand, *row[6:]])

    # Apply function to each row
    df = df.apply(extend_upstream, axis=1)

    # Print the first few entries to see how they are modified
    print(df.head())
    # Print the genome dictionary to ensure correct chromosome lengths
    print(genome)

    # Sort DataFrame by chromosome and start position for consistent order
    df.sort_values([0, 1], inplace=True)

    # Write the modified DataFrame to a file
    df.to_csv(upstream_bed, sep="\t", header=False, index=False)

    logging.info(f"Upstream regions annotated at: {upstream_bed}")
    return upstream_bed

def duplicate_to_minus_strand(input_bed, output_dir):
    output_bed = os.path.join(output_dir, "all_features_with_antisense.bed")
    # TODO: verify this command
    command = f"awk 'BEGIN{{OFS=\"\\t\"}} {{print $1, $2, $3, $4, $5, ($6 == \"+\" ? \"-\" : \"+\"), $7, $8, $9; print $0}}' {input_bed} > {output_bed}"
    run_command(command)
    logging.info(f"Features duplicated to opposite strand at: {output_bed}")
    return output_bed



def clip_features(base_bed, features_to_remove, output_dir, feature_type):
   
    clipped_bed = os.path.join(output_dir, f"{feature_type}_clipped_features.bed")
    
    command = f"bedtools subtract -a {base_bed} -b {features_to_remove} > {clipped_bed}"
    run_command(command)
    
    logging.info(f"Features clipped for {feature_type} at: {clipped_bed}")
    return clipped_bed

def merge_all_features(feature_beds, output_dir):
    merged_bed = os.path.join(output_dir, "merged_all_features.bed")
    bed_list = " ".join(feature_beds)
    command = f"cat {bed_list} | sort -k1,1 -k2,2n | bedtools merge -i - > {merged_bed}"
    run_command(command)
    logging.info(f"All features merged into: {merged_bed}")
    return merged_bed



def main(args):

    setup_logging(logging.DEBUG if args.verbose else logging.INFO)

    output_dir = os.path.dirname(args.output_bed) 
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logging.info("Output directory created.")
    
    # create genome index for bedtools from the FASTA index file 
    genome_file = create_genome_file_from_fasta_index(args.fasta_index_file, output_dir)
    logging.info(f"Genome size file created at: {genome_file}")

    # extract gene and exon features from gtf 
    feature_types = ['gene', 'exon']
    bed_files = {}

    for feature_type in feature_types:
        bed_file = gtf_to_bed(args.gtf_file, feature_type, output_dir)
        bed_files[feature_type] = bed_file
        logging.info(f"BED file for {feature_type} created: {bed_file}")

        # we still have overlapping exons in cases where two exons share an exon region 
        # next, merge exons 
        if feature_type == 'exon':
            merge_feature(bed_file, output_dir, "exons")

        if feature_type == 'gene':
            merge_feature(bed_file, output_dir, "gene")


    # create downstream of gene regions 
    if 'gene' in bed_files:
        extended_bed_file = extend_gene_bed(bed_files['gene'], output_dir, genome_file, args.extend_bases)
        logging.info(f"Extended gene regions saved at: {extended_bed_file}")
        merge_feature(extended_bed_file, output_dir, "dog")
        

    # create introns by subtracting exons from genes
    if 'gene' in bed_files and 'exon' in bed_files:
        intron_bed = create_introns(bed_files['gene'], bed_files['exon'], output_dir, args.verbose)

    # create non-intersecting DOGs
    if 'gene' in bed_files:
        dog_bed = create_non_intersecting_dogs(bed_files['gene'], extended_bed_file, output_dir, args.verbose)

    # merge all BED files for introns, exons, dogs 
    merge_features(output_dir, bed_files['exon'], intron_bed, dog_bed, args.verbose)

    # color the bed file 
    input_bed = os.path.join(output_dir, "all_features.bed")
    output_bed = os.path.join(output_dir, "all_features_colored.bed")
    color_features(input_bed, output_bed)

    # extend genes upstream to get 'upstream_of_gene' features and clip by existing annotations
    upstream_bed = annotate_upstream_regions(bed_files['gene'], output_dir, genome_file, 2000)
    clipped_upstream_bed = clip_features(upstream_bed, output_bed, output_dir, "upstream")

    # merge all BED files for introns, exons, dogs, clipped_upstream_bed 
    final_merged_bed = merge_all_features([os.path.join(output_dir, "merged_all_features.bed"), clipped_upstream_bed], output_dir)

    # Handle antisense features
    antisense_bed = duplicate_to_minus_strand(final_merged_bed, output_dir)
    clipped_antisense_bed = clip_features(antisense_bed, final_merged_bed, output_dir, "antisense")
    final_merged_bed = merge_all_features([final_merged_bed, clipped_antisense_bed], output_dir)

    logging.info("Genome indexing and feature extension complete.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create comprehensive genome annotation from a GTF file.")
    parser.add_argument("-g", "--gtf_file", required=True, help="Input GTF file path.")
    parser.add_argument("-f", "--fasta_index_file", required=True, help="Path containing the fasta index for the genome.")
    parser.add_argument("-o", "--output_bed", required=True, help="Path to store output bed files.")
    parser.add_argument("--extend_bases", type=int, default=5000, help="Number of bases to extend gene regions to create DOGs.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose output for debugging purposes.")


    args = parser.parse_args()
    main(args)
