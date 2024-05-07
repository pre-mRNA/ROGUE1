from tempfile import NamedTemporaryFile
from run_command import run_command

# generate a genome file for bedtools intersect 
def generate_genome_file(bam_file, output_dir):
    genome_file = NamedTemporaryFile(dir=output_dir, delete=False, suffix="_genome.txt").name
    cmd = f"samtools view -H {bam_file} | grep @SQ | sed 's/@SQ\tSN://' | sed 's/LN://' | sort -k1,1 --parallel=104 --buffer-size=80G > {genome_file}"
    run_command(cmd, "Generating genome file")
    return genome_file

# Read list of chromosomes from genome_file
def read_chromosomes_from_genome_file(genome_file):
    with open(genome_file, 'r') as f:
        chromosomes = [line.split()[0] for line in f]
    return chromosomes

# Filter bed file based on chromosome list
def filter_bed_by_chromosome_inplace(bed_file, chromosomes):
    with open(bed_file, 'r') as f:
        lines = f.readlines()

    filtered_lines = [line for line in lines if line.split()[0] in chromosomes]

    with open(bed_file, 'w') as f:
        f.writelines(filtered_lines)