#!/bin/bash 

# test for SRSF3 
export annotation="/home/150/as7425/R1/test/data/SRSF3.gtf"
export out_dir="/home/150/as7425/R1/test/outputs/SRSF3/create_index"; mkdir -p ${out_dir} 2>/dev/null/
export fai="/g/data/lf10/as7425/genomes/human_genome/ensembl_release_110/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai"

python3 ~/R1/utils/create_genome_index.py -g ${annotation} -f ${fai} -o ${out_dir}/SRSF3_index.bed

##################################################

# test for hg38 ensembl 110 

export annotation="/g/data/lf10/as7425/genomes/human_genome/ensembl_release_110/Homo_sapiens.GRCh38.110.chr.gtf"
export out_dir="/g/data/lf10/as7425/2023_mrna-biogenesis-maps/analysis/2024-05-10_R1-genome-index"; mkdir -p ${out_dir} 2>/dev/null/
export fai="/g/data/lf10/as7425/genomes/human_genome/ensembl_release_110/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai"

time python3 ~/R1/utils/create_genome_index.py -g ${annotation} -f ${fai} -o ${out_dir}/hg38_index.bed
