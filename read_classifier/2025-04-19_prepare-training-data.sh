#!/bin/bash 

export bam_in="/g/data/lf10/as7425/2023_mrna-biogenesis-maps/analysis/2024-03-28_ASR012_HeLa-total-RNA004-rep1/ASR012_HeLa-total-RNA004-rep1_primary_genome_alignments_modCalled.bam"
export anno_in="/g/data/lf10/as7425/genomes/human_genome/ensembl_release_110/Homo_sapiens.GRCh38.110.chr.gtf"

# Output directory and files
work_dir="/g/data/lf10/as7425/2023_mrna-biogenesis-maps/analysis/2024-03-28_ASR012_HeLa-total-RNA004-rep1/ASR012_HeLa-total-RNA004/classify_reads_tests/"; mkdir -p ${work_dir} 2>/dev/null
bam_out="${work_dir}/filtered_reads.bam"
gtf_out="${work_dir}/gene_annotation.gtf"

# Gene of interest
target_gene="RPL41"

# Extract GTF entry for the target gene and ensure it is only one unique entry
gtf_line=$(awk '$3 == "gene" && /gene_name "'$target_gene'";/' $anno_in)
echo "$gtf_line" > "$gtf_out"
if [ $(echo "$gtf_line" | wc -l) -ne 1 ]; then
    echo "Error: The gene name query returned multiple entries or none; please check the gene name."
    exit 1
fi

# Parse GTF entry to get coordinates
read chr start end strand <<< $(echo $gtf_line | awk '{print $1, $4, $5, $7}')

# Adjust coordinates for ±1000 nt
start=$((start - 1000))
end=$((end + 1000))

# Filter BAM file for reads overlapping the extended gene region
samtools view -b "$bam_in" $chr:$start-$end -o $bam_out

echo "GTF and BAM files for $target_gene have been successfully created in $work_dir."

cat $anno_in | awk '/gene_name "'$target_gene'";/' > "$gtf_out"