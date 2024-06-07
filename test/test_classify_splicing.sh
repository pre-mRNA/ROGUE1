#!/bin/bash 

export target_gene="SRSF3"
export bam="/home/150/as7425/R1/test/data/${target_gene}.bam"
export gtf="/home/150/as7425/R1/test/data/${target_gene}.gtf"
export out_dir="/home/150/as7425/R1/test/outputs/${target_gene}/"; mkdir -p ${out_dir} 2>/dev/null
export output="${out_dir}/output_table.tsv"

python3 ~/R1/R1.py -b ${bam} -g ${gtf} -o ${output}

python3 ~/R1/classify_splicing/classify_splicing.py ${output} "${out_dir}/classify_splicing_srsf3.txt" ${gtf} 

# Using thresholds - Splice: 1, Intron Positive: 75, Intron Negative: 25
# Initial read count: 945
# Filtered read count (protein-coding with introns): 945
# Spliced reads: 655
# Intron retained reads: 18
# Ambiguous reads: 272
# Discarded reads: 0

export anno_in="/g/data/lf10/as7425/genomes/human_genome/ensembl_release_110/Homo_sapiens.GRCh38.110.chr.gtf"


python3 ~/R1/classify_splicing/classify_splicing.py ${output} "${out_dir}/classify_splicing_srsf3.txt" ${anno_in} 

python3 ~/R1/color_bam/color_bam_from_list.py -ibam $bam -obam ${out_dir}/SRSF3_splicing)_color.bam -class_file 