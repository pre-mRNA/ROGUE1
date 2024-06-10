#!/bin/bash 

export target_gene="SRSF3"
export bam="/home/150/as7425/R1/test/data/${target_gene}.bam"
export gtf="/home/150/as7425/R1/test/data/${target_gene}.gtf"
export out_dir="/home/150/as7425/R1/test/outputs/${target_gene}/"; mkdir -p ${out_dir} 2>/dev/null
export output="${out_dir}/output_table.tsv"

time python3 ~/R1/R1.py -b ${bam} -g ${gtf} -o ${output}
# 0m19.332s

# classify splicing status from ROGUE1 table
time python3 ~/R1/classify_splicing/classify_splicing.py ${output} "${out_dir}/classify_splicing_srsf3.txt" ${gtf} 
# real    0m0.925s

# Using thresholds - Splice: 1, Intron Positive: 75, Intron Negative: 25
# Initial read count: 945
# Filtered read count (protein-coding with introns): 945
# Spliced reads: 655
# Intron retained reads: 18
# Ambiguous reads: 272
# Discarded reads: 0

# test the process using the entire genome annotation 
export anno_in="/g/data/lf10/as7425/genomes/human_genome/ensembl_release_110/Homo_sapiens.GRCh38.110.chr.gtf"
time python3 ~/R1/classify_splicing/classify_splicing.py ${output} "${out_dir}/classify_splicing_srsf3.txt" ${anno_in} 

# color the bam file using the splicing classification table 
# we can set the read group tag instead of overwriting the YC tag?
python3 ~/R1/color_bam/color_bam_from_list.py -ibam $bam -obam ${out_dir}/SRSF3_splicing_color.bam -class_file "${out_dir}/classify_splicing_srsf3.txt"

# split the bam from splicing status
time python3 ~/R1/classify_splicing/split_bam_from_splicing_class.py -i "${out_dir}/SRSF3_splicing_color.bam" -o "${out_dir}" -s "SRSF3_split" 