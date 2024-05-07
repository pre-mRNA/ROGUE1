#!/bin/bash 

export target_gene="SRSF3"
export bam="/home/150/as7425/R1/test/data/${target_gene}.bam"
export gtf="/home/150/as7425/R1/test/data/${target_gene}.gtf"
export out_dir="/home/150/as7425/R1/test/outputs/${target_gene}/"; mkdir -p ${out_dir} 2>/dev/null
export output="${out_dir}/output_table.tsv"

python3 ~/R1/R1.py -b ${bam} -g ${gtf} -o ${output}