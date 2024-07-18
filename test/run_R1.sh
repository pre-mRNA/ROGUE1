#!/bin/bash 

export target_gene="SRSF3"
export bam="/home/150/as7425/R1/test/data/${target_gene}.bam"
export gtf="/home/150/as7425/R1/test/data/${target_gene}.gtf"
export out_dir="/home/150/as7425/R1/test/outputs/${target_gene}/"; mkdir -p ${out_dir} 2>/dev/null
export output="${out_dir}/output_table.tsv"


# bam="/home/150/as7425/R1/test/data/merged.bam"
# gtf="/home/150/as7425/R1/test/data/merged.gtf"

# run R1 without calculating modifications
time python3 ~/R1/R1.py -b ${bam} -g ${gtf} -o ${output} -p -j > ${out_dir}/ROGUE1-log.txt 2>&1


# run R1 with index 
time python3 ~/R1/R1.py -b ${bam} -g ${gtf} -o ${output} -p -j --index /g/data/lf10/as7425/2023_mrna-biogenesis-maps/analysis/204-06-17_ROGUE1-splicing-order/annotation_with_clusters > ${out_dir}/ROGUE1-log.txt 2>&1





# run R1 while calculating modifications 
time python3 ~/R1/R1.py -b ${bam} -g ${gtf} -o ${output} -m -j 
wc -l ${output}
# 946 lines 
# real    0m7.534s 

# run R1 while calculate polyA tail 
time python3 ~/R1/R1.py -b ${bam} -g ${gtf} -o ${output} -p
wc -l ${output}
# 946 lines 
# real    0m7.534s 