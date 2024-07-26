#!/bin/bash 

# export target_gene="SRSF3"
# export bam="/home/150/as7425/R1/test/data/${target_gene}.bam"
# export gtf="/home/150/as7425/R1/test/data/${target_gene}.gtf"
export out_dir="${scratch}/R1/test/outputs/${target_gene}/"; mkdir -p ${out_dir} 2>/dev/null
export output="${out_dir}/output_table.tsv"
export index="/g/data/lf10/as7425/2023_mrna-biogenesis-maps/analysis/204-06-17_ROGUE1-splicing-order/annotation_with_clusters"

# bam="/home/150/as7425/R1/test/data/merged.bam"
# gtf="/home/150/as7425/R1/test/data/merged.gtf"

# run R1 without calculating modifications
# time python3 ~/R1/R1.py -b ${bam} -g ${gtf} -o ${output} -p -j > ${out_dir}/ROGUE1-log.txt 2>&1

# myc
export bam="/home/150/as7425/R1/test/data/myc_point.bam"
export gtf="/home/150/as7425/R1/test/data/myc.gtf"

# run R1 with index 
time python3 ~/R1/R1.py -b ${bam} -g ${gtf} -o ${output} -p -j --index ${index} --record_exons #> ${out_dir}/ROGUE1-log.txt 2>&1
wc -l $output

# classify the splicing from acceptor 
time python3 ~/R1/classify_splicing/run_classify_splicing_from_acceptor.py ${output} "${output%.*}_classify.tsv"
wc -l "${output%.*}_classify.tsv"

# add tags for 3' features 
time python3 ~/R1/color_bam/color_bam_from_3prime_feature.py -ibam $bam -obam ${bam%.*}_color_threeprime.bam -table ${output} # use original R1 output 

# color bam from splicing status 
time python3 ~/R1/color_bam/color_bam_from_list.py -ibam "${bam%.*}_color_threeprime.bam" -obam "${bam%.*}_tag_threeprime_color_splicing.bam" -class_file "${output%.*}_classify.tsv"

# add RE and AL tags 
time python3 ~/R1/igv/create_bam_re_tag.py -ibam "${bam%.*}_tag_threeprime_color_splicing.bam" -obam "${bam%.*}_tag_threeprime_color_splicing_RE-AL.bam"

# color by distance to donor 
time python3 ~/R1/color_bam/color_bam_from_donor_distance.py -q 50,100 -ibam "${bam%.*}_tag_threeprime_color_splicing_RE-AL.bam" -obam "${bam%.*}_tag_threeprime_color_splicing_RE-AL_pA-dist.bam" --class-file ${output}

##########################################################################################
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