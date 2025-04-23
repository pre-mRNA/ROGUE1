
export bam="/home/150/as7425/R1/test/data/myc_point_modcalled_0.995.bam"
export script="/home/150/as7425/R1/parse_modifications_data/filter_empty_MM_tags.py"

# filter for one read 
# samtools view -b <(cat <(samtools view -H $bam) <(samtools view $bam | grep "1bb32a64-48f1-48b6-8a0d-702af474ff4c")) > /home/150/as7425/R1/test/data/myc_point_modcalled_0.995_one_read.bam
# export bam="/home/150/as7425/R1/test/data/myc_point_modcalled_0.995_one_read.bam"

time python3 $script -ibam "$bam" -obam "${bam%.*}_filt_MM.bam" -t 48 

# samtools view $bam
# samtools view "${bam%.*}_filt_MM.bam" 