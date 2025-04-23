




export bam="/home/150/as7425/R1/test/data/ACTG1.bam"
export gtf="/home/150/as7425/R1/test/data/ACTG1.gtf"
export out_dir="/home/150/as7425/R1/test/outputs/ACTG1/"; mkdir -p ${out_dir} 2>/dev/null
export output="${out_dir}/output_table.tsv"

export bam="/home/150/as7425/R1/test/data/ACTG1_nascent.bam"
cat <(samtools view -H $bam) <(samtools view $bam | grep "3638fdcf") | samtools view -b > ~/R1/test/data/tmp.bam
samtools index  ~/R1/test/data/tmp.bam

# single read test 
python3 ~/R1/R1.py -b ~/R1/test/data/tmp.bam -g ${gtf} -o ${output}


# run for a single gene 
python3 ~/R1/R1.py -b ${bam} -g ${gtf} -o ${output}

cat modifications.csv | cut -f2-4




cat <(samtools view -H $bam) <(samtools view $bam | grep "0fc2829e-c896") | samtools view -b > ~/R1/test/data/tmp.bam




# 2056f743-ec24-41e4-9d3c-4581cd8595a6    16      17      81513001        60      86S36M1D5M1I1M1I22M1D102M1I50M2I10M1D30M1D7M2D1M1I20M1I7M     \
#  *       0       0       GGGGAAGTTAGGGATGTATGTAGGTAGAGGGGGATGGATGATATGAATAAAATAAATAAGGGATGGGATGGGAGGGAGGGGATTTTCCACCCCCGCGCGCTCGCGACCCTGCGCGGGCCGGCGGCGGAGCTCCGAGTTGGGGCGCCCTCCGAGGCCGCCGGGGAGGCCGAAGGGCTGACGGGCCTGGCCCCTCCCCGGGACTGCTGCGCCGTGGGGAGGGCCCTGCTGCGCCCCGAAACTGCCCGACCCGGGGCGGGGGCCGCGCCGGAGCTGGGGTGGGTCCCCGAGTCCCCGGCCACGCTGCGAAGGGCTTTGCTCCTGGGACGTCCCTTGCAATCTTTCCCCTCGGCTCCAATGGATCCCGGGCGCCAGTGCCGGGGCTGA      \
# $%&--''&&%%*,)%%$$%%%$%%$####$$#%&&&&%&()(,+(&%$&)+**/10.-+,+,(''$&+-00/0//-++,**)'',.<:96621--/../6854)((((@98888<;989:870.-22('))))+*))/.3=>ACEC<?:97665%%))+1.10-7;<=BD934776>>6555AA??;;:::>8647452/?975'%&''==<==E@?<:>.88@8>==>>>F@EB@996666:*(''(<<;;;8I<>G=C>ABFGGA:9:7:<<6=ABBADFFB@86---68004326669;A;91/-...0@>;=@@?+,44353002156?AC?=<?9>>>:58833444>4++****5820///;;;:9/---.:9))0/&   
# NM:i:20 ms:i:236        AS:i:235        nn:i:0  tp:A:P  cm:i:46 s1:i:228        s2:i:0  de:f:0.0596     rl:i:0  qs:i:10 du:f:3.5895     ns:i:14358      ts:i:2150  mx:i:1  ch:i:2586       st:Z:2024-03-20T05:30:06.350+00:00      rn:i:3390       fn:Z:PAU42445_ff6565ed_35211b33_44.pod5 sm:f:73.7878    sd:f:21.393     sv:Z:quantile   dx:i:0  RG:Z:35211b33233c60d557e96d0f60fe3368d79d635d_rna004_130bps_sup@v3.0.1  pt:i:30 MN:i:384  
# MM:Z:A+a?,19,0;  ML:B:C,0,0      YC:Z:255,0,112

####################################################################################################

# test 2024-05-14

export wd="/scratch/lf10/as7425/test_mods"
export in_bam="/g/data/lf10/as7425/2023_mrna-biogenesis-maps/analysis/2024-03-20_ASR011_HeLa-POINT-RNA004-rep1/ASR011_HeLa-POINT-RNA004-rep1_primary_genome_alignments_modCalled.bam"
export gtf="/g/data/lf10/as7425/genomes/human_genome/ensembl_release_110/Homo_sapiens.GRCh38.110.chr_patch_hapl_scaff.gtf"

export target_read="dce8fcd4-e99c-4df4-bd0e-bd25f7183481"
cat <(samtools view -H ${in_bam}) <(samtools view -@ 104 ${in_bam} | grep "${target_read}") | samtools view -b > ${wd}/input.bam; samtools index ${wd}/input.bam; samtools view ${wd}/input.bam

python3 ~/R1/R1.py -b ${wd}/input.bam -g ${gtf} -o ${wd}/out.txt 
