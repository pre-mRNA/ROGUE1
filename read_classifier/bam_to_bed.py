import os 
import logging
from tempfile import NamedTemporaryFile
import subprocess
import pandas as pd 
from concurrent.futures import ThreadPoolExecutor, as_completed, ProcessPoolExecutor

import pysam 
from run_command import run_command

# convert a bam file to bed 
def bam_to_bed(bam_file, output_dir, num_files):
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    temp_bam_files = [NamedTemporaryFile(dir=output_dir, delete=False, suffix=f"_temp_{i}.bam").name for i in range(num_files)]
    temp_bed_files = [NamedTemporaryFile(dir=output_dir, delete=False, suffix=f"_temp_{i}.bed").name for i in range(num_files)]
    
    seq_lengths = {}
    splice_counts = {}
    bam_contents = [[] for _ in range(num_files)]

    # Single pass over the bam file
    for i, read in enumerate(bam):
        seq_name = read.query_name
        seq_length = len(read.query_sequence)
        cigar_string = read.cigarstring
        
        # Populate seq_lengths
        seq_lengths[seq_name] = seq_length

        # count splice events 
        splice_counts[seq_name] = cigar_string.count('N')

        # Distribute alignments
        bam_contents[i % num_files].append(read)

    print([len(bam_content) for bam_content in bam_contents])

    # Create temporary BAM files
    for i in range(num_files):
        with pysam.AlignmentFile(temp_bam_files[i], "wb", header=bam.header) as out_bam:
            for read in bam_contents[i]:
                out_bam.write(read)

    bam.close()

    # Convert BAM files to BED in parallel
    with ThreadPoolExecutor(max_workers=48) as executor:
        for i in range(num_files):
            executor.submit(run_command, f"bedtools bamtobed -i {temp_bam_files[i]} -splitD > {temp_bed_files[i]}", f"Converting BAM to BED: {i + 1}/{num_files}")

    # Concatenate all BED files
    final_bed_file = NamedTemporaryFile(dir=output_dir, delete=False, suffix="_final.bed").name
    cmd = "cat " + " ".join(temp_bed_files) + f" > {final_bed_file}"
    
    logging.info("Running: Merging and sorting all BED files")
    process = subprocess.Popen(cmd, shell=True)
    process.wait()

    logging.info("Calculating alignment lengths")

    # Compute total aligned length for each read
    bed_df = pd.read_csv(final_bed_file, sep="\t", header=None, low_memory=False)
    bed_df.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand']
    bed_df['aligned_length'] = bed_df['end'] - bed_df['start']

    ## fetch the 3' end coordinate for each read 
    # first, split into positive and negative strand dataframes
    positive_strand_bed_df = bed_df[bed_df['strand'] == '+'].sort_values('end', ascending=False)
    negative_strand_bed_df = bed_df[bed_df['strand'] == '-'].sort_values('start')

    # dict to store the result
    read_end_coordinate = {}

    # process the rows of positive_strand_bed_df
    for row in positive_strand_bed_df.itertuples(index=False):
        # if this name is not yet in the dictionary, add it with the 'chrom:end:strand' value
        if getattr(row, 'name') not in read_end_coordinate:
            read_name = getattr(row, 'name')
            chrom = getattr(row, 'chrom')
            end = getattr(row, 'end')
            strand = getattr(row, 'strand')
            read_end_coordinate[read_name] = f'{chrom}:{end}:{strand}'

    # process the rows of negative_strand_bed_df
    for row in negative_strand_bed_df.itertuples(index=False):
        # if this name is not yet in the dictionary, add it with the 'chrom:start:strand' value
        if getattr(row, 'name') not in read_end_coordinate:
            read_name = getattr(row, 'name')
            chrom = getattr(row, 'chrom')
            start = getattr(row, 'start')
            strand = getattr(row, 'strand')
            read_end_coordinate[read_name] = f'{chrom}:{start}:{strand}'

    # calculate the total aligned length for each read 
    with ThreadPoolExecutor(max_workers=48) as executor:
        total_aligned_length = bed_df.groupby('name')['aligned_length'].sum().to_dict()

    logging.info("Calculating sequence lengths")

    # Assign sequence length, total aligned length, and splice counts to column 5
    bed_df['s-a_length'] = bed_df['name'].map(seq_lengths).astype(str) + ',' + bed_df['name'].map(total_aligned_length).astype(str) + ',' + bed_df['name'].map(splice_counts).astype(str)
    bed_df = bed_df[['chrom', 'start', 'end', 'name', 's-a_length', 'strand']]
    bed_df.to_csv(final_bed_file, sep="\t", index=False, header=False)

    # Remove temporary files
    for temp_bam in temp_bam_files:
        os.remove(temp_bam)
    for temp_bed in temp_bed_files:
        os.remove(temp_bed)

    return final_bed_file, read_end_coordinate