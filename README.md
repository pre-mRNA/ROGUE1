# ROGUE1

ROGUE1 performs single-molecule integration, classification and quantification of splicing, transcription, and RNA modification data from long-read sequencing data. 

# Quick-start guide 

```

# map long reads to genome and retain primary mappings 
minimap2 -a -x splice "${reference.fa}" "${reads.fq}" > ./alignments.sam
samtools view -b -F 260 ./alignments.sam | samtools sort > ./primary.bam
samtools index primary.bam 

# run ROGUE1 in 1st pass-mode 
python3 ./R1/R1.py integrate -b $primary.bam -g ${annotation.gtf} -o ./R1_table.tsv

# refine transcript models from multiple R1 results
python3 ./R1/R1.py refine-annotation 



```

# File Structure
TODO: update file structure 
```
r1.py                         # ROGUE1 master script
read_classifier/
├── bam_to_bed.py             # Module to convert BAM to BED and process alignments
├── gtf_to_bed.py             # Module for converting GTF files to BED files
├── parse_output.py           # Script to parse output files and summarize data
├── process_genome.py         # Utilities to generate and process genome-related files
├── run_bedtools.py           # Main script to run BEDTools operations
└── utilities/
    ├── logging_config.py     # Configuration for logging across scripts
    └── run_command.py        # Utility to run shell commands

color_bam/
    └── color_bam_from_lis.py # Set alignment YC color bam tag from list of read IDs 

utils/
    └── run_mod_extraction.py # Extract the genome positions and probabilities of single-read MM and ML tags 
```

# Running Modification Extraction

R1 convers modification positions embedded in MM/ML tags to their single-read genome positions. 

```
python3 utils/run_mod_extraction.py --bam_file path/to/bamfile.bam --output_file path/to/output_modifications.txt
```

__Output format__: A tab-delimited table of: 
``` 
read_id     read_chromosoome        read_strand     mod_position:probability,mod_position:probability
```
Column 4 contains tuples of modification_site:probability, separated by comma delimiters 

Currently, this only works for m6A predictions generated by Dorado v0.6
TODO: test and implement for other mod types 

# Dependencies

### TODO: update dependencies 
```
Python 3.8 or higher
Pandas
Pysam
NumPy
BEDTools
Samtools
```

# Input Data:

### TODO: update input data requirements
BAM File: Prepare a BAM file containing RNA sequencing reads.
GTF File: A GTF file containing gene annotations.

# Usage: 

### TODO: update usage
Run the main script run_bedtools.py with the necessary arguments:

```
python r1.py --bam_file path/to/bamfile.bam --gtf_file path/to/annotations.gtf --output_dir path/to/output_table.txt
```

# Output:

### TODO: update outputs 


## Creating a genome undex 
``` 
python3 ./utils/create_genome_index.py
```
__output__: A bed file that annotes each part of the genome as intron, exon, or downstream of gene 

## filtering GTF annotaiton: 
```
python3 ./filter_annotation/
```
TODO: Test and document this code 

## Setting Color Tags in BAM Files

To set color tags for specific reads in a BAM file, use the `bam_set_color.py` script. This utility requires a list of read IDs and the RGB color code that will be applied as a tag to these reads.


