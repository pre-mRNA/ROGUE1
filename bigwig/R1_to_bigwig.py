import pandas as pd
import pyBigWig
import os
import argparse

def create_bigwig(df, strand, output_file):

    strand_df = df[df['read_end_strand'] == strand]
    
    grouped = strand_df.groupby(['read_end_chromosome', 'read_end_position']).size().reset_index(name='count')

    grouped = grouped.sort_values(['read_end_chromosome', 'read_end_position'])
    
    chromosomes = grouped['read_end_chromosome'].unique()
    
    chrom_sizes = {chrom: grouped[grouped['read_end_chromosome'] == chrom]['read_end_position'].max() for chrom in chromosomes}
    
    # create the bigWig file
    bw = pyBigWig.open(output_file, 'w')
    bw.addHeader(list(chrom_sizes.items()))
    
    for chrom in chromosomes:
        chrom_data = grouped[grouped['read_end_chromosome'] == chrom]
        starts = chrom_data['read_end_position'].values - 1  # 0-based start
        values = chrom_data['count'].values.astype(float)  
        
        if len(starts) > 0:
            try:
                # span size can be changed for visibility
                bw.addEntries(chrom, starts.tolist(), values=values.tolist(), span=1)
            except RuntimeError as e:
                print(f"Error adding entries for chromosome {chrom}: {e}")
                print(f"Data shape: starts {starts.shape}, values {values.shape}")
                print(f"First few entries: starts {starts[:5]}, values {values[:5]}")
        else:
            print(f"Skipping chromosome {chrom} due to insufficient data.")
    
    bw.close()

def main():
    parser = argparse.ArgumentParser(description="Create bigWig files of read end position from R1 bulkfile.")
    parser.add_argument("-i", "--input", required=True, help="input R1 bulkfile (TSV)")
    parser.add_argument("-o", "--output_dir", required=True, help="output directory")
    args = parser.parse_args()

    try:
        df = pd.read_csv(args.input, sep='\t')
        print(f"Successfully read input file. Shape: {df.shape}")
    except Exception as e:
        print(f"Error reading input file: {e}")
        return

    if 'read_end_chromosome' not in df.columns or 'read_end_position' not in df.columns or 'read_end_strand' not in df.columns:
        print("Error: Input file must contain 'read_end_chromosome', 'read_end_position', and 'read_end_strand' columns")
        return

    os.makedirs(args.output_dir, exist_ok=True)

    prefix = os.path.splitext(os.path.basename(args.input))[0]

    plus_file = os.path.join(args.output_dir, f'{prefix}_plus_strand.bw')
    minus_file = os.path.join(args.output_dir, f'{prefix}_minus_strand.bw')

    print("Creating bigWig file for plus strand...")
    create_bigwig(df, '+', plus_file)
    print("Creating bigWig file for minus strand...")
    create_bigwig(df, '-', minus_file)

    print(f"bigWig files created successfully in '{args.output_dir}':")
    print(f"  {os.path.basename(plus_file)}")
    print(f"  {os.path.basename(minus_file)}")

if __name__ == "__main__":
    main()