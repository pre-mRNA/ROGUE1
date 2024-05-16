import pandas as pd
import argparse

# classify splicing based on given thresholds
def classify_splicing(df, splice_threshold, intron_positive_threshold, intron_negative_threshold):
    
    # define conditions for classification
    conditions = [
        (df['splice_count'] >= splice_threshold) & (df['intronic_alignment'] <= intron_negative_threshold),
        (df['intronic_alignment'] > intron_positive_threshold)
    ]
    
    # define corresponding classifications
    choices = ['spliced', 'intron_retained']
    
    # default to 'ambiguous' if none of the above conditions are met
    df['splicing_classification'] = pd.Series(['ambiguous'] * len(df))
    
    df.loc[conditions[0], 'splicing_classification'] = choices[0]
    df.loc[conditions[1], 'splicing_classification'] = choices[1]
    
    return df

# load the dataset
def load_data(file_path):
    # read the file using pandas
    return pd.read_csv(file_path, sep='\t')

def save_data(df, file_path):
    df.to_csv(file_path, sep='\t', index=False)


def main(args):

    df = load_data(args.input_table)
    
    # classify splicing
    df = classify_splicing(df, args.splice_threshold, args.intron_positive_threshold, args.intron_negative_threshold)
    
    # save updated data
    save_data(df, args.output_table)
    
    print(f"Data processed and saved to {args.output_table}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Classify splicing status based on alignment metrics")
    parser.add_argument("-i", "--input_table", type=str, required=True, help="input table containing splicing data")
    parser.add_argument("-o", "--output_table", type=str, required=True, help="output table to save classified data")
    parser.add_argument("-st", "--splice_threshold", type=int, default=1, help="threshold for splice count to classify as spliced")
    parser.add_argument("-ipt", "--intron_positive_threshold", type=int, default=40, help="threshold for intronic alignment to classify as intron retained")
    parser.add_argument("-int", "--intron_negative_threshold", type=int, default=10, help="upper limit of intronic alignment to classify as spliced")
    
    args = parser.parse_args()
    
    main(args)
