import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

from classify_splicing_from_acceptor import classify_splicing

def setup_logging(log_file, log_level):
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )

def validate_file(file_path, expected_extension):
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    if path.suffix.lower() != expected_extension:
        raise ValueError(f"The file {file_path} should have a {expected_extension} extension.")
    return path

def process_file(input_file, output_file):
    logging.info(f"Reading input file: {input_file}")
    try:
        df = pd.read_csv(input_file, sep='\t')
    except Exception as e:
        logging.error(f"Error reading input file: {e}")
        raise

    required_columns = ['splice_count', 'exon_ids', 'intron_ids']
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Input file is missing required columns: {', '.join(missing_columns)}")

    logging.info("Applying splicing classification from acceptor")
    df = classify_splicing(df)

    logging.info(f"Writing output to: {output_file}")
    df.to_csv(output_file, sep='\t', index=False)
    logging.info("Processing complete")

def main():
    parser = argparse.ArgumentParser(description="Classify splicing in a TSV file")
    parser.add_argument("input_file", help="Path to the input TSV file")
    parser.add_argument("output_file", help="Path to the output TSV file")
    parser.add_argument("--log-file", help="Path to the log file (default: output_file_name.log)")
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="Set the logging level")
    args = parser.parse_args()

    try:
        input_file = validate_file(args.input_file, '.tsv')
        output_file = Path(args.output_file)
        
        # default log file 
        if args.log_file is None:
            log_file = output_file.with_suffix('.log')
        else:
            log_file = Path(args.log_file)

        setup_logging(log_file, args.log_level)
        
        logging.info(f"Log file: {log_file}")
        process_file(input_file, output_file)
        
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()