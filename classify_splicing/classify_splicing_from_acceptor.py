import pandas as pd
import numpy as np
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def classify_splicing(df: pd.DataFrame) -> pd.DataFrame:
    """
    Classify RNA splicing events based on exon and intron IDs within sequencing reads.

    Parameters:
        df (pd.DataFrame): DataFrame containing the following columns:
            - 'exon_ids' (str): Comma-separated exon IDs.
            - 'intron_ids' (str): Comma-separated intron IDs.
            - 'splice_count' (int): Number of splicing events.
            - 'read_id' (str): Unique identifier for each read.

    Returns:
        pd.DataFrame: Original DataFrame with an additional column 'splicing_classification'.
    """
    required_columns = {'exon_ids', 'intron_ids', 'splice_count', 'read_id'}
    missing = required_columns - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    def parse_ids(id_string: str) -> set:
        """
        Parse a comma-separated string of IDs into a set of integers.

        Parameters:
            id_string (str): Comma-separated string of IDs.

        Returns:
            set: Set of integer IDs.
        """
        try:
            return set(int(x.strip()) for x in id_string.split(',') if x.strip().isdigit()) if pd.notna(id_string) else set()
        except ValueError as ve:
            logger.warning(f"ValueError parsing IDs '{id_string}': {ve}")
            return set()
        except Exception as e:
            logger.error(f"Unexpected error parsing IDs '{id_string}': {e}")
            return set()

    def classify_row(row):
        """
        Classify a single row based on exon_ids, intron_ids, and splice_count.

        Parameters:
            row (pd.Series): A row from the DataFrame.

        Returns:
            str or np.nan: Splicing classification or NaN if invalid.
        """
        try:
            exons = parse_ids(row['exon_ids'])
            introns = parse_ids(row['intron_ids'])
            splice_count = row['splice_count']
            read_id = row['read_id']

            # Ensure splice_count is an integer
            if not isinstance(splice_count, (int, np.integer)):
                logger.warning(f"Read ID: {read_id}, Invalid splice_count: {splice_count}")
                return np.nan

            if splice_count == 0:
                if len(exons) == 0 and len(introns) == 0: # no features 
                    return 'ambiguous'
                elif len(exons) == 1 and len(introns) == 0: # single exon only 
                    return 'ambiguous'
                elif len(introns) == 1 and len(exons) <= 1: # single intron only 
                    intron = next(iter(introns))
                    exon = next(iter(exons)) if exons else None
                    if exon is not None and intron + 1 == exon:
                        return 'fully-unspliced' # intron-exon boundary has been transcribed 
                    else:
                        return 'ambiguous' # single intron only 
                elif any(i + 1 in exons for i in introns): # no splicing in read; and intron-exon boundary has been transcribed 
                    return 'fully-unspliced' 
                elif len(exons) > 1 and len(introns) == 0: # error case where we see two introns but splice_count = 0 
                    return 'ambiguous'
                else:
                    logger.warning(f"Read ID: {read_id}, Unexpected case: splice_count=0, exons={exons}, introns={introns}")
                    return 'ambiguous'
            
            elif splice_count > 0:
                if len(exons) >= 2:
                    return 'spliced'
                elif len(exons) == 1:
                    return 'ambiguous'
                elif len(exons) == 0 and len(introns) == 0:
                    return 'ambiguous'
                elif len(exons) == 0 and len(introns) > 0:
                    # **Modification:** Directly classify as 'ambiguous' without raising a warning
                    return 'ambiguous'
                else:
                    logger.warning(f"Read ID: {read_id}, Unexpected case: splice_count>0, exons={exons}, introns={introns}")
                    return 'ambiguous'
            
            else:
                logger.warning(f"Read ID: {read_id}, Unexpected splice_count: {splice_count}")
                return 'ambiguous'

        except Exception as e:
            logger.error(f"Read ID: {row.get('read_id', 'Unknown')}, Exception during classification: {e}")
            return np.nan

    # apply the classification function to each row
    df = df.copy()  # To avoid modifying the original DataFrame
    df['splicing_classification'] = df.apply(classify_row, axis=1)
    return df
