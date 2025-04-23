import pandas as pd
import numpy as np
import warnings
from tqdm import tqdm

tqdm.pandas()

def parse_ids(id_string):
    """
    Convert a comma-separated string of integer IDs into a sorted list of ints.
    Returns an empty list if it's missing or invalid.
    """
    if pd.isna(id_string) or id_string.strip() == '':
        return []
    return sorted(int(x) for x in id_string.split(',') if x.strip().isdigit())


def classify_previous_splicing(
    df, 
    gene_to_introns,     
    verbose=False
):
    """
    Classify each read as 'cryptic' or 'canonical' based on exon_ids, intron_ids, splice_count, 
    and three_prime_feature. Also sub-classify 'unspliced' as 'completely_unspliced' if the read 
    spans ALL introns in that gene, or 'fully_unspliced' otherwise.

    Parameters:
    -----------
    df (pd.DataFrame): Must contain:
        - 'gene_id' (str)
        - 'exon_ids' (str): Comma-separated exon IDs
        - 'intron_ids' (str): Comma-separated intron IDs
        - 'splice_count' (int): Number of splice events
        - 'three_prime_feature' (str): 'exon_*', 'intron_*', or 'region_DOG'
        - 'read_id' (str): Unique identifier for the read

    gene_to_introns (dict): 
        A dict mapping gene_id -> set of all intron numbers for that gene.

    verbose (bool): If True, include the specific cryptic reason in the final output.

    Returns:
    --------
    pd.DataFrame: 
      The original df plus added columns:
        - 'cryptic_splicing' ('cryptic' or 'canonical')
        - 'cryptic_reason'   (if verbose=True)
        - 'previous_splicing_status' 
          which can be 'fully_spliced', 'partially_spliced', 'completely_unspliced',
          'fully_unspliced', or 'ambiguous'
    """

    def classify_row(row):
        exon_ids = parse_ids(row['exon_ids'])
        intron_ids = parse_ids(row['intron_ids'])
        splice_count = row['splice_count']
        three_prime_feature = row['three_prime_feature']
        read_id = row.get('read_id', 'unknown_read')
        
        # get the set of introns for the gene
        gene_id = row.get('gene_id', None)
        # default to empty set if missing
        all_introns_for_gene = gene_to_introns.get(gene_id, set())

        expected_splice_count = 0
        possible_splice_count = 0
        
        # gather features in ascending numeric order
        all_features = []
        for exon_id in exon_ids:
            all_features.append(('exon', exon_id))
        for intron_id in intron_ids:
            all_features.append(('intron', intron_id))
        all_features = sorted(all_features, key=lambda x: (x[1], x[0] == 'intron'))

        if three_prime_feature == 'region_DOG':
            all_features.append(('region_DOG', float('inf')))
        
        # parse the 3' feature
        if three_prime_feature.startswith('exon_'):
            current_type = 'exon'
            current_id = int(three_prime_feature.split('_')[1])
        elif three_prime_feature.startswith('intron_'):
            current_type = 'intron'
            current_id = int(three_prime_feature.split('_')[1])
        elif three_prime_feature == 'region_DOG':
            current_type = 'region_DOG'
            current_id = float('inf')
        else:
            reason = 'Unknown three_prime_feature type.'
            return ('cryptic', reason if verbose else None, None)
        
        # find this feature in the sorted list
        current_index = None
        for idx, (ftype, fid) in enumerate(all_features):
            if ftype == current_type and fid == current_id:
                current_index = idx
                break
        
        if current_index is None:
            reason = "3' end feature not found among exons/introns."
            return ('cryptic', reason if verbose else None, None)
        
        # confirm highest feature
        highest_ftype, highest_fid = all_features[-1]
        if current_type != highest_ftype or current_id != highest_fid:
            reason = '3\' end feature is not the highest-ranked feature.'
            return ('cryptic', reason if verbose else None, None)
        
        # backtrack to check for cryptic splicing
        while current_index > 0:
            ctype, cid = all_features[current_index]
            ptype, pid = all_features[current_index - 1]
            
            if ctype == 'exon':
                if ptype == 'intron':
                    # intron must be that exon - 1
                    if pid != cid - 1:
                        reason = f"Exon {cid} not preceded by intron {cid-1}."
                        return ('cryptic', reason if verbose else None, None)
                    current_index -= 1
                    possible_splice_count += 1
                elif ptype == 'exon':
                    expected_splice_count += 1
                    possible_splice_count += 1
                    current_index -= 1
                else:
                    reason = f"Unexpected preceding feature {ptype} for exon {cid}."
                    return ('cryptic', reason if verbose else None, None)
            
            elif ctype == 'intron':
                if ptype == 'exon':
                    # intron must be preceded by that same ID exon
                    if pid != cid:
                        reason = f"Intron {cid} not preceded by exon {cid}."
                        return ('cryptic', reason if verbose else None, None)
                    current_index -= 1
                else:
                    reason = f"Intron {cid} not preceded by an exon."
                    return ('cryptic', reason if verbose else None, None)
            
            elif ctype == 'region_DOG':
                # must be preceded by highest exon
                if ptype != 'exon' or pid != max(exon_ids) if exon_ids else -9999:
                    reason = "region_DOG not preceded by highest exon."
                    return ('cryptic', reason if verbose else None, None)
                current_index -= 1
            
            else:
                reason = f"Unknown feature type {ctype} in backtrack."
                return ('cryptic', reason if verbose else None, None)
        
        # final check: expected vs actual splices
        if expected_splice_count != splice_count:
            reason = f"Splice count mismatch: expected {expected_splice_count}, got {splice_count}."
            return ('cryptic', reason if verbose else None, None)
        
        # decide final classification
        if possible_splice_count > 0 and splice_count == 0:
            # the read is unspliced
            read_intron_set = set(intron_ids)
            # if the gene has >=1 intron and read_intron_set matches all_introns_for_gene => completely_unspliced

            # print(f"all_introns_for_gene: {all_introns_for_gene}")
            # print(f"read_intron_set: {read_intron_set}")
            # print(f"gene_id: {gene_id}")
            # print(f"read_id: {read_id}")
            # print(f"exon_ids: {exon_ids}")
            # print(f"intron_ids: {intron_ids}")
            # print(f"splice_count: {splice_count}")
            # print(f"possible_splice_count: {possible_splice_count}")
            # print(f"three_prime_feature: {three_prime_feature}")
            # print(f"current_index: {current_index}")
            # print(f"expected_splice_count: {expected_splice_count}")
            # print(f"possible_splice_count: {possible_splice_count}")
            # print(f"current_type: {current_type}")
            # print(f"current_id: {current_id}")
            # print(f"highest_ftype: {highest_ftype}")
            # print(f"highest_fid: {highest_fid}")
            # print(f"current_index: {current_index}")
            # print(f"all_features: {all_features}")
            # print(f"current_index: {current_index}")


            if len(all_introns_for_gene) >= 1 and read_intron_set == all_introns_for_gene:
                return ('canonical', None, 'completely_unspliced')
            else:
                return ('canonical', None, 'fully_unspliced')

        elif possible_splice_count > 0 and splice_count < possible_splice_count:
            return ('canonical', None, 'partially_spliced')
        
        elif possible_splice_count > 0 and splice_count == possible_splice_count:
            return ('canonical', None, 'fully_spliced')
        
        else:
            # possible_splice_count == 0 => ambiguous
            return ('canonical', None, 'ambiguous')

    # use tqdm to classify
    classification = df.progress_apply(classify_row, axis=1, result_type='expand')
    classification.columns = ['cryptic_splicing', 'cryptic_reason', 'previous_splicing_status']
    
    # attach these columns
    df = pd.concat([df, classification], axis=1)
    
    if not verbose:
        df.drop(columns=['cryptic_reason'], inplace=True)
    
    return df

# test cases for the classify_previous_splicing function
def test_classify_splicing():
    """
    Test the classify_previous_splicing function with a variety of test cases to ensure correctness.
    """
    test_cases = [
        # 1. Canonical: Proper splicing with exons and introns
        {'gene_id': 'GENE_01', 'splice_count': 1, 'exon_ids': '2,3', 'intron_ids': '1', 'three_prime_feature': 'exon_3', 'read_id': 'read_1', 'expected': 'canonical'},

        # 2. Cryptic: splice_count exceeds possible based on exons
        {'gene_id': 'GENE_01', 'splice_count': 3, 'exon_ids': '1,2', 'intron_ids': '', 'three_prime_feature': 'exon_2', 'read_id': 'read_2', 'expected': 'cryptic'},

        # 3. Canonical: A read that has exactly 1 intron in gene GENE_02, is unspliced, and the gene has 1 intron total => 'completely_unspliced'
        {'gene_id': 'GENE_02', 'splice_count': 0, 'exon_ids': '', 'intron_ids': '1', 'three_prime_feature': 'intron_1', 'read_id': 'read_3', 'expected': 'canonical'},  

        # 4. Cryptic: Spans 2 introns but splice_count=0 and missing exons between introns
        {'splice_count': 0, 'exon_ids': '', 'intron_ids': '2,3', 'three_prime_feature': 'intron_3', 'read_id': 'read_4', 'expected': 'cryptic'},
        
        # 5. Cryptic: Spans multiple exons and introns with wrong splice_count
        {'splice_count': 2, 'exon_ids': '4,5,6', 'intron_ids': '4,5', 'three_prime_feature': 'exon_6', 'read_id': 'read_5', 'expected': 'cryptic'},
        
        # 6. Cryptic: Spans 1 intron but splice_count=2 exceeds exons spanned
        {'splice_count': 2, 'exon_ids': '', 'intron_ids': '4', 'three_prime_feature': 'intron_4', 'read_id': 'read_6', 'expected': 'cryptic'},
        
        # 7. Cryptic: Spans 3 exons and with all introns retained but splice_count=1 
        {'splice_count': 1, 'exon_ids': '7,8,9', 'intron_ids': '7,8', 'three_prime_feature': 'exon_9', 'read_id': 'read_7', 'expected': 'cryptic'},
        
        # 8. Canonical: Single exon with no splicing
        {'splice_count': 0, 'exon_ids': '10', 'intron_ids': '', 'three_prime_feature': 'exon_10', 'read_id': 'read_8', 'expected': 'canonical'},
        
        # 9. Cryptic: Read ends in DOG but the highest ranked internal feature is an intron 
        {'splice_count': 0, 'exon_ids': '', 'intron_ids': '1', 'three_prime_feature': 'region_DOG', 'read_id': 'read_9', 'expected': 'cryptic'},
        
        # 10. Cryptic: Read in region_DOG with splice_count >0 and no exons/introns
        {'splice_count': 2, 'exon_ids': '', 'intron_ids': '', 'three_prime_feature': 'region_DOG', 'read_id': 'read_10', 'expected': 'cryptic'},
        
        # 11. Canonical: read in intron with splice_count=0 (unspliced)
        {'splice_count': 0, 'exon_ids': '', 'intron_ids': '5', 'three_prime_feature': 'intron_5', 'read_id': 'read_11', 'expected': 'canonical'},
        
        # 12. Cryptic: splice_count exceeds possible with exons and introns
        {'splice_count': 5, 'exon_ids': '11,12', 'intron_ids': '11,12,13', 'three_prime_feature': 'exon_12', 'read_id': 'read_12', 'expected': 'cryptic'},
        
        # 13. Cryptic: Single intron spanned with inappropriate splice_count
        {'splice_count': 1, 'exon_ids': '', 'intron_ids': '6', 'three_prime_feature': 'intron_6', 'read_id': 'read_13', 'expected': 'cryptic'},
        
        # 14. Cryptic: Single exon but splice_count=1 indicates a splice where none can occur
        {'splice_count': 1, 'exon_ids': '14', 'intron_ids': '', 'three_prime_feature': 'exon_14', 'read_id': 'read_14', 'expected': 'cryptic'},
        
        # 15. Canonical: Multiple exons with correct splice_count
        {'splice_count': 3, 'exon_ids': '15,16,17,18', 'intron_ids': '15,16,17', 'three_prime_feature': 'exon_18', 'read_id': 'read_15', 'expected': 'cryptic'},
        
        # 16. Cryptic: Spans multiple introns with exons in between and incorrect splice_count
        {'splice_count': 3, 'exon_ids': '19,20,21,22', 'intron_ids': '19,20,21', 'three_prime_feature': 'exon_22', 'read_id': 'read_16', 'expected': 'cryptic'},
        
        # 17. Cryptic: Spans multiple introns but missing exons between them
        {'splice_count': 2, 'exon_ids': '', 'intron_ids': '22,23', 'three_prime_feature': 'intron_23', 'read_id': 'read_17', 'expected': 'cryptic'},
        
        # 18. Cryptic: gap between intron 25 and intron 26
        {'splice_count': 2, 'exon_ids': '23,24,25', 'intron_ids': '22,26,27', 'three_prime_feature': 'intron_27', 'read_id': 'read_18', 'expected': 'cryptic'},

        # 19: Wrong 3' end feature (DOG to intron)
        {'splice_count': 2, 'exon_ids': '23,24,25', 'intron_ids': '22,26,27', 'three_prime_feature': 'exon_25', 'read_id': 'read_19', 'expected': 'cryptic'},
        
        # 20. Cryptic: Spans exons and introns but splice_count does not align
        {'splice_count': 4, 'exon_ids': '26,27,28', 'intron_ids': '26,27', 'three_prime_feature': 'exon_28', 'read_id': 'read_20', 'expected': 'cryptic'},
        
        # 21. Canonical: No splicing events, with exons spanned correctly
        {'splice_count': 0, 'exon_ids': '29,30', 'intron_ids': '29', 'three_prime_feature': 'exon_30', 'read_id': 'read_21', 'expected': 'canonical'},

        # 22. Canonical: No splicing events, with exons spanned correctly
        {'splice_count': 1, 'exon_ids': '1,2,4,5', 'intron_ids': '1,2,3', 'three_prime_feature': 'exon_5', 'read_id': 'read_21', 'expected': 'cryptic'},
    ]
    
    df_tests = pd.DataFrame(test_cases)

    # a minimal dictionary:
    #   GENE_01 => {1,2,3} (meaning it has introns #1, #2, #3)
    #   GENE_02 => {1}     (just one intron)
    #   anything else => empty set
    mock_gene_to_introns = {
        'GENE_01': {1, 2, 3},
        'GENE_02': {1},
    }

    # classify 
    classified_df = classify_previous_splicing(
        df_tests, 
        gene_to_introns=mock_gene_to_introns,  
        verbose=True
    )
    
    # verify results
    for i, row in classified_df.iterrows():
        actual_splicing = row['cryptic_splicing']
        expected = row['expected']
        read_id = row['read_id']
        if actual_splicing == expected:
            print(f"Test case {i+1} passed ({read_id}).")
        else:
            reason = row.get('cryptic_reason')
            if pd.isna(reason):
                reason = 'No reason provided.'
            print(f"Test case {i+1} failed ({read_id}): Expected '{expected}', got '{actual_splicing}'. Reason: {reason}")

RUN_TESTS = False 

if RUN_TESTS:
    test_classify_splicing()