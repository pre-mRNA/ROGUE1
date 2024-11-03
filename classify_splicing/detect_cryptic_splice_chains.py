import pandas as pd
import numpy as np
import warnings

# parse_ids function to convert comma-separated string of IDs to a list of integers
def parse_ids(id_string):

    if pd.isna(id_string) or id_string.strip() == '':
        return []
    return sorted(int(x) for x in id_string.split(',') if x.strip().isdigit())

# classify_splicing function to classify reads as 'cryptic' or 'canonical' based on exon_ids, intron_ids, splice_count, and three_prime_feature
def classify_splicing(df, verbose=False):
    """
    Classify each read as 'cryptic' or 'canonical' based on exon_ids, intron_ids, splice_count, and three_prime_feature.

    Parameters:
    df (pd.DataFrame): DataFrame containing the necessary columns:
        - 'exon_ids' (str): Comma-separated exon IDs.
        - 'intron_ids' (str): Comma-separated intron IDs.
        - 'splice_count' (int): Number of splice events.
        - 'three_prime_feature' (str): Feature where the read ends ('exon_*', 'intron_*', or 'region_DOG').
        - 'read_id' (str): Unique identifier for the read (used for warnings).
    verbose (bool): If True, include reasons for 'cryptic' classifications.

    Returns:
    pd.DataFrame: DataFrame with additional 'cryptic_splicing' and 'cryptic_reason' columns.
    """
    
    def classify_row(row):
        exon_ids = parse_ids(row['exon_ids'])
        intron_ids = parse_ids(row['intron_ids'])
        splice_count = row['splice_count']
        three_prime_feature = row['three_prime_feature']
        read_id = row['read_id']
        
        expected_splice_count = 0
        
        # sort features 
        all_features = []
        for exon_id in exon_ids:
            all_features.append(('exon', exon_id))
        for intron_id in intron_ids:
            all_features.append(('intron', intron_id))
        all_features = sorted(all_features, key=lambda x: (x[1], x[0] == 'intron'))
        if three_prime_feature == 'region_DOG':
            all_features.append(('region_DOG', float('inf')))
        
        # determine the read 3' feature 
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
            return ('cryptic', reason if verbose else None)
        
        # find the index of the starting feature in the sorted list
        current_index = None
        for idx, (feature_type, feature_id) in enumerate(all_features):
            if feature_type == current_type and feature_id == current_id:
                current_index = idx
                break
        
        if current_index is None:
            reason = 'Read 3\' end feature feature not found in the list of introns or exons.'
            return ('cryptic', reason if verbose else None)
        
        # check if three_prime_feature matches the highest-ranked feature
        highest_feature_type, highest_feature_id = all_features[-1]
        if current_type != highest_feature_type or current_id != highest_feature_id:
            reason = 'Three_prime_feature does not match the highest-ranked feature.'
            return ('cryptic', reason if verbose else None)
        
        # decrement through features to check for cryptic splicing
        while current_index > 0:
            current_type, current_id = all_features[current_index]
            previous_type, previous_id = all_features[current_index - 1]
            
            if current_type == 'exon':
                if previous_type == 'intron':
                    
                    # the previous intron must be the immediately preceding one
                    if previous_id != current_id - 1:
                        reason = f'Exon {current_id} is not immediately preceded by the corresponding intron.'
                        return ('cryptic', reason if verbose else None)
                    current_index -= 1
                elif previous_type == 'exon':
                    
                    # if moving to a lower-ranked exon, increment expected splice count
                    expected_splice_count += 1
                    current_index -= 1
                else:
                    reason = f'Unexpected preceding feature type {previous_type} for exon {current_id}.'
                    return ('cryptic', reason if verbose else None)
            elif current_type == 'intron':
                if previous_type == 'exon':
                    
                    # the previous exon must be the exon of the immediately lower rank
                    if previous_id != current_id:
                        reason = f'Intron {current_id} is not immediately preceded by the corresponding exon.'
                        return ('cryptic', reason if verbose else None)
                    current_index -= 1
                else:
                    
                    # if moving from one intron to another without an intervening exon, it is cryptic
                    reason = f'Intron {current_id} is not immediately preceded by an exon.'
                    return ('cryptic', reason if verbose else None)
            elif current_type == 'region_DOG':
                
                # DOG regions should only be preceded by the highest-ranked exon
                if previous_type != 'exon' or previous_id != max(exon_ids):
                    reason = 'region_DOG is not properly preceded by the highest-ranked exon.'
                    return ('cryptic', reason if verbose else None)
                current_index -= 1
            
            else:
                reason = 'Unknown feature type during backtracking.'
                return ('cryptic', reason if verbose else None)
        
        # compare expected splice count with actual splice count
        if expected_splice_count != splice_count:
            reason = f'Splice count {splice_count} does not match expected {expected_splice_count}.'
            return ('cryptic', reason if verbose else None)
        
        return ('canonical', None)
    
    # classify each read 
    classification = df.apply(classify_row, axis=1, result_type='expand')
    classification.columns = ['cryptic_splicing', 'cryptic_reason']
    
    # assign to original df 
    df = pd.concat([df, classification], axis=1)
    
    # if verbose is False, drop the 'cryptic_reason' column
    if not verbose:
        df = df.drop(columns=['cryptic_reason'])
    
    return df

# test cases for the classify_splicing function
def test_classify_splicing():
    """
    Test the classify_splicing function with a variety of test cases to ensure correctness.
    """
    test_cases = [

        # 1. Canonical: Proper splicing with exons and introns
        {'splice_count': 1, 'exon_ids': '2,3', 'intron_ids': '1', 'three_prime_feature': 'exon_3', 'read_id': 'read_1', 'expected': 'canonical'},
        
        # 2. Cryptic: splice_count exceeds possible based on exons
        {'splice_count': 3, 'exon_ids': '1,2', 'intron_ids': '', 'three_prime_feature': 'exon_2', 'read_id': 'read_2', 'expected': 'cryptic'},
        
        # 3. Cryptic: splice_count >0 but no exons or introns spanned in DOG region 
        {'splice_count': 1, 'exon_ids': '', 'intron_ids': '', 'three_prime_feature': 'region_DOG', 'read_id': 'read_3', 'expected': 'cryptic'},
        
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
        {'splice_count': 1, 'exon_ids': '1,2,4,5', 'intron_ids': '1,2,3', 'three_prime_feature': 'exon_5', 'read_id': 'read_21', 'expected': 'canonical'},
    ]
    
    # test df
    df_tests = pd.DataFrame(test_cases)
    
    # classify with verbosity 
    classified_df = classify_splicing(df_tests, verbose=True)
    
    # verify results
    for i, row in classified_df.iterrows():
        actual = row['cryptic_splicing']
        expected = row['expected']
        read_id = row['read_id']
        if actual == expected:
            print(f"Test case {i+1} passed.")
        else:
            reason = row.get('cryptic_reason')
            if pd.isna(reason):
                reason = 'No reason provided.'
            print(f"Test case {i+1} failed: Expected '{expected}', got '{actual}'. Reason: {reason}")

# run test 
test_classify_splicing()
