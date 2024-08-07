import pandas as pd
import numpy as np
import warnings

def classify_splicing(df):
    def parse_ids(id_string):
        return set(int(x) for x in id_string.split(',') if x.strip().isdigit()) if pd.notna(id_string) else set()

    def classify_row(row):
        exons = parse_ids(row['exon_ids'])
        introns = parse_ids(row['intron_ids'])
        splice_count = row['splice_count']
        read_id = row['read_id']

        # handling reads where the alignment does not contain a mapped splice-junction 
        if splice_count == 0:

            # no exons, no introns, no splicing
            if len(exons) == 0 and len(introns) == 0:
                return 'ambiguous' 
            
            # single exon and no splicing 
            elif len(exons) == 1 and len(introns) == 0:
                return 'ambiguous'
            
            # single intron, single exon, no splicing 
            elif len(introns) == 1 and len(exons) <= 1:

                # intron is directly upstream of exon and the read is unspliced 
                # we span donor
                # the read is 'fully unspliced'
                if len(exons) == 1 and list(introns)[0] + 1 == list(exons)[0]:
                    return 'fully-unspliced'
                else:

                    # in any other case, we span an intron and exon but we don't span the intron 3' ss 
                    # so the read is ambigous 
                    return 'ambiguous'
            
                
            elif any(i + 1 in exons for i in introns):
                return 'fully-unspliced' # intron running into exon, therefore unspliced 
            

            elif len(exons) > 1 and len(introns) == 0:
                return 'ambiguous'  # multiple exons without splicing or introns 
            else:
                warnings.warn(f"Read ID: {read_id}, Unexpected case: splice_count=0, exons={exons}, introns={introns}")
                return 'ambiguous'  
        elif splice_count > 0:
            if len(exons) >= 2:
                return 'spliced'
            elif len(exons) == 1 or (len(exons) == 0 and len(introns) == 0):
                return 'ambiguous'  # splicing event that doesn't span 2 annotated exons or other noise 
            else:
                warnings.warn(f"Read ID: {read_id}, Unexpected case: splice_count>0, exons={exons}, introns={introns}")
                return 'ambiguous'  
        else:
            warnings.warn(f"Read ID: {read_id}, Unexpected splice_count: {splice_count}")
            return 'ambiguous'  

    df['splicing_classification'] = df.apply(classify_row, axis=1)
    return df

# simple unit tests 
def test_classify_splicing():
    test_cases = [
        {'read_id': '0', 'splice_count': 0, 'exon_ids': '1', 'intron_ids': '1', 'expected': 'ambiguous'},
        {'read_id': '1', 'splice_count': 0, 'exon_ids': '2', 'intron_ids': '1', 'expected': 'fully-unspliced'},
        # {'read_id': '1', 'splice_count': 0, 'exon_ids': '1', 'intron_ids': '', 'expected': 'ambiguous'}, # single exon with no splicing or intron span, therefore ambiguous
        # {'read_id': '2', 'splice_count': 1, 'exon_ids': '1,2,3', 'intron_ids': '1', 'expected': 'spliced'}, # intron 1 unspliced, intron 2 spliced, therefore, a spliced event
        # {'read_id': '3', 'splice_count': 0, 'exon_ids': '', 'intron_ids': '1', 'expected': 'ambiguous'}, # single intron with no splicing or exon span, therefore ambiguous
        # {'read_id': '4', 'splice_count': 0, 'exon_ids': '2', 'intron_ids': '1', 'expected': 'fully-unspliced'}, # spans the intron donor with no evidence of splicing, therefore unspliced 
        # {'read_id': '5', 'splice_count': 1, 'exon_ids': '1,2', 'intron_ids': '', 'expected': 'spliced'}, # two adjacent exons linked by a single junction, therefore spliced 
        # {'read_id': '6', 'splice_count': 0, 'exon_ids': '1,2', 'intron_ids': '', 'expected': 'ambiguous'}, # two adjacent exons with no splicing or intron span, therefore ambiguous
        # {'read_id': '7', 'splice_count': 1, 'exon_ids': '1', 'intron_ids': '', 'expected': 'ambiguous'}, # splicing event that doesn't span 2 exons   
        # {'read_id': '8', 'splice_count': 0, 'exon_ids': '', 'intron_ids': '', 'expected': 'ambiguous'}, # no data 
        # {'read_id': '9', 'splice_count': 1, 'exon_ids': '', 'intron_ids': '', 'expected': 'ambiguous'}, # splicing event that doesn't span 2 exons 
        # {'read_id': '10', 'splice_count': 0, 'exon_ids': '2', 'intron_ids': '3', 'expected': 'ambiguous'}, # doesn't span intron acceptor, therefore ambiguous
    ]
    for i, case in enumerate(test_cases):
        df = pd.DataFrame([case])
        result = classify_splicing(df)
        expected = case['expected']
        actual = result.loc[0, 'splicing_classification']
        
        if pd.isna(expected) and pd.isna(actual):
            print(f"Test case {i+1} passed")
        elif expected == actual:
            print(f"Test case {i+1} passed")
        else:
            raise AssertionError(f"Test case {i+1} failed: expected {expected}, got {actual}")


test_classify_splicing()