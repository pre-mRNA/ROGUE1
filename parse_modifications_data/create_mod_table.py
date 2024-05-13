import pysam
import re
import logging
from Bio.Seq import Seq

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# regex pattern for performance when parsing CIGAR strings
CIGAR_PATTERN = re.compile(r'(\d+)([MIDNSHPX=])')

def parse_cigar(cigar_string):
    try:
        parsed_cigar = [(int(length), op) for length, op in CIGAR_PATTERN.findall(cigar_string)]
        logging.info(f"Parsed CIGAR: {parsed_cigar}")
        return parsed_cigar
    except Exception as e:
        logging.error(f"Error parsing CIGAR string: {cigar_string} | Error: {e}")
        return []

def parse_mm_ml_tags(mm_tag, ml_probs):
    modifications = []
    prob_index = 0
    try:
        for mod in mm_tag.split(';'):
            if not mod:
                continue
            numbers_start = re.search(r'\d', mod)
            if numbers_start:
                elements = mod[numbers_start.start():].split(',')
                base_position = 0
                for skip in elements:
                    base_position += int(skip)
                    if prob_index < len(ml_probs):
                        modifications.append((base_position, ml_probs[prob_index] / 255))
                        prob_index += 1
    except Exception as e:
        logging.error(f"Error parsing MM tag {mm_tag}: {e}")
    return modifications

def process_modifications(read, mm_tag, ml_tag):
    modifications = []
    seq = read.seq
    start_pos = read.reference_start
    cigar = read.cigarstring
    is_reverse = read.is_reverse

    
    mod_positions = parse_mm_ml_tags(mm_tag, ml_probs=list(ml_tag))


    adjusted_positions = adjust_positions_for_orientation(mm_tag, seq, is_reverse)

    print(f"adjusted positions are {adjusted_positions}")

    for pos, prob in adjusted_positions:
        genomic_pos = calculate_genomic_position(pos, cigar, start_pos, is_reverse)
        if genomic_pos is not None:
            modifications.append(f"{genomic_pos}:{prob:.4f}")

    if not modifications:
        logging.info(f"No valid modifications detected for read starting at {read.query_name}")

    return modifications


# return the index of A in the original read sequence 
def adjust_positions_for_orientation(mm_tag, seq, is_reverse):
    adjusted_positions = []

    # revcomp sequences on the minus strand
    if is_reverse:
        seq = str(Seq(seq).reverse_complement())
        print(f"reversed string is {seq}")

    # get indices of A in the original sequence 
    a_positions = [i for i, base in enumerate(seq) if base == 'A']
    print(f"positions are {a_positions}")

    # get the deltas from the MM tag 
    delta_positions = []
    for mod in mm_tag.split(';'):
        if 'A+a' in mod:
            parts = mod.split(',')
            for part in parts[1:]:  
                if part.isdigit():
                    delta_positions.append(int(part))

    # calculate the positions of the As in the read 
    current_index = -1  

    for i, delta in enumerate(delta_positions):
        print(f"i is {i}, delta is {delta}")
        if i == 0:
            current_index += delta + 1  
        else:
            if delta == 0:
                current_index += 1
            else:
                current_index += delta
        
        # check that the index is within the bounds of the positions list
        if current_index < len(a_positions):
            actual_pos = a_positions[current_index]
            if is_reverse:

                # adjust the position for the reverse-complemented sequence
                actual_pos = len(seq) - actual_pos - 1
            adjusted_positions.append((actual_pos, 1))  # Use dummy probability 1, replace as needed
            logging.info(f"Actual position for 'A' at index {current_index}: {actual_pos}")
        else:
            logging.warning(f"'A' at index {current_index} is out of range.")

    return adjusted_positions


def calculate_genomic_position(pos, cigar, start_pos, is_reverse):
   
    current_pos = start_pos # position on reference 
    seq_pos = 0 # position on the read 

    print(f"genomic pos is {current_pos}, distance to target is {pos - seq_pos}")

    for length, operation in parse_cigar(cigar):
        print(f"genomic pos is {current_pos} and distance to target is {pos - seq_pos}")
        print(f"iterating {length} nt due to {operation}")

        if operation in {'M', 'X', '='}:  # match
            print("match")
            if seq_pos <= pos < seq_pos + length:
                print(f"we have a hit at {current_pos + (pos - seq_pos)}")
                return current_pos + (pos - seq_pos) if not is_reverse else current_pos + (length - (pos - seq_pos) - 1)
            
            seq_pos += length
            current_pos += length

        elif operation in {'D', 'N'}:  # deletion / splice 
            print("deletion")    
            current_pos += length

        elif operation in {'I', 'S'}:  
            print("insertion")  
            # if the target position is within the softclip, return none since it's not in the reference 
            if seq_pos <= pos < seq_pos + length:
                print("Modification is softclipped")
                return None  # insertions are not seen in the reference 
            
            # otherwise, move down the sequence
            seq_pos += length
            print(f"seq_pos is {seq_pos}")

    return None  

# main function 
def extract_modifications(bam_file):
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            
            # grab mod tags 
            mm_tag = read.get_tag("MM") if read.has_tag("MM") else None
            ml_tag = read.get_tag("ML") if read.has_tag("ML") else None
            
            if mm_tag and ml_tag:
                modifications = process_modifications(read, mm_tag, ml_tag)
                if modifications:
                    yield modifications
                else:
                    logging.info(f"No modifications found for read {read.query_name}")
            else:
                logging.info("MM or ML tags missing for some reads.")
