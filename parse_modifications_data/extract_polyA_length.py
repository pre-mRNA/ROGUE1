import pysam

def fetch_polyA_pt(bam_file):
    """
    Extracts polyA tail lengths from a BAM file.

    Args:
    bam_file (str): Path to the BAM file.

    Returns:
    dict: A dictionary of read IDs and their polyA tail lengths.
    """
    polyA_lengths = {}

    with pysam.AlignmentFile(bam_file, "rb") as bam:

        for read in bam:

            # check for dorado pt tag 
            if read.has_tag('pt'):
                # convert the polyA length to integer if present
                polyA_lengths[read.query_name] = int(read.get_tag('pt'))
            
            else:
                # store None if the polyA length is not available
                polyA_lengths[read.query_name] = None
    
    return polyA_lengths