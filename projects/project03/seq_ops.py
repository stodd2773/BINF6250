COMP_TRANS = str.maketrans('ACGTacgt', 'TGCAtcga')

def reverse_complement(seq):
    """Get the reverse complement of a nucleotide sequence

    Returns the reverse complement of the input string representing a DNA 
    sequence. Works only with DNA sequences consisting solely of  A, C, G, T or N 
    characters. Preserves the case of the input sequence.

    Args:
        seq (str): a DNA sequence string

    Returns:
        (str): The reverse complement of the input DNA sequence string.
    """

    # Translate then reverse seq
    return seq.translate(COMP_TRANS)[::-1]

def get_seq(seq, start, end, strand, size):
    """Get the desired sequence
    
    Args:
        seq (str): nucleotide sequence
        start (int): left-most desired nucleotide sequence
        end (int): right-most desired nucleotide sequence
        strand (str): Whether the entry is on the forward (+), backward (-) strand or N/A (.)
        size (int): how far upstream to get extra sequence (default: 50)
    
    Returns:
        promoter_seq (str): the desired GFF entry and its upstream promoter sequence
    """
    promoter_seq = ''

    if strand == '-':
        promoter_seq = reverse_complement(seq[end:(end+size)])
    else :
        promoter_seq = seq[(start-size-1):(start-1)]

    return promoter_seq



