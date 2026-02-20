import os
import gzip

#Functions from class 3

#GFF READER
from functools import total_ordering

@total_ordering
class GffEntry:
    """Main class for handling GFF entries
    
    Truncates the GFF entry to required data only. Also, class is totally ordered.
    This means that comparison operators can be used on GffEntry objects
    
    Attributes:
        seqid (str): The contig that the GFF entry is associated with
        type (str): The type of object that the GFF entry classifies as
        start (int): The left-most nucleotide position of the GFF entry relative to seqid (1-indexed)
        end (int): The right-most nucleotide position of the GFF entry relative to seqid (1-indexed)
        strand (str): Whether the entry is on the forward (+), backward (-) strand or N/A (.)
    """
    
    slots = 'seqid type start end strand'.split()
    
    def __init__(self, args):
        """Initialize the object.
        
        Aggregates all GFF entry columns, and selectively assigns them to attributes
        
        Args:
            args (list): the complete stripped and split GFF entry line
        """
        self.seqid = args[0]
        self.type = args[2]
        self.start = int(args[3])
        self.end = int(args[4])
        self.strand = args[6]
    
    def __str__(self):
        """Determines how GffEntry appear when `print()` is called on them"""
        return f'{self.seqid}\t{self.type}\t{self.start}\t{self.end}\t{self.strand}'
    
    def __len__(self):
        """Determines how GffEntry reports when `len()` is called on them"""
        return self.end - self.start

    def __eq__(self, other):
        self_check = (self.seqid, self.type, self.start, self.end, self.strand)
        other_check = (other.seqid, other.type, other.start, other.end, other.strand)
        if  self_check == other_check:
            return True
    
    def __lt__(self, other):
        if self.seqid < other.seqid:
            return True
        elif self.seqid == other.seqid:
            if self.start < other.start:
                return True
            elif self.start == other.start:
                if self.end < other.end:
                    return True
        

def get_gff(gff_file):
    """Generator that lazily reports each of the GFF entries within the GFF file
    
    Args:
        gff_file (str): /path/and/name/to.gff[.gz]
    
    Yields:
        (GffEntry): A GFF entry object with the attributes of seqid, type, start, end, and strand
    """
    if '.gz' in gff_file:
        gff_file = gzip.open(gff_file, 'rb')
    else:
        gff_file = open(gff_file, 'rb')
    for entry in gff_file:
        entry = entry.decode('ascii')
        if entry.startswith('#') or not entry:
            continue
        yield GffEntry(entry.strip().split('\t'))

#FASTA READER

def get_fasta(file):
    """Generator to lazily get all the fasta entries from a fasta file

    Args:
        fasta_file (str): /path/and/name/to.fasta

    Yields:
        header (str): header sequence of fasta entry (includes '>')
        seq (str): concatenated string sequence of the fasta entry
    """
    
    if '.gz' in file:
        file_object = gzip.open(file, 'rt')
    else:
        file_object = open(file, 'rt')
        
    name=''
    seq=''
    for line in file_object:
        
        # Capture the next header, report what we have, and update
        if line.startswith('>') and seq: #not first seq
            name = name[1:] #removes the carrot
            yield name, seq
            name=line.strip()
            seq=''
            
        # Get to the first header
        elif line.startswith('>'):  #first seq
            name=line.strip()
            
        # Just add sequence if it is the only thing there
        else:
            seq+=line.strip()
            
    # At the end, return the last entries
    if name and seq: #last seq
            name = name[1:]
            yield name, seq
