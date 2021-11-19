"""
I actualy don't know what it does
"""
from Bio.Seq import Seq
from Bio.SeqUtils import GC, molecular_weight


def translate_dna(sequence, frames, reverse=False):
    """ Translates a given DNA sequence in the given frames
    :param sequence: String; a DNA sequence
    :param frames: List; a numeric list given frames (1-3 for forward, 4-6 for reverse)
    :param reverse: Bool; whether or not to reverse-complement the sequence before translation
    :return: Dict; a dictionary with frame number as key and the translated sequence as value
    """
    translated = {}

    sequence += 'N' * (len(sequence) % 3)

    for frame in frames:
        if reverse:
            dna = Seq(sequence[:-(frame-3)]).reverse_complement()
        else:
            dna = Seq(sequence[frame:])
        translated[frame] = str(dna.translate())

    return translated

def generate_random_dna(nseq, length, coding):
    ## TODO: implement functionality (see webapp.generate_dna)
    pass

def mutate_dna(seq, prob):
    ## TODO: implement functionality (see webapp.mutate_dna)
    pass

def get_fasta_stats(filename):
    ## TODO: implement functionality (see webapp.fasta_statistics)
    # Hint: use Biopython (see imports)
    return dict()

def calculate_molecular_weight(sequence):
    ## TODO: implement functionality (see webapp.molecular_weight)
    # Hint: use Biopython (see imports)
    return 0

def calculate_gc_percentage(sequence):
    ## TODO: implement functionality (see webapp.gc_percentage)
    # Hint: use Biopython (see imports)
    return 0
