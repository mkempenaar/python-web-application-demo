from Bio.Seq import Seq


def translate_dna(sequence, frames, reverse=False):
    """ Translates a given DNA sequence in the given frames
    :param sequence: String; a DNA sequence
    :param frames: List; a numeric list given frames (1-3 for forward, 4-6 for reverse)
    :param reverse: Bool; whether or not to reverse-complement the sequence before translation
    :return: Dict; a dictionary with frame number as key and the translated sequence as value
    """
    translated = {}

    for f in frames:
        if reverse:
            dna = Seq(sequence[:-(f-3)]).reverse_complement()
        else:
            dna = Seq(sequence[f:])
        translated[f] = str(dna.translate())

    return translated
