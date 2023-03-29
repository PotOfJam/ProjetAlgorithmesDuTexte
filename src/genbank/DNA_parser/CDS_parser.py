import logging
import sequence_parser as sp

def CDSParser(DNA, DNA_length, feature):

    # Find sequence location
    sequence_location = sp.sequenceLocation(feature, DNA_length)

    # Recreate DNA sequence
    DNA_sequence = ""
    if sequence_location == []:
        return
    elif len(sequence_location) == 1:
        sequence_start = sequence_location[0][0]
        sequence_end = sequence_location[0][1]
        DNA_sequence += DNA[sequence_start : sequence_end]
    else:
        DNA_sequence = sp.defragmentSequence(DNA, sequence_location)

    # Check for invalid DNA sequence
    if sp.incorrectSequence(DNA_sequence):
        return

    # Write result in file
    return