import logging
import sequence_parser as sp

def CDSParser(path, id, organism, DNA, DNA_length, feature):

    sequence_info = {
        "path": path,
        "id": id,
        "organism": organism.replace(" ", "_"),
        "type": "CDS"
    }

    # Find sequence location
    sequence_info["location"] = sp.sequenceLocation(feature, DNA_length)

    # Recreate DNA sequence
    if sequence_info["location"] == []:
        return
    elif len(sequence_info["location"]) == 1:
        sequence_info["start"] = sequence_info["location"][0][0]
        sequence_info["end"] = sequence_info["location"][0][1]
        sequence_info["DNA_sequence"] += DNA[sequence_info["start"] : sequence_info["end"]]
    else:
        sequence_info["DNA_sub_sequence"] = []
        for sub_sequence_location in sequence_info["location"]:
            sequence_info["DNA_sub_sequence"].append(DNA[sub_sequence_location[0] : sub_sequence_location[1]])
        sequence_info["DNA_sequence"] = sp.defragmentSequence(DNA, sequence_info["location"])

    # Check for invalid DNA sequence
    if sp.incorrectSequence(sequence_info["DNA_sequence"]):
        return

    # Write result in file
    sp.writeSequence(sequence_info)