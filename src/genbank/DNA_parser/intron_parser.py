import logging
import genbank.DNA_parser.sequence_parser as sp

def IntronParser(path, id, organism, DNA, DNA_length, feature):

    # Create dictionnary containing informations relative to the sequence
    sequence_info = {
        "path": path,
        "id": id,
        "organism": organism.replace(" ", "_"),
        "type": "intron",
        "DNA_sequence": "",
        "DNA_sub_sequence": []
    }
    print(sequence_info)

    # Find sequence location
    sequence_info["location"] = sp.sequenceLocation(feature, DNA_length)

    # Recreate DNA sequence
    if sequence_info["location"] == []:
        logging.warning("Incorrect sequence location: (empty location)")
        return
    elif len(sequence_info["location"]) == 1:
        sequence_info["start"] = sequence_info["location"][0][0]
        sequence_info["end"] = sequence_info["location"][0][1]
        logging.debug("location = " + str(sequence_info["start"]) + "," + str(sequence_info["end"]))
        sequence_info["DNA_sequence"] = DNA[sequence_info["start"] : sequence_info["end"]]
    else:
        sequence_info["DNA_sub_sequence"] = []
        for sub_sequence_location in sequence_info["location"]:
            sequence_info["DNA_sub_sequence"].append(DNA[sub_sequence_location[0] : sub_sequence_location[1]])
        sequence_info["DNA_sequence"] = sp.defragmentSequence(DNA, sequence_info["location"])

    # Reverse completement
    if feature.strand == -1:
        sequence_info["DNA_sequence"] = sequence_info["DNA_sequence"].reverse_complement()
        for sub_sequence in sequence_info["DNA_sub_sequence"]:
            sub_sequence = sub_sequence.reverse_complement()            
    
    # Check for invalid DNA sequence
    if sp.incorrectSequence(sequence_info["DNA_sequence"]):
        logging.warning("Incorrect sequence")
        return

    # Write result in file
    sp.writeSequence(sequence_info)