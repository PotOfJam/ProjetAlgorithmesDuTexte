import os, logging

def incorrectSequenceLocation(start_location, end_location, DNA_length):

    # Invalid location (not int)
    if ("<" in str(start_location)) or (">" in str(start_location)) or ("<" in str(end_location)) or (">" in str(end_location)):
        logging.error("Invalid sequence start/end location (%s,%s)" % (start_location, end_location))
        return True
    # Incorrect start location
    if start_location < 0:
        logging.error("Incorrect sequence start location (%d,%d) (start_location < 0)" % (start_location, end_location))
        return True
    if end_location > DNA_length:
        logging.error("Incorrect sequence end location (%d,%d) (end_location > DNA_seq_len)" % (start_location, end_location))
        return True
    if end_location <= start_location:
        logging.error("Incorrect sequence start/end location (%d,%d) (end_location <= start_location)" % (start_location, end_location))
        return True
    return False


def sequenceLocation(feature, DNA_length):

    # Initialise variable
    sequence_location = []

    # Fragmented DNA sequence
    if str(feature.location_operator) == "join":
        for part_location in feature.location.parts:
            start_location = part_location.start
            end_location = part_location.end
            if incorrectSequenceLocation(start_location, end_location, DNA_length):
                return []
            sequence_location.append((int(start_location), int(end_location)))
        # Reorder sub-sequence locations
        def getFirst(tuple):
            return tuple[0]
        sequence_location.sort(key=getFirst)
    # Contiguous DNA sequence
    else:
        start_location = part_location.start
        end_location = part_location.end
        if incorrectSequenceLocation(start_location, end_location, DNA_length):
                return []
        sequence_location.append((int(start_location), int(end_location)))
    
    return sequence_location


def defragmentSequence(DNA, sequence_location):

    # Note: sequence locations have been ordered in sequenceLocation function

    # Recreate sequence
    DNA_sequence = ""
    for sub_sequence_location in sequence_location:
        sub_sequence_start = sub_sequence_location[0]
        sub_sequence_end = sub_sequence_location[1]
        DNA_sequence += DNA[sub_sequence_start : sub_sequence_end]

    return DNA_sequence


def incorrectSequence(DNA_sequence):

    # Initialise constants
    valid_start_codon = ["ATG", "CTG", "TTG", "GTG", "ATA", "ATC", "ATT", "TTA"]
    valid_end_codon = ["TAA", "TAG", "TGA"]

    start_codon = DNA_sequence[0:3]
    end_codon = DNA_sequence[len(DNA_sequence)-3:]

    # Invalid start codon
    if start_codon not in valid_start_codon:
        logging.error("Invalid start codon (%s)" % start_codon)
        return True
    #Invalid end codon
    if end_codon not in valid_end_codon:
        logging.error("Invalid end codon (%s)" % end_codon)
        return True
    # Invalid DNA base
    if not all(base in "ATGC" for base in DNA_sequence):
        logging.error("Invalid base in sequence")
        return True

    return False

def writeSequence(sequence_info):

    # File path
    file_path = os.path.join(sequence_info["path"], sequence_info["type"] + "_" + sequence_info["organism"] + "_NC_" + str(sequence_info["NC"]) + ".txt" )
    
    # Sequence description
    sequence_description_text = sequence_info["type"] + " " + sequence_info["organism"] + " " + sequence_info["NC"] + ": "
    if len(sequence_info["locations"]) == 1:
        sequence_description_text += str(sequence_info["location"][0][0]) + ".." + str(sequence_info["location"][0][1])
    else:
        sequence_description_text += "join("
        for sub_sequence_location in sequence_info["locations"]:
            sequence_description_text += str(sub_sequence_location[0]) + ".." + str(sub_sequence_location[1]) + ","
        sequence_description_text -= ","
        sequence_description_text += ")"

    try:
        with open(file_path) as file:
            # Write full sequence
            file.write(sequence_description_text)
            file.write(sequence_info["DNA_sequence"])

            # Write sub-sequences
            exon_id = 0
            for subsequence in sequence_info["DNA_sub_sequence"]:
                exon_id += 1
                file.write(sequence_description_text + " Exon " + str(exon_id))
                file.write(subsequence)
    except:
        logging.error("Unable to write in file: " + file_path)