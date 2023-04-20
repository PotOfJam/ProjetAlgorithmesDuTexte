import os, logging, traceback
from ...app.logger import emitLog, Log

def incorrectSequenceLocation(start_location, end_location, DNA_length, worker=None):
    """
    Check if the sequence location is correct.

    Args:
        start_location (_type_): _description_
        end_location (_type_): _description_
        DNA_length (_type_): _description_
        worker (_type_, optional): _description_. Defaults to None.

    Returns:
        _type_: _description_
    """
    # Invalid location (not int)
    if ("<" in str(start_location)) or (">" in str(start_location)) or ("<" in str(end_location)) or (">" in str(end_location)):
        emitLog(Log.WARNING, "Invalid sequence start/end location (%s,%s)" % (start_location, end_location), worker)
        return True
    # Incorrect start location
    if start_location < 0:
        emitLog(Log.WARNING, "Incorrect sequence start location (%d,%d) (start_location < 0)" % (start_location, end_location), worker)
        return True
    if end_location > DNA_length:
        emitLog(Log.WARNING, "Incorrect sequence end location (%d,%d) (end_location > DNA_seq_len)" % (start_location, end_location), worker)
        return True
    if end_location <= start_location:
        emitLog(Log.WARNING, "Incorrect sequence start/end location (%d,%d) (end_location <= start_location)" % (start_location, end_location), worker)
        return True
    return False


def sequenceLocation(feature, DNA_length, worker=None):
    """


    Args:
        feature (_type_): _description_
        DNA_length (_type_): _description_
        worker (_type_, optional): _description_. Defaults to None.

    Returns:
        _type_: _description_
    """
    # Initialise variable
    sequence_location = []

    try:
        # Fragmented DNA sequence
        if str(feature.location_operator) == "join":
            for part_location in feature.location.parts:
                start_location = part_location.start
                end_location = part_location.end
                emitLog(Log.INFO, "location = " + str(start_location) + "," + str(end_location), worker)
                if incorrectSequenceLocation(start_location, end_location, DNA_length, worker=worker):
                    return []
                sequence_location.append((int(start_location), int(end_location)))
            # Reorder sub-sequence locations
            def getFirst(tuple):
                return tuple[0]
            sequence_location.sort(key=getFirst)
        # Contiguous DNA sequence
        else:
            start_location = feature.location.start
            end_location = feature.location.end
            emitLog(Log.INFO, "location = " + str(start_location) + "," + str(end_location), worker)
            if incorrectSequenceLocation(start_location, end_location, DNA_length, worker=worker):
                return []
            sequence_location.append((int(start_location), int(end_location)))
    except:
        emitLog(Log.ERROR, "Unable to find sequence location", worker)
    
    return sequence_location


def intronLocation(CDS_info_location):

    intron_location = []

    if len(CDS_info_location) > 1:
        for i in range(len(CDS_info_location) - 1):
            start_location = CDS_info_location[i][1] + 1
            end_location = CDS_info_location[i+1][0] - 1
            intron_location.append((start_location, end_location))
    
    return intron_location


def defragmentSequence(DNA, sequence_location, worker=None):

    # Note: sequence locations have been ordered in sequenceLocation function

    # Recreate sequence
    DNA_sequence = ""
    for sub_sequence_location in sequence_location:
        sub_sequence_start = sub_sequence_location[0]
        sub_sequence_end = sub_sequence_location[1]
        DNA_sequence += DNA[sub_sequence_start : sub_sequence_end]

    return DNA_sequence


def incorrectSequence(DNA_sequence, sequence_type, worker=None):

    # Initialise constants
    valid_start_codon = ["ATG", "CTG", "TTG", "GTG", "ATA", "ATC", "ATT", "TTA"]
    valid_end_codon = ["TAA", "TAG", "TGA"]

    start_codon = DNA_sequence[0:3]
    end_codon = DNA_sequence[len(DNA_sequence)-3:]
    
    # Invalid start codon
    if sequence_type in ["CDS"] and start_codon not in valid_start_codon:
        emitLog(Log.WARNING, "Invalid start codon (%s)" % start_codon, worker)
        return True
    # Invalid end codon
    if sequence_type in ["CDS"] and end_codon not in valid_end_codon:
        emitLog(Log.WARNING, "Invalid end codon (%s)" % end_codon, worker)
        return True
    # Invalid length
    if sequence_type in ["CDS"] and len(DNA_sequence) % 3 != 0:
        emitLog(Log.WARNING, "Invalid sequence length", worker)
        return True
    # Invalid DNA base
    if not all(base in "ATGC" for base in DNA_sequence):
        emitLog(Log.WARNING, "Invalid base in sequence", worker)
        return True

    return False

def writeSequence(sequence_info, worker=None):

    # File path
    try:
        file_path = os.path.join(sequence_info["path"], sequence_info["type"] + "_" + sequence_info["organism"] + "_" + sequence_info["file_name"] + ".txt" )
    except:
        emitLog(Log.ERROR, "Invalid file path", worker)

    # Sequence description
    try:
        sequence_description_text = sequence_info["type"] + " " + sequence_info["organism"].replace("_", " ") + " " + sequence_info["file_name"] + ": "
        if sequence_info["strand"] == -1:
            sequence_description_text += "complement("
        if len(sequence_info["location"]) == 1:
            sequence_description_text += str(sequence_info["location"][0][0]) + ".." + str(sequence_info["location"][0][1])
        else:
            sequence_description_text += "join("
            for sub_sequence_location in sequence_info["location"]:
                sequence_description_text += str(sub_sequence_location[0]) + ".." + str(sub_sequence_location[1]) + ","
            sequence_description_text = sequence_description_text[:len(sequence_description_text)-1]
            sequence_description_text += ")"
        if sequence_info["strand"] == -1:
            sequence_description_text += ")"
    except:
        emitLog(Log.ERROR, "Invalid sequence description", worker)

    try:
        # Write full sequence
        with open(file_path, "a") as file:
            if sequence_info["DNA_sequence"] != "":
                file.write(sequence_description_text + "\n")
                file.write(str(sequence_info["DNA_sequence"]) + "\n")
    except:
        emitLog(Log.ERROR, traceback.format_exc(), worker)
        emitLog(Log.ERROR, "Unable to write in file: " + file_path, worker)