import os, logging, traceback
from .sequence_parser_utils import *
from ...app.logger import emitLog, Log

def parseCDS(path, file_name, id, organism, DNA, DNA_length, feature, CDS_flag, intron_flag, worker=None):

    # Create dictionnary containing informations relative to the CDS sequence
    CDS_info = {
        "path": path,
        "file_name": file_name,
        "id": id,
        "organism": organism.replace(" ", "_"),
        "type": "CDS",
        "location": [],
        "DNA_sequence": "",
        "DNA_sub_sequence": [],
        "strand": feature.strand
    }

    # Create dictionnary containing informations relative to the intron(s) sequence(s)
    intron_info = {
        "path": path,
        "file_name": file_name,
        "id": id,
        "organism": organism.replace(" ", "_"),
        "type": "intron",
        "location": [],
        "CDS_location": [],
        "DNA_sub_sequence": [],
        "strand": feature.strand
    }

    # Find sequence location
    CDS_info["location"] = sequenceLocation(feature, DNA_length, worker=worker)
    intron_info["CDS_location"] = CDS_info["location"]

    # Find intron(s) location(s)
    if intron_flag:
        intron_info["location"] = intronLocation(CDS_info["location"], worker=worker)

    # Recreate CDS sequence
    if CDS_info["location"] == []:
        emitLog(Log.ERROR, "Incorrect sequence location: (empty location)", worker)
        return
    elif len(CDS_info["location"]) == 1:
        intron_flag = False
        CDS_info["start"] = CDS_info["location"][0][0]
        CDS_info["end"] = CDS_info["location"][0][1]
        emitLog(Log.INFO, "location = " + str(CDS_info["start"]) + "," + str(CDS_info["end"]), worker)
        CDS_info["DNA_sequence"] = DNA[CDS_info["start"] : CDS_info["end"]]
    else:
        CDS_info["DNA_sub_sequence"] = []
        for sub_sequence_location in CDS_info["location"]:
            CDS_info["DNA_sub_sequence"].append(DNA[sub_sequence_location[0] : sub_sequence_location[1]])
        CDS_info["DNA_sequence"] = defragmentSequence(DNA, CDS_info["location"], worker=worker)

    # CDS reverse completement
    if CDS_info["strand"] == -1:
        CDS_info["DNA_sequence"] = CDS_info["DNA_sequence"].reverse_complement()
        for sub_sequence in CDS_info["DNA_sub_sequence"]:
            sub_sequence = sub_sequence.reverse_complement()
    # Check for invalid DNA sequence
    if incorrectSequence(CDS_info["DNA_sequence"], "CDS", worker=worker):
        emitLog(Log.ERROR, "Incorrect sequence", worker)
        return
    
    # Recreate intron(s) sequence(s)
    if intron_flag:
        for start_location, end_location in intron_info["location"]:
            intron_info["DNA_sub_sequence"].append(DNA[start_location : end_location])
    
        # Intron reverse completement
        if CDS_info["strand"] == -1:
            for sub_sequence in intron_info["DNA_sub_sequence"]:
                sub_sequence = sub_sequence.reverse_complement()

    # Write CDS sequence in CDS file
    if CDS_flag:
        writeCDS(CDS_info, worker=worker)

    # Write introns sequence in intron file
    if intron_flag:
        writeIntron(intron_info, worker=worker)


def writeCDS(CDS_info, worker=None):
    """
    Write CDS in result file.

    Args:
        CDS_info (dict): Information relative to the CDS.
        worker (_type_, optional): _description_. Defaults to None.
    """
    # File path
    try:
        file_path = os.path.join(CDS_info["path"], CDS_info["type"] + "_" + CDS_info["organism"] + "_" + CDS_info["file_name"] + ".txt" )
    except:
        emitLog(Log.ERROR, "Invalid file path", worker)

    # Sequence description
    try:
        sequence_description_text = CDS_info["type"] + " " + CDS_info["organism"].replace("_", " ") + " " + CDS_info["file_name"] + ": "
        if CDS_info["strand"] == -1:
            sequence_description_text += "complement("
        if len(CDS_info["location"]) == 1:
            sequence_description_text += str(CDS_info["location"][0][0]) + ".." + str(CDS_info["location"][0][1])
        else:
            sequence_description_text += "join("
            for sub_sequence_location in CDS_info["location"]:
                sequence_description_text += str(sub_sequence_location[0]) + ".." + str(sub_sequence_location[1]) + ","
            sequence_description_text = sequence_description_text[:len(sequence_description_text)-1]
            sequence_description_text += ")"
        if CDS_info["strand"] == -1:
            sequence_description_text += ")"
    except:
        emitLog(Log.ERROR, "Invalid sequence description", worker)

    try:
        with open(file_path, "a") as file:
            # Write full sequence
            if CDS_info["DNA_sequence"] != "":
                file.write(sequence_description_text + "\n")
                file.write(str(CDS_info["DNA_sequence"]) + "\n")

            # Write sub-sequences
            exon_id = 0
            for subsequence in CDS_info["DNA_sub_sequence"]:
                exon_id += 1
                file.write(sequence_description_text + " Exon " + str(exon_id) + "\n")
                file.write(str(subsequence) + "\n")
    except:
        emitLog(Log.ERROR, traceback.format_exc(), worker)
        emitLog(Log.ERROR, "Unable to write in file: " + file_path, worker)


def writeIntron(intron_info, worker=None):
    """
    Write intron in file.

    Args:
        intron_info (dict): Information relative to the intron.
        worker (_type_, optional): _description_. Defaults to None.
    """
    # File path
    try:
        file_path = os.path.join(intron_info["path"], intron_info["type"] + "_" + intron_info["organism"] + "_" + intron_info["file_name"] + ".txt" )
    except:
        emitLog(Log.ERROR, "Invalid file path", worker)

    # Sequence description
    try:
        sequence_description_text = intron_info["type"] + " " + intron_info["organism"].replace("_", " ") + " " + intron_info["file_name"] + ": "
        if intron_info["strand"] == -1:
            sequence_description_text += "complement("
        sequence_description_text += "join("
        for sub_sequence_location in intron_info["CDS_location"]:
            sequence_description_text += str(sub_sequence_location[0]) + ".." + str(sub_sequence_location[1]) + ","
        sequence_description_text = sequence_description_text[:len(sequence_description_text)-1]
        sequence_description_text += ")"
        if intron_info["strand"] == -1:
            sequence_description_text += ")"
    except:
        emitLog(Log.ERROR, "Invalid sequence description", worker)

    try:
        with open(file_path, "a") as file:
            # Write sub-sequences
            intron_id = 0
            for subsequence in intron_info["DNA_sub_sequence"]:
                intron_id += 1
                file.write(sequence_description_text + " Intron " + str(intron_id) + "\n")
                file.write(str(subsequence) + "\n")
    except:
        emitLog(Log.ERROR, traceback.format_exc(), worker)
        emitLog(Log.ERROR, "Unable to write in file: " + file_path, worker)