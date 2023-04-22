import logging
from .sequence_parser_utils import *
from ...app.logger import emitLog, Log

def parseSequence(path, file_name, id, organism, DNA, DNA_length, feature, feature_type, worker=None):

    # Create dictionnary containing informations relative to the CDS sequence
    sequence_info = {
        "path": path,
        "file_name": file_name,
        "id": id,
        "organism": organism.replace(" ", "_"),
        "type": feature_type,
        "location": [],
        "DNA_sequence": "",
        "DNA_sub_sequence": [],
        "strand": feature.strand
    }

    # Find sequence location
    sequence_info["location"] = sequenceLocation(feature, DNA_length, worker=worker)

    # Recreate CDS sequence
    if sequence_info["location"] == []:
        emitLog(Log.WARNING, "Incorrect sequence location: (empty location)", worker)
        return
    elif len(sequence_info["location"]) == 1:
        sequence_info["start"] = sequence_info["location"][0][0]
        sequence_info["end"] = sequence_info["location"][0][1]
        emitLog(Log.DEBUG, "location = " + str(sequence_info["start"]) + "," + str(sequence_info["end"]), worker)
        sequence_info["DNA_sequence"] = DNA[sequence_info["start"] : sequence_info["end"]]
    else:
        sequence_info["DNA_sub_sequence"] = []
        for sub_sequence_location in sequence_info["location"]:
            sequence_info["DNA_sub_sequence"].append(DNA[sub_sequence_location[0] : sub_sequence_location[1]])
        sequence_info["DNA_sequence"] = defragmentSequence(DNA, sequence_info["location"], worker=worker)

    # CDS reverse completement
    if sequence_info["strand"] == -1:
        sequence_info["DNA_sequence"] = sequence_info["DNA_sequence"].reverse_complement()
        for sub_sequence in sequence_info["DNA_sub_sequence"]:
            sub_sequence = sub_sequence.reverse_complement()
    # Check for invalid DNA sequence
    if incorrectSequence(sequence_info["DNA_sequence"], sequence_info["type"]):
        emitLog(Log.WARNING, "Incorrect sequence", worker)
        return

    # Write CDS sequence in CDS file
    writeSequence(sequence_info)