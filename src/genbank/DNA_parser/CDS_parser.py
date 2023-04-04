import os, logging, traceback
import genbank.DNA_parser.sequence_parser_utils as spu

def parseCDS(path, id, organism, DNA, DNA_length, feature):

    # Create dictionnary containing informations relative to the CDS sequence
    CDS_info = {
        "path": path,
        "id": id,
        "organism": organism.replace(" ", "_"),
        "type": "CDS",
        "location": [],
        "DNA_sequence": "",
        "DNA_sub_sequence": [],
        "strand": feature.strand
    }

    # Create dictionnary containing informations relative to the intron(s) sequence(s)
    intron_flag = False
    intron_info = {
        "path": path,
        "id": id,
        "organism": organism.replace(" ", "_"),
        "type": "intron",
        "location": [],
        "CDS_location": [],
        "DNA_sub_sequence": [],
        "strand": feature.strand
    }

    # Find sequence location
    CDS_info["location"] = spu.sequenceLocation(feature, DNA_length)
    intron_info["CDS_location"] = CDS_info["location"]

    # Find intron(s) location(s)
    intron_info["location"] = spu.intronLocation(CDS_info["location"])

    # Recreate CDS sequence
    if CDS_info["location"] == []:
        logging.warning("Incorrect sequence location: (empty location)")
        return
    elif len(CDS_info["location"]) == 1:
        CDS_info["start"] = CDS_info["location"][0][0]
        CDS_info["end"] = CDS_info["location"][0][1]
        logging.debug("location = " + str(CDS_info["start"]) + "," + str(CDS_info["end"]))
        CDS_info["DNA_sequence"] = DNA[CDS_info["start"] : CDS_info["end"]]
    else:
        intron_flag = True
        CDS_info["DNA_sub_sequence"] = []
        for sub_sequence_location in CDS_info["location"]:
            CDS_info["DNA_sub_sequence"].append(DNA[sub_sequence_location[0] : sub_sequence_location[1]])
        CDS_info["DNA_sequence"] = spu.defragmentSequence(DNA, CDS_info["location"])

    # CDS reverse completement
    if CDS_info["strand"] == -1:
        CDS_info["DNA_sequence"] = CDS_info["DNA_sequence"].reverse_complement()
        for sub_sequence in CDS_info["DNA_sub_sequence"]:
            sub_sequence = sub_sequence.reverse_complement()
    # Check for invalid DNA sequence
    if spu.incorrectSequence(CDS_info["DNA_sequence"]):
        logging.warning("Incorrect sequence")
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
    writeCDS(CDS_info)

    # Write introns sequence in intron file
    if intron_flag:
        writeIntron(intron_info)


def writeCDS(CDS_info):

    # File path
    try:
        file_path = os.path.join(CDS_info["path"], CDS_info["type"] + "_" + CDS_info["organism"] + "_" + str(CDS_info["id"]) + ".txt" )
    except:
        logging.debug("PROBLEME 1")

    # Sequence description
    try:
        sequence_description_text = CDS_info["type"] + " " + CDS_info["organism"] + " " + CDS_info["id"] + ": "
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
        logging.debug("PROBLEME 2")

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
        print(traceback.format_exc())
        logging.error("Unable to write in file: " + file_path)


def writeIntron(intron_info):

    # File path
    try:
        file_path = os.path.join(intron_info["path"], intron_info["type"] + "_" + intron_info["organism"] + "_" + str(intron_info["id"]) + ".txt" )
    except:
        logging.debug("PROBLEME 1")

    # Sequence description
    try:
        sequence_description_text = intron_info["type"] + " " + intron_info["organism"] + " " + intron_info["id"] + ": "
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
        logging.debug("PROBLEME 2")

    try:
        with open(file_path, "a") as file:
            # Write sub-sequences
            intron_id = 0
            for subsequence in intron_info["DNA_sub_sequence"]:
                intron_id += 1
                file.write(sequence_description_text + " Intron " + str(intron_id) + "\n")
                file.write(str(subsequence) + "\n")
    except:
        print(traceback.format_exc())
        logging.error("Unable to write in file: " + file_path)