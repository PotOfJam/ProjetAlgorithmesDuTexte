import logging
from .DNA_parser.CDS_parser import parseCDS
from .DNA_parser.sequence_parser import parseSequence
from ..app.parser_thread import emitLog

def parseFeatures(region_type, path, id, organism, record, worker=None):
    """
    Parse a feature found in a record (GenBank file).

    Args:
        region_type (list): Types of DNA region to parse.
        path (string): Path of the organism in the "Result" folder.
        id (string): ID of the GenBank file.
        organism (string): Name of the organism.
        record (...): Record (content of the GenBank file).
    """
    emitLog(worker, "Start parsing id = " + str(id))

    # Initialize variables
    file_name = id.split(".")[0]
    DNA = None
    DNA_length = -1

    # Read DNA sequence from record
    try:
        DNA = record.seq
        DNA_length = len(DNA)
    except:
        emitLog(worker, "Unable to read DNA sequence from id = " + str(id))

    # Parse feature according to its type
    try:
        for feature in record.features:
            # if feature.type in region_type:
            if feature.type == "CDS" and (("CDS" in region_type) or ("intron" in region_type)):
                emitLog(worker, "Found CDS or intron in id = " + str(id))
                parseCDS(path, file_name, id, organism, DNA, DNA_length, feature, "CDS" in region_type, "intron" in region_type, worker=worker)
            elif feature.type == "centromere" and "centromere" in region_type:
                emitLog(worker, "Found centromere in id = " + str(id))
                parseSequence(path, file_name, id, organism, DNA, DNA_length, feature, "centromere", worker=worker)
            elif feature.type == "mobile_element" and "mobile_element" in region_type:
                emitLog(worker, "Found mobile_element in id = " + str(id))
                parseSequence(path, file_name, id, organism, DNA, DNA_length, feature, "mobile_element", worker=worker)
            elif feature.type == "telomere" and "telomere" in region_type:
                emitLog(worker, "Found telomere in id = " + str(id))
                parseSequence(path, file_name, id, organism, DNA, DNA_length, feature, "telomere", worker=worker)
            elif feature.type == "3'UTR" and "3'UTR" in region_type:
                emitLog(worker, "Found 3'UTR in id = " + str(id))
                parseSequence(path, file_name, id, organism, DNA, DNA_length, feature, "3'UTR", worker=worker)
            elif feature.type == "5'UTR" and "5'UTR" in region_type:
                emitLog(worker, "Found 5'UTR in id = " + str(id))
                parseSequence(path, file_name, id, organism, DNA, DNA_length, feature, "5'UTR", worker=worker)
            elif feature.type == "tRNA" and "tRNA" in region_type:
                emitLog(worker, "Found tRNA in id = " + str(id))
                parseSequence(path, file_name, id, organism, DNA, DNA_length, feature, "tRNA", worker=worker)
            elif feature.type == "rRNA" and "rRNA" in region_type:
                emitLog(worker, "Found rRNA in id = " + str(id))
                parseSequence(path, file_name, id, organism, DNA, DNA_length, feature, "rRNA", worker=worker)
            elif feature.type == "ncRNA" and "ncRNA" in region_type:
                emitLog(worker, "Found ncRNA in id = " + str(id))
                parseSequence(path, file_name, id, organism, DNA, DNA_length, feature, "ncRNA", worker=worker)
    except:
        emitLog(worker, "Unable to read feature(s) from id = " + str(id))