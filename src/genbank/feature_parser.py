import logging
from genbank.DNA_parser.CDS_parser import parseCDS
from genbank.DNA_parser.sequence_parser import parseSequence

def parseFeatures(region_type, path, id, organism, record):
    """
    Parse a feature found in a record (GenBank file).

    Args:
        region_type (list): Types of DNA region to parse.
        path (string): Path of the organism in the "Result" folder.
        id (string): ID of the GenBank file.
        organism (string): Name of the organism.
        record (...): Record (content of the GenBank file).
    """
    logging.info("Start parsing id = " + str(id))

    # Initialize variables
    file_name = ""
    DNA = None
    DNA_length = -1

    # Read DNA sequence from record
    try:
        file_name = record.id
        logging.debug("file_name = " + file_name)
        DNA = record.seq
        DNA_length = len(DNA)
    except:
        logging.error("Unable to read DNA sequence from id = " + str(id))

    # Parse feature according to its type
    try:
        for feature in record.features:
            # if feature.type in region_type:
            if feature.type == "CDS" and (("CDS" in region_type) or ("intron" in region_type)):
                logging.info("Found CDS in id = " + str(id))
                parseCDS(path, file_name, id, organism, DNA, DNA_length, feature, "CDS" in region_type, "intron" in region_type)
            elif feature.type == "centromere":
                logging.info("Found centromere in id = " + str(id))
                parseSequence(path, file_name, id, organism, DNA, DNA_length, feature, "centromere")
            elif feature.type == "mobile_element":
                logging.info("Found mobile_element in id = " + str(id))
                parseSequence(path, file_name, id, organism, DNA, DNA_length, feature, "mobile_element")
            elif feature.type == "telomere":
                logging.info("Found telomere in id = " + str(id))
                parseSequence(path, file_name, id, organism, DNA, DNA_length, feature, "telomere")
            elif feature.type == "3'UTR":
                logging.info("Found 3'UTR in id = " + str(id))
                parseSequence(path, file_name, id, organism, DNA, DNA_length, feature, "3'UTR")
            elif feature.type == "5'UTR":
                logging.info("Found 5'UTR in id = " + str(id))
                parseSequence(path, file_name, id, organism, DNA, DNA_length, feature, "5'UTR")
            elif feature.type == "tRNA":
                logging.info("Found tRNA in id = " + str(id))
                parseSequence(path, file_name, id, organism, DNA, DNA_length, feature, "tRNA")
            elif feature.type == "rRNA":
                logging.info("Found rRNA in id = " + str(id))
                parseSequence(path, file_name, id, organism, DNA, DNA_length, feature, "rRNA")
            elif feature.type == "ncRNA":
                logging.info("Found ncRNA in id = " + str(id))
                parseSequence(path, file_name, id, organism, DNA, DNA_length, feature, "ncRNA")
    except:
        logging.error("Unable to read feature(s) from id = " + str(id))