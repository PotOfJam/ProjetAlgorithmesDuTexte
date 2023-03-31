import logging
from genbank.DNA_parser.CDS_parser import CDSParser

def parseFeatures(path, id, organism, record):

    logging.debug("#### id = " + str(id) + " ####")

    # Initialize variables
    DNA = None
    DNA_length = -1

    # Read DNA sequence from record
    try:
        DNA = record.seq
        DNA_length = len(DNA)
    except:
        logging.error("Unable to read DNA sequence from id = " + str(id))

    # Parse feature according to its type
    try:
        for feature in record.features:
            # if feature.type in region_type
            if feature.type == "CDS":
                logging.info("Found CDS in id = " + str(id))
                CDSParser(path, id, organism, DNA, DNA_length, feature)
            elif feature.type == "centromere":
                logging.info("Found centromere in id = " + str(id))
            elif feature.type == "intron":
                logging.info("Found intron in id = " + str(id))
            elif feature.type == "mobile_element":
                logging.info("Found mobile_element in id = " + str(id))
            elif feature.type == "telomere":
                logging.info("Found telomere in id = " + str(id))
            elif feature.type == "3'UTR":
                logging.info("Found 3'UTR in id = " + str(id))
            elif feature.type == "5'UTR":
                logging.info("Found 5'UTR in id = " + str(id))
            elif feature.type == "tRNA":
                logging.info("Found tRNA in id = " + str(id))
            elif feature.type == "rRNA":
                logging.info("Found rRNA in id = " + str(id))
            elif feature.type == "ncRNA":
                logging.info("Found ncRNA in id = " + str(id))
    except:
        logging.error("Unable to read feature(s) from id = " + str(id))