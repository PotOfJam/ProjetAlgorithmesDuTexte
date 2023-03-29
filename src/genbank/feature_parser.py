import logging

from DNA_parser.CDS_parser import CDSParser

def parseFeatures(id, record):

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
            if feature.type == "CDS":
                CDSParser(DNA, DNA_length, feature)
            elif feature.type == "centromere":
                logging.warning("Found centromere in id = " + str(id))
            elif feature.type == "intron":
                logging.warning("Found intron in id = " + str(id))
            elif feature.type == "mobile_element":
                logging.warning("Found mobile_element in id = " + str(id))
            elif feature.type == "telomere":
                logging.warning("Found telomere in id = " + str(id))
            elif feature.type == "3'UTR":
                logging.warning("Found 3'UTR in id = " + str(id))
            elif feature.type == "5'UTR":
                logging.warning("Found 5'UTR in id = " + str(id))
            elif feature.type == "tRNA":
                logging.warning("Found tRNA in id = " + str(id))
            elif feature.type == "rRNA":
                logging.warning("Found rRNA in id = " + str(id))
            elif feature.type == "ncRNA":
                logging.warning("Found ncRNA in id = " + str(id))
    except:
        logging.error("Unable to read feature(s) from id = " + str(id))