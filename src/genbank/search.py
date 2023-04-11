import logging
from Bio import Entrez

def searchID(organism, search_db = "nucleotide"):

    logging.info("Looking for organism ids...")

    # Set-up for request
    Entrez.email = "fabien.allemand@etu.unistra.fr"

    if organism == "":
        logging.warning("Unable to retrieve ID(s) from search_db = " + search_db + " (no organism provided)")
        return []
    
    handle = Entrez.esearch(db=search_db, term="(" + organism + "[Organism] AND NC_*[Accession])", retmax ="999999999", usehistory='y')
    record = Entrez.read(handle)
    handle.close()

    logging.info("Found %d ids to analyse:" % len(record["IdList"]))
    for id in record["IdList"]:
        logging.info("-> %s" % id)

    return record["IdList"]