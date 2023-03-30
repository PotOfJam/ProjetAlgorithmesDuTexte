import logging
from Bio import Entrez

def searchID(organism, search_db = "nucleotide"):

    if organism == "":
        logging.warning("Unable to retrieve ID(s) from search_db = " + search_db + " (no organism provided)")
        return []
    
    handle = Entrez.esearch(db=search_db, term="(" + organism + "[Organism] AND NC_*[Accession])", retmax ="999999999", usehistory='y')
    record = Entrez.read(handle)
    handle.close()

    return record["IdList"]