import logging
from Bio import Entrez
from ..app.parser_thread import emitLog

def searchID(organism, search_db = "nucleotide", worker=None):

    emitLog(worker, "Looking for organism ids...")

    # Set-up for request
    Entrez.email = "fabien.allemand@etu.unistra.fr"

    if organism == "":
        emitLog(worker, "Unable to retrieve ID(s) from search_db = " + search_db + " (no organism provided)")
        return []
    
    handle = Entrez.esearch(db=search_db, term="(" + organism + "[Organism] AND NC_000001:NC_999999[ACCN])", retmax ="999999999", usehistory='y', idtype="acc")
    record = Entrez.read(handle)
    handle.close()

    emitLog(worker, "Found %d ids to analyse:" % len(record["IdList"]))
    for id in record["IdList"]:
        emitLog(worker, "-> %s" % id)

    return record["IdList"]