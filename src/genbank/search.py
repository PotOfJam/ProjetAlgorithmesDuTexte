import logging
from Bio import Entrez
from ..app.logger import emitLog, Log

def searchID(organism, search_db = "nucleotide", worker=None):

    emitLog(Log.INFO, "Looking for organism ids...", worker)

    # Set-up for request
    Entrez.email = "fabien.allemand@etu.unistra.fr"

    if organism == "":
        emitLog(Log.ERROR, "Unable to retrieve ID(s) from search_db = " + search_db + " (no organism provided)", worker)
        return []
    
    handle = Entrez.esearch(db=search_db, term="(" + organism + "[Organism] AND NC_000001:NC_999999[ACCN])", retmax ="999999999", usehistory='y', idtype="acc")
    record = Entrez.read(handle)
    handle.close()

    emitLog(Log.INFO, "Found %d ids to analyse:" % len(record["IdList"]), worker)
    for id in record["IdList"]:
        emitLog(Log.INFO, "-> %s" % id, worker)

    return record["IdList"]