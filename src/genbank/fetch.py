import logging, traceback
from Bio import Entrez, SeqIO

def fetchFromID(id, fetch_db="nuccore", rettype="gbwithparts"):

    logging.info("Fetching data from GenBank database...(please wait)")

    # Set-up for request
    Entrez.email = "fabien.allemand@etu.unistra.fr"

    # Fetch data from database
    handle = None
    try:
        handle = Entrez.efetch(db=fetch_db, id=id, rettype=rettype, retmode="text")
    except:
        print(traceback.format_exc())
        logging.error("Unable to fetch id = " + str(id) + " from fetch_db = " + fetch_db)
        return

    # Read data
    record = None
    try:
        record = SeqIO.read(handle, "genbank")
    except:
        logging.error("Unable to read id = " + str(id))
        return

    # Close record
    try:
        handle.close()
    except:
        logging.error("Unable to close handle")
        return

    return record