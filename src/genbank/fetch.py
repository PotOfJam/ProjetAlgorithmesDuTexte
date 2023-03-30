import logging
from Bio import Entrez, SeqIO

def fetchFromID(id, fetch_db="nuccore"):

    # Set-up for request
    Entrez.email = "fabien.allemand@etu.unistra.fr"

    # Fetch data from database
    try:
        handle = Entrez.efetch(db=fetch_db, id=id, rettype="gbwithparts", retmode="text")
    except:
        logging.error("Unable to fetch id = " + str(id) + " from fetch_db = " + fetch_db)

    # Read data
    try:
        record = SeqIO.read(handle, "genbank")
    except:
        logging.error("Unable to read id = " + str(id))

    # Close record
    handle.close()

    return record