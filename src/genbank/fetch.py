import logging, traceback
from Bio import Entrez, SeqIO
import random

def fetchFromID(id, fetch_db="nuccore", rettype="gbwithparts"):

    logging.info("Fetching data from GenBank database...(please wait)")

    # Set-up for request
    Entrez.email = "fabien.allemand@etu.unistra.fr"

    handle = None

    fetched = False
    delay = 0.5 # In seconds
    while not fetched:
        # Fetch data from database
        
        try:
            handle = Entrez.efetch(db=fetch_db, id=id, rettype=rettype, retmode="text")
        except IOError:
            delay += random.uniform(0, 0.5)
            logging.warning("Unable to fetch id = " + str(id) + " from fetch_db = " + fetch_db + ", retrying in " + str(delay) + " seconds")
        except:
            print(traceback.format_exc())
            logging.error("Unable to fetch id = " + str(id) + " from fetch_db = " + fetch_db)
            if(handle!=None):
                # Close record
                try:
                    handle.close()
                except:
                    logging.error("Unable to close handle")
                    return
            return
        else:
            fetched = True

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