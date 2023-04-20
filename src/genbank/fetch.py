import logging, traceback
from Bio import Entrez, SeqIO
import time, random

from ..app.logger import emitLog, Log

def fetchFromID(id, fetch_db="nuccore", rettype="gbwithparts", worker=None):

    emitLog(Log.INFO, "Fetching data from GenBank database...(please wait)", worker)

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
            emitLog(Log.ERROR, "Unable to fetch id = " + str(id) + " from fetch_db = " + fetch_db + ", retrying in " + str(delay) + " seconds", worker)
            time.sleep(delay)
        except:
            emitLog(Log.ERROR, traceback.format_exc(), worker)
            emitLog(Log.ERROR, "Unable to fetch id = " + str(id) + " from fetch_db = " + fetch_db, worker)
            if(handle!=None):
                # Close record
                try:
                    handle.close()
                except:
                    emitLog(Log.ERROR, "Unable to close handle", worker)
                    return
            return
        else:
            fetched = True

    # Read data
    record = None
    try:
        record = SeqIO.read(handle, "genbank")
    except:
        emitLog(Log.ERROR, "Unable to read id = " + str(id), worker)
        return

    # Close record
    try:
        handle.close()
    except:
        emitLog(Log.ERROR, "Unable to close handle", worker)
        return

    return record