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


# def searchID(organism, search_db = "nucleotide"):
#     # Set your NCBI email and API key
#     Entrez.email = "fabien.allemand@etu.unistra.fr"

#     if organism == "":
#         logging.warning("Unable to retrieve ID(s) from search_db = " + search_db + " (no organism provided)")
#         return []

#     # Define your search query
#     search_term = "(" + organism + "[Organism] AND NC_*[Accession])"

#     # Fetch the first 600 records
#     handle = Entrez.esearch(db="nuccore", term=search_term, retmax=600)
#     record = Entrez.read(handle)
#     handle.close()

#     # Extract the list of IDs from the first set of results
#     id_list = record["IdList"]

#     # Fetch additional sets of results iteratively, if available
#     while "WebEnv" in record and "QueryKey" in record:
#         # Fetch the next set of results using the WebEnv and QueryKey
#         handle = Entrez.esearch(db="nuccore", term=search_term, retmax=600, webenv=record["WebEnv"], query_key=record["QueryKey"])
#         record = Entrez.read(handle)
#         handle.close()
#         # Extend the list of IDs with the new set of results
#         id_list.extend(record["IdList"])

#     # Print the list of IDs
#     logging.info("Found %d ids to analyse:" % len(record["IdList"]))
#     for id in record["IdList"]:
#         logging.info("-> %s" % id)

#     return record["IdList"]