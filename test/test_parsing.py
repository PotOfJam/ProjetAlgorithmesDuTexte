import os, sys
import pytest
from src import *
from filecmp import cmp

# get absolute path of test_parsing.py
script_dir = os.path.dirname(__file__)

module_dir = os.path.join(script_dir, "..", "src", "genbank")
sys.path.append(module_dir)

path_reference_file = os.path.join(os.path.dirname(__file__), "test_files", "reference_files")

#==== Search =================================================================#

def searchAndCompare(organism, reference_ids):
    """*
    For a species, get its GenBank file IDs with searchID and compare them to the reference IDs.

    Args:
        organism (string): the species for which we want to get the reference ids
        reference_ids (string): the supposed to be obtained reference ids

    Returns:
        bool: returns true if the supposed to be obtained reference ids are the same as the obtained reference ids, otherwise returns false
    """
    ids = search.searchID(organism)
    print("Got: ", ids)
    return ids == reference_ids


def test_searchID():
    """
    Test of the function searchID.
    Note: parameters contains tuples of organism name and list of expected GenBank files IDs.
    """
    parameters = [("Homo sapiens", ["NC_013993.1", "NC_012920.1", "NC_011137.1", "NC_060925.1", "NC_060926.1", "NC_060927.1", "NC_060928.1", "NC_060929.1", "NC_060930.1", "NC_060931.1", "NC_060932.1", "NC_060933.1", "NC_060934.1", "NC_060935.1", "NC_060936.1", "NC_060937.1", "NC_060938.1", "NC_060939.1", "NC_060940.1", "NC_060941.1", "NC_060942.1", "NC_060943.1", "NC_060944.1", "NC_060945.1", "NC_060946.1", "NC_060947.1", "NC_060948.1", "NC_000001.11", "NC_000002.12", "NC_000003.12", "NC_000004.12", "NC_000005.10", "NC_000006.12", "NC_000007.14", "NC_000008.11", "NC_000009.12", "NC_000010.11", "NC_000011.10", "NC_000012.12", "NC_000013.11", "NC_000014.9", "NC_000015.10", "NC_000016.10", "NC_000017.11", "NC_000018.10", "NC_000019.10", "NC_000020.11", "NC_000021.9", "NC_000022.11", "NC_000023.11", "NC_000024.10"]), ("Candidatus Carsonella ruddii", ["NC_018418.1", "NC_018417.1", "NC_008512.1", "NC_021894.1", "NC_018416.1", "NC_018415.1", "NC_018414.1"])]
    for organism, reference_ids in parameters:
        print("Organism: ", organism)
        print("Expected: ", reference_ids)
        if not searchAndCompare(organism, reference_ids):
            assert False

#==== Fetch ==================================================================#

def test_fetchFromID():
    """
    Test for a given reference id, if we get a record when fetching. 
    """
    reference_ids = ["NC_018416"]
    for id in reference_ids:
        record = fetch.fetchFromID(id)

#==== Parse ==================================================================#

def parseAndCompare(region_type, path, id, organism, record, path_reference_file, file_name):
    """
    For a species, one of its id and some region types, parse the associated record.
    The record is stored in a path and compared with what the record should look like,
    which is a file in a given path.

    Args:
        region_type (list of strings): list of region to extract.
        path (string): path in which we will temporarly store the parsed record
        id (string): id of the species
        organism (string): species for which we want to get the record
        record (SeqRecord): the record we want to parse
        path_reference_file (string): the path in which the file of reference is located
        file_name (string): name of the file of reference

    Returns:
        bool: true if the parsed file is the same as the reference file, false otherwise
    """
    if not os.path.exists(path):
        os.makedirs(path)
    
    parseFeatures(region_type, path, id, organism, record)
    path_file_name = os.path.join(path, file_name)
    path_ref = os.path.join(path_reference_file, file_name)
    return cmp(path_file_name, path_ref)


def test_parseAndCompare():
    """
    Test of the function parseFeatures.
    Note: parameters contains tuples containing a list of region types, two paths a reference id, a species name and a record to parse
    """
    test_path = os.path.join(os.path.dirname(__file__), "test_files", "parsed_files")
    parameters = [
        #(["CDS"], test_path, "NC_000021.9", "Homo sapiens", fetch.fetchFromID("NC_000021.9")),
                  (["CDS"], test_path, "NC_018418.1", "Candidatus Carsonella ruddii", fetch.fetchFromID("NC_018418.1"), "CDS_Candidatus_Carsonella_ruddii_NC_018418.txt"),
                  (["CDS"], test_path, "NC_018416", "Candidatus Carsonella ruddii", fetch.fetchFromID("NC_018416"), "CDS_Candidatus_Carsonella_ruddii_NC_018416.txt")]

    for region_type, path, id, organism, record, file_name in parameters:
        if parseAndCompare(region_type, path, id, organism, record, path_reference_file, file_name) == False:
            assert False
        # else, file will be concatened with itself at next textlaunch
        os.remove(os.path.join(path, file_name))

    assert True
