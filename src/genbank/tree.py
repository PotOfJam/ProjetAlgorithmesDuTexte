import os, sys, logging
import csv
import requests
import pytz
from dateutil.parser import parse as parsedate
import datetime

from genbank.search import searchID
from genbank.fetch import fetchFromID

def updateTree(overview_file="../overview.txt"):
    """
    Create or update the tree in the results folder located in the main folder.

    Args:
        overview_file (str, optional): Path to the overview.txt file. Defaults to "../overview.txt".
    """
    # If overview file exists, check for update
    if os.path.exists(overview_file):
        utc = pytz.UTC

        # URL
        url = 'https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/overview.txt'
        r = requests.head(url)
        url_time = r.headers['last-modified']
        url_date = parsedate(url_time)
        # url_date = utc.localize(url_date)

        # File
        file_date = datetime.datetime.fromtimestamp(os.path.getmtime(overview_file))
        file_date = utc.localize(file_date)

        # Check for update
        # print(url_date)
        # print(file_date)
        if url_date > file_date: # NEED TO BE TESTED
            downloadAndUpdateTree(overview_file)
    # Download overview file and create tree
    else:
        downloadAndUpdateTree(overview_file)


def downloadAndUpdateTree(overview_file):
    """
    Download latest version of the overview.txt file and update the tree in the Results folder located in the main folder.

    Args:
        overview_file (str): Path to the overview.txt file.
    """
    # Donwload tree decription file from genbank
    url = 'https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/overview.txt'
    r = requests.get(url, allow_redirects=True)
    try:
        open(overview_file, 'wb').write(r.content)
    except:
        logging.error("No such file: " + overview_file)
        sys.exit(1)

    # Parse tree description file and create tree
    with open(overview_file, newline='') as csvfile:
        file_tree = csv.reader(csvfile, delimiter='\t')
        next(file_tree, None)
        for row in file_tree:
            organism, kingdom, group, subgroup, *_ = row
            # organism = organism.replace(" ", "_").replace("?", "").replace("/", "").replace("[", "").replace("]", "")
            organism = organism.replace("?", "").replace("/", "").replace("[", "").replace("]", "")
            os.makedirs(os.path.join("../Results/Organisme", kingdom, group, subgroup, organism), exist_ok=True)  # Ignore existing folders


def findOrganisms(folder):
    """
    Find organisms name located in a folder.

    Args:
        folder (string): Path of the folder to explore.

    Returns:
        List: List of organisms name and path.
    """
    # Look for sub folders
    sub_folders = []
    try:
        sub_folders = [sub_folder for sub_folder in os.listdir(folder) if os.path.isdir(os.path.join(folder, sub_folder))]
    except Exception as e:
        logging.error(e)
        sys.exit(1)

    # If no sub folder, folder is an organism
    if sub_folders == []:
        return [os.path.basename(os.path.normpath(folder))]
    # Else, recursively look for organisms
    else:
        organisms = []
        for sub_folder in sub_folders:
            organisms.append(findOrganisms(os.path.join(folder, sub_folder)))
        return [organism for sublist in organisms for organism in sublist] # Flatten result


def findSubFolders(path):

    # Look for sub folders
    sub_folders = []
    try:
        sub_folders = [sub_folder for sub_folder in os.listdir(path) if os.path.isdir(os.path.join(path, sub_folder))]
    except Exception as e:
        logging.error(e)
        return []
    
    for sub_folder in sub_folders:
        sub_folder = os.path.join(path, sub_folder)
    return sub_folders

def findOrganismsAndPath(selected_folder_path):

    # Already done
    # if not os.path.isdir(selected_folder_path):
    #     logging.error("...")

    organisms = []
    sub_folders = findSubFolders(selected_folder_path)

    while sub_folders != []:
        for sub_folder in sub_folders:
            sub_folders += findSubFolders(sub_folder)
    


def convertRecordDate(modification_date):
    day, month, year = modification_date.split("-")

    day = int(day)
    months = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]
    month = months.index(month) + 1
    year = int(year)

    modification_date = datetime.datetime(year, month, day)
    
    return modification_date


def findLastUpdateDate(ids):

    last_modification_date = -1
    for id in ids:
        record = fetchFromID(id, rettype="genbank") # Create custom function
        modification_date = record.annotations["date"]
        modification_date = convertRecordDate(modification_date)

        if last_modification_date == -1:
            last_modification_date = modification_date
        elif last_modification_date < modification_date:
            last_modification_date = modification_date   

    return last_modification_date


def findLastParsingDate(path):

    # Find files in directory
    files = []
    try:
        files = [file for file in os.listdir(path) if os.path.isdir(os.path.join(path, file))]
        print("$ files =", files)
    except Exception as e:
        logging.error(e)
        sys.exit(1)

    # Find date
    parsing_date = -1
    for file in files:
        file_date = datetime.datetime.fromtimestamp(os.path.getmtime(file))
        print("$", file, file_date)
        if parsing_date < file_date:
            parsing_date = file_date

    return parsing_date


def genbankDate(path, organism):

    # MULTITHREADING!!

    # Genbank date
    ids = searchID(organism)
    last_modification_date = findLastUpdateDate(ids)
    print("# last_modification_date =", last_modification_date)

    # Parsing date
    last_parsing_date = findLastParsingDate(path)
    print("# last_parsing_date =", last_parsing_date)