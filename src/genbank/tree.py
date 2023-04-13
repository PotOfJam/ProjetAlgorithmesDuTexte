import os, sys, logging, threading
import csv
import requests
import pytz
from dateutil.parser import parse as parsedate
import datetime

import genbank.fetch

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
            logging.info("Updating Results file tree...")
            downloadAndUpdateTree(overview_file)
    # Download overview file and create tree
    else:
        logging.info("Creating Results file tree...")
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


def findSubFolders(path):
    """
    Find sub-folders.

    Args:
        path (string): Path of the folder to look inside.

    Returns:
        list: Paths of the sub-folders.
    """
    # Look for sub-folders
    sub_folders = []
    try:
        sub_folders = [sub_folder for sub_folder in os.listdir(path) if os.path.isdir(os.path.join(path, sub_folder))]
    except Exception as e:
        logging.error(e)
        return []
    
    for sub_folder in sub_folders:
        sub_folder = os.path.join(path, sub_folder)
    return sub_folders


def findLastSubFolders(selected_folder_path):
    """
    Find folders at the bottom layer of the hierachy starting at a given path.

    Args:
        selected_folder_path (string): Path of the folder at the top of the hierarchy.


    Returns:
        list: Paths of the sub-folders located at the bottom of the hierarchy.
    """
    sub_folders = findSubFolders(selected_folder_path)
    
    # Selected folder does not contain any sub-folder
    if sub_folders == []:
        return [selected_folder_path]

    # Selected folder contains sub-folder(s)
    while sub_folders != []:
        for sub_folder in sub_folders:
            sub_sub_folders = findSubFolders(sub_folder)
            if sub_sub_folders != []:
                sub_folders += sub_sub_folders
                sub_folders.remove(sub_folder)

    return sub_folders


def findOrganisms(selected_folder_path):
    """
    Find organisms contained in a given folder.

    Args:
        selected_folder_path (string): Path of the folder to look inside.

    Returns:
        list: Tuples containing the name of the organism and the path to its folder.
    """
    logging.info("Looking for organism(s) in the selected folder...")

    organisms = []

    organisms_paths = findLastSubFolders(selected_folder_path)
    for path in organisms_paths:
        organisms.append((os.path.basename(os.path.normpath(path)), path))

    logging.info("Found %d organisms to analyse:" % len(organisms))
    for organism, _ in organisms:
        logging.info("-> %s" % organism)

    return organisms


def convertRecordDate(modification_date):
    """
    Convert record (GenBank file) modification date to datetime object.

    Args:
        modification_date (string): Last modification date of the GenBank file.

    Returns:
        datetime.datetime: Last modification date of the GenBank file.
    """
    day, month, year = modification_date.split("-")

    day = int(day)
    months = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]
    month = months.index(month) + 1
    year = int(year)

    modification_date = datetime.datetime(year, month, day)
    
    return modification_date


def findLastUpdateDate(ids):
    """
    Find the date at wich the organism (all GenBank files) has been last modified in the GenBank database.

    Args:
        ids (list): List of GenBank IDs related to an organism.

    Returns:
        datetime.datetime: Last modification date.
    """

    logging.info("Looking for last update on GenBank...")

    last_modification_date = None
    for id in ids:
        record = genbank.fetch.fetchFromID(id, rettype="genbank")
        modification_date = record.annotations["date"]
        modification_date = convertRecordDate(modification_date)

        if last_modification_date is None:
            last_modification_date = modification_date
        elif last_modification_date < modification_date:
            last_modification_date = modification_date   

    return last_modification_date


def findLastParsingDate(path):
    """
    Find the date at wich the organism (all parsing results files) has been last modified in the local "Result" folder.

    Args:
        path (string): Path of the organism's folder.

    Returns:
        datetime.datetime: Last modification date.
    """

    logging.info("Looking for last local update...")

    # Find files in directory
    files = []
    try:
        files = [file for file in os.listdir(path) if os.path.isfile(os.path.join(path, file))]
    except Exception as e:
        logging.error(e)
        return

    # Find date
    parsing_date = None
    for file in files:
        file_date = datetime.datetime.fromtimestamp(os.path.getmtime(os.path.join(path, file)))
        if parsing_date == None:
            parsing_date = file_date
        if parsing_date < file_date:
            parsing_date = file_date

    return parsing_date


def needParsing(organism_path, ids):
    """
    Check if an organism needs to be parsed.

    Args:
        organism_path (str): Path of the organism's folder.
        ids (list): List of GenBank IDs related to an organism.

    Returns:
        int: Number of GenBank files to parse.
    """
    # global last_update_date, last_parsing_date
    organism_files = [file for file in os.listdir(organism_path)]

    # Never parsed
    if organism_files == []:
        logging.info("Organism was never parsed, all files need to be parsed...")
        return len(ids)
    
    # Already parsed, check for update...    
    # Genbank date
    last_update_date = findLastUpdateDate(ids)
    logging.info("Last GenBank update: %s" % last_update_date)

    # Parsing date
    last_parsing_date = findLastParsingDate(organism_path)
    logging.info("Last local update: %s" % last_parsing_date)

    # Parsing needs to be updated
    if last_parsing_date is None or last_update_date is None:
        logging.warning("Invalid date")
        return 0
    if last_parsing_date < last_update_date:
        return len(ids)
    
    # Parsing is up to date!
    return 0
