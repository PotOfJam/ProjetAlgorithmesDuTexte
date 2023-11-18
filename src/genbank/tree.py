import os
import csv
import requests
import pytz
from dateutil.parser import parse as parsedate
import datetime

from .fetch import *
from ..app.logger import emitLog, Log
import os

def updateTree(overview_file="overview.txt"):
    """
    Create or update the tree in the results folder located in the main folder.

    Args:
        overview_file (str, optional): Path to the overview.txt file. Defaults to "../overview.txt".
    """
    # If overview file exists, check for update
    if os.path.exists(overview_file):
        utc = pytz.UTC

        # URL
        url = "https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/overview.txt"
        r = requests.head(url)
        url_time = r.headers["last-modified"]
        url_date = parsedate(url_time)
        # url_date = utc.localize(url_date)

        # File
        file_date = datetime.datetime.fromtimestamp(os.path.getmtime(overview_file))
        file_date = utc.localize(file_date)

        # Check for update
        # print(url_date)
        # print(file_date)
        if url_date > file_date: # NEED TO BE TESTED
            emitLog(Log.INFO, "Updating Results file tree...")
            downloadAndUpdateTree(overview_file)
    # Download overview file and create tree
    else:
        emitLog(Log.INFO, "Creating Results file tree...")
        downloadAndUpdateTree(overview_file)


def downloadAndUpdateTree(overview_file, worker=None):
    """
    Download latest version of the overview.txt file and update the tree in the Results folder located in the main folder.

    Args:
        overview_file (str): Path to the overview.txt file.
    """
    # Donwload tree decription file from genbank
    url = "https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/overview.txt"
    try:
        r = requests.get(url, allow_redirects=True)
    except:
        emitLog(Log.ERROR, "Unable to download file: " + overview_file, worker)
        return

    try:
        open(overview_file, "wb").write(r.content)
    except:
        emitLog(Log.ERROR, "No such file: " + overview_file, worker)
        return

    # Parse tree description file and create tree
    with open(overview_file, newline="") as csvfile:
        file_tree = csv.reader(csvfile, delimiter="\t")
        next(file_tree, None)
        for row in file_tree:
            try:
                organism, kingdom, group, subgroup, *_ = row
                organism = organism.replace("?", "").replace("/", "").replace("[", "").replace("]", "").replace(":", "")
                os.makedirs(os.path.join("Results", "Organisme", kingdom, group, subgroup, organism), exist_ok=True)  # Ignore existing folders
            except:
                emitLog(Log.WARNING, "Invalid overview line: %s" % row, worker)


def findSubFolders(path, worker=None):
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
        sub_folders = [os.path.join(path, sub_folder) for sub_folder in os.listdir(path) if os.path.isdir(os.path.join(path, sub_folder))]
    except Exception as e:
        return []
    
    return sub_folders


def findLastSubFolders(selected_folder_path, worker=None):
    """
    Find folders at the bottom layer of the hierachy starting at a given path.

    Args:
        selected_folder_path (string): Path of the folder at the top of the hierarchy.

    Returns:
        list: Paths of the sub-folders located at the bottom of the hierarchy.
    """
    last_sub_folders = []
    sub_folders = findSubFolders(selected_folder_path, worker)
    
    # Selected folder does not contain any sub-folder
    if sub_folders == []:
        return [selected_folder_path]

    # Selected folder contains sub-folder(s)
    while sub_folders != []:
        for sub_folder in sub_folders:
            last_sub_folders += findLastSubFolders(sub_folder, worker)
            sub_folders.remove(sub_folder)

    return last_sub_folders


def findOrganisms(selected_folder_path, worker=None):
    """
    Find organisms contained in a given folder.

    Args:
        selected_folder_path (string): Path of the folder to look inside.

    Returns:
        list: Tuples containing the name of the organism and the path to its folder.
    """
    emitLog(Log.INFO, "Looking for organism(s) in the selected folder...", worker)

    organisms = []

    organisms_paths = findLastSubFolders(selected_folder_path, worker)
    for path in organisms_paths:
        organisms.append((os.path.basename(os.path.normpath(path)), path))

    emitLog(Log.INFO, "Found %d organisms to analyse:" % len(organisms), worker)
    for organism, _ in organisms:
        emitLog(Log.INFO, "-> %s" % organism, worker)

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


def findLastUpdateDate(ids, worker=None):
    """
    Find the date at which the organism (all GenBank files) has been last modified in the GenBank database.

    Args:
        ids (list): List of GenBank file IDs related to an organism.

    Returns:
        datetime.datetime: Last modification date.
    """

    emitLog(Log.INFO, "Looking for last update on GenBank...", worker)

    last_modification_date = None
    for id in ids:
        record = fetchFromID(id, rettype="genbank", worker=worker)
        modification_date = record.annotations["date"]
        modification_date = convertRecordDate(modification_date)

        if last_modification_date is None:
            last_modification_date = modification_date
        elif last_modification_date < modification_date:
            last_modification_date = modification_date   

    return last_modification_date


def findLastParsingDate(path, worker=None):
    """
    Find the date at which the organism (all parsing results files) has been last modified in the local "Result" folder.

    Args:
        path (string): Path of the organism's folder.

    Returns:
        datetime.datetime: Last modification date.
    """

    emitLog(Log.INFO, "Looking for last local update...", worker)

    # Find files in directory
    files = []
    try:
        files = [file for file in os.listdir(path) if os.path.isfile(os.path.join(path, file))]
    except Exception as e:
        emitLog(Log.ERROR, e, worker)
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


def needParsing(organism_path, ids, worker=None, region_type=['CDS']):
    """
    Check if an organism needs to be parsed.

    Args:
        organism_path (str): Path of the organism's folder.
        ids (list): List of GenBank file IDs related to an organism.
        worker ():

    Returns:
        int: Number of GenBank files to parse.
    """
    organism_files = [file for file in os.listdir(organism_path)]
    # Never parsed
    if organism_files == []:
        emitLog(Log.INFO, "Organism was never parsed, all files need to be parsed...", worker)
        return len(ids)
    # Partially parsed
    # Following lines do not work because there will be many files if user parse for CDS and intron and...
    # if len(organism_files) < len(ids):
    #     emitLog(Log.INFO, "Organism was partially parsed, all files need to be parsed...", worker)
    #     emitLog(Log.WARNING, "Local files are deleted", worker)
    #     for file in organism_files:
    #         if os.path.exists(file):
    #             os.remove(file)
    #     return len(ids)
    
    # Already parsed, check for update...    
    # Genbank date
    last_update_date = findLastUpdateDate(ids, worker)
    emitLog(Log.INFO, "Last GenBank update: %s" % last_update_date, worker)

    # Parsing date
    last_parsing_date = findLastParsingDate(organism_path, worker)
    emitLog(Log.INFO, "Last local update: %s" % last_parsing_date, worker)

    # Parsing needs to be updated
    if last_parsing_date is None or last_update_date is None:
        emitLog(Log.WARNING, "Invalid date", worker)
        return 0
    if last_parsing_date < last_update_date:
        emitLog(Log.INFO, "Parsing needs to be updated...", worker)
        return len(ids)
    else:
        for region in region_type:
            region_found = False
            for file in organism_files:
                if region in file:
                    region_found = True
                    break
            if not region_found:
                emitLog(Log.INFO, f"Region '{region}' needs to be parsed...", worker)
                return len(ids)
    # Parsing is up to date!
    return 0

def findRegionsToParse(organism_path, region_type, worker=None):
    organism_files = [file for file in os.listdir(organism_path)]

    regions_to_parse = []
    if organism_files == []:
        regions_to_parse = region_type
    else:
        for region in region_type:
                region_found = False
                for file in organism_files:
                    if region in file:
                        region_found = True
                if not region_found:
                    regions_to_parse.append(region)
    emitLog(Log.INFO,f"Regions that really need to be parse: {regions_to_parse}", worker)
    return regions_to_parse
