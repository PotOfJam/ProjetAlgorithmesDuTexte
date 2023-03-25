import sys, os
import numpy as np
import csv
import requests
import pytz
from dateutil.parser import parse as parsedate
import datetime
from enum import Enum


class FolderType(Enum):
    """
    Type of folder.
    """
    RESULTS = 0
    ORGANISME = 1
    KINGDOM = 2
    GROUP = 3
    SUBGROUP = 4
    ORGANISM = 5


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
        print(url_date)
        print(file_date)
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
    open(overview_file, 'wb').write(r.content)

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
        folder (string_): Path of the folder to explore.

    Returns:
        List: List of organisms name.
    """
    # Look for sub folders
    sub_folders = []
    try:
        # sub_folders = list(filter(os.path.isdir, os.listdir(folder)))
        sub_folders = [sub_folder for sub_folder in os.listdir(folder) if os.path.isdir(os.path.join(folder, sub_folder))]
    except Exception as e:
        print(e)
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
    
# print(findOrganisms("../../Results/Organisme/Archaea/Candidatus Hydrothermarchaeota"))

