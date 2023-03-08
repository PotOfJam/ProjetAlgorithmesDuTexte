import os
import csv
import requests
import pytz
from dateutil.parser import parse as parsedate
import datetime


def updateTree(overview_file):
    """
    Update tree in the Results folder located in the main folder
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
        file_date = datetime.datetime.fromtimestamp(os.path.getmtime("overview.txt"))
        file_date = utc.localize(file_date)

        # Check for update
        print(url_date)
        print(file_date)
        if url_date > file_date:
            downloadAndUpdateTree()
    # Download overview file and create tree
    else:
        downloadAndUpdateTree()


def downloadAndUpdateTree():
    # Donwload tree decription file from genbank
    url = 'https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/overview.txt'
    r = requests.get(url, allow_redirects=True)
    open('overview.txt', 'wb').write(r.content)

    # Parse tree description file and create tree
    with open('overview.txt', newline='') as csvfile:
        file_tree = csv.reader(csvfile, delimiter='\t')
        next(file_tree, None)
        for row in file_tree:
            organism, kingdom, group, subgroup, *_ = row
            organism = organism.replace(" ", "_").replace(
                "?", "").replace("/", "").replace("[", "").replace("]", "")
            os.makedirs(os.path.join("Organisme", kingdom, group, subgroup,
                        organism), exist_ok=True)  # Ignore existing folders
