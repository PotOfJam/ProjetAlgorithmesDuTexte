from Bio import Entrez
import requests
import csv
import os

url = 'https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/overview.txt'
r = requests.get(url, allow_redirects=True)
open('overview.txt', 'wb').write(r.content)

with open('overview.txt', newline='') as csvfile:
    file_tree = csv.reader(csvfile, delimiter='\t')
    next(file_tree,None)
    for row in file_tree:
        organism, kingdom , group, subgroup, *_ = row
        organism = organism.replace(" ","_").replace("?","").replace("/","").replace("[","").replace("]","")
        os.makedirs(os.path.join("Organisme",kingdom,group,subgroup,organism),exist_ok=True)
