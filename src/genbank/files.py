from Bio import Entrez
from Bio import SeqIO

# Entrez.email = "fabien.allemand@etu.unistra.fr"

# orga = "Homo sapiens"
# cible = "gene"
# db = "nuccore"
# db = "nucleotide"
# handle = Entrez.esearch(db=db, term="(" + orga + "[Organism] AND NC_*[Accession])", retmax ="99999999", usehistory='y')
# record = Entrez.read(handle)
# handle.close()
# ids = record["IdList"]
# print(len(ids))
# print(ids)
# for id in ids:
#     handle = Entrez.efetch(db=db, id=id, rettype="gb", retmode="text")
#     record = SeqIO.read(handle, "genbank")
#     handle.close()
#     print(record)
#     Dna_seq = record.seq
#     print(Dna_seq)
#     cible_positions = []
#     for feature in record.features:
#         print(feature.type)
#         if feature.type == cible:
#             cible_positions.append((feature.location.start, feature.location.end))
#             print(feature.location.start, feature.location.end)

from enum import Enum
class Color(Enum):
    RED = 1
    GREEN = 2     
    BLUE = 3

"""
Color.RED
RED
1
Color.RED
"""
print(Color.RED)
print(Color.RED.name)
print(Color.RED.value)
print(Color["RED"])

val = "RED"
print(f"Value is {Color[val].value}")

import os
print(os.listdir("/home/fabien/TPS/test"))