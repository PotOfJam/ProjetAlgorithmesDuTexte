from Bio import Entrez

Entrez.email = "fabien.allemand@eyu;unistra.fr"
handle = Entrez.egquery(term="Opuntia AND rpl16")
record = Entrez.read(handle)
for row in record["eGQueryResult"]:
    if row["DbName"] == "nuccore":
        print(row["Count"])

handle = Entrez.esearch(db="nuccore", term="Opuntia AND rpl16")
record = Entrez.read(handle)
gi_list = record["IdList"]
print(gi_list)

handle = Entrez.efetch(db="nuccore", id=gi_list, rettype="gb", retmode="text")

text = handle.read()
print(text)

print("END")