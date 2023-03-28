from Bio import Entrez
from Bio import SeqIO

Entrez.email = "fabien.allemand@etu.unistra.fr"

search_db = "nucleotide"
fetch_db = "nuccore"

region = ["CDS", "centrometre", "intron", "mobile_element", "telomere", "3'UTR", "5'UTR"] # tRNA, rRNA, ncRNA ???
valid_codon_start = ["ATG", "CTG", "TTG", "GTG", "ATA", "ATC", "ATT", "TTA"]
valid_codon_stop = ["TAA", "TAG", "TGA"]

id = "NC_000021" # For testing purpose, very small organism

print("#### id =", id, "####")
handle = Entrez.efetch(db=fetch_db, id=id, rettype="gbwithparts", retmode="text")
try:
    record = SeqIO.read(handle, "genbank")
except:
    print("ERROR: Unable to read from id =", id)
handle.close()

DNA_seq = None
DNA_seq_len = -1
try:
    DNA_seq = record.seq
    DNA_seq_len = len(DNA_seq)
    print(DNA_seq)
except:
    print("ERROR: Unable to read sequence(s) from id =", id)

try:
    for feature in record.features:
        if feature.type in region:
            print("############")
            print("TYPE =", feature.type)
            print("START =", feature.location.start)
            print("STOP =", feature.location.end)
            print("STRAND =", feature.location.strand)
            print("QUALIFIER =", feature.qualifiers)
            print("LOCATION_OPERATOR =", feature.location_operator)
            if str(feature.location_operator) == "join": # str??
                print("PARTS = ", feature.location.parts)
            if ("<" in str(feature.location.start)) or (">" in str(feature.location.start)) or ("<" in str(feature.location.end)) or (">" in str(feature.location.end)):
                print("ERROR: Invalid sequence start/stop (%s,%s)" % (feature.location.start, feature.location.end))
                continue
            sequence_start = int(feature.location.start)
            sequence_end = int(feature.location.end)
            print("DEBUG =", sequence_start, sequence_end)
            if sequence_start < 0:
                print("ERROR: Invalid sequence start/stop (%d,%d)(sequence_start < 0)" % (sequence_start, sequence_end))
                continue
            if sequence_end > DNA_seq_len:
                print("ERROR: Invalid sequence start/stop (%d,%d)(sequence_end > DNA_seq_len)" % (sequence_start, sequence_end))
                continue
            if sequence_end <= sequence_start:
                print("ERROR: Invalid sequence start/stop (%d,%d)(sequence_end <= sequence_start)" % (sequence_start, sequence_end))
                continue
            feature_DNA_seq = DNA_seq[sequence_start : sequence_end] 
            if len(feature_DNA_seq) % 3 != 0:
                print("ERROR: Invalid sequence length (%d)" % len(feature_DNA_seq))
                continue
            # print(type(feature_DNA_seq))
            if feature.strand == -1:
                feature_DNA_seq = feature_DNA_seq.reverse_complement()
            print("SEQUENCE =", str(feature_DNA_seq))
            codon_start = feature_DNA_seq[0:3]
            codon_stop = feature_DNA_seq[len(feature_DNA_seq)-3:]
            if codon_start not in valid_codon_start:
                print("ERROR: Invalid codon start (%s)" % codon_start)
                continue
            if codon_stop not in valid_codon_stop:
                print("ERROR: Invalid codon stop (%s)" % codon_stop)
                continue
            if not all(base in "ATGC" for base in feature_DNA_seq):
                print("ERROR: Invalid base found in sequence")
                continue
except:
    print("ERROR: Unable to read feature(s) from id =", id)
