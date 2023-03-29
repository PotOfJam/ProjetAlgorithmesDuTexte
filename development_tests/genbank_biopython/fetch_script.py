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
            print("STRAND =", feature.location.strand)
            print("QUALIFIER =", feature.qualifiers)
            print("START =", feature.location.start)
            print("STOP =", feature.location.end)
            print("LOCATION_OPERATOR =", feature.location_operator)

            # Load parts positions
            sequence_pos = []
            error = False
            if str(feature.location_operator) == "join": # str??
                # print("PARTS =", feature.location.parts)
                for part_pos in feature.location.parts:
                    start_pos = 0
                    end_pos = 0
                    print("PART_POS =", part_pos)
                    print(type(part_pos))
                    start_pos = part_pos.start
                    end_pos = part_pos.end
                    print("DEBUG =", start_pos, end_pos)
                    if ("<" in str(start_pos)) or (">" in str(start_pos)) or ("<" in str(end_pos)) or (">" in str(end_pos)):
                        print("WARNING: Invalid sequence start/stop (%s,%s) [DISCARDED]" % (start_pos, end_pos))
                        error = True
                        break
                    sequence_pos.append((int(start_pos), int(end_pos)))
                if error:
                    continue
            else:
                if ("<" in str(feature.location.start)) or (">" in str(feature.location.start)) or ("<" in str(feature.location.end)) or (">" in str(feature.location.end)):
                    print("ERROR: Invalid sequence start/stop (%s,%s) [DISCARDED]" % (feature.location.start, feature.location.end))
                    continue
                sequence_pos.append((int(feature.location.start), int(feature.location.end)))

            # Test parts positions
            for position in sequence_pos:
                error = False
                sequence_start = position[0]
                sequence_end = position[1]
                if sequence_start < 0:
                    print("ERROR: Invalid sequence start/stop (%d,%d)(sequence_start < 0)" % (sequence_start, sequence_end))
                    error = True
                if sequence_end > DNA_seq_len:
                    print("ERROR: Invalid sequence start/stop (%d,%d)(sequence_end > DNA_seq_len)" % (sequence_start, sequence_end))
                    error = True
                if sequence_end <= sequence_start:
                    print("ERROR: Invalid sequence start/stop (%d,%d)(sequence_end <= sequence_start)" % (sequence_start, sequence_end))
                    error = True

            def getFirst(tuple):
                return tuple[0]
            
            sequence_pos.sort(key=getFirst)
            print("SEQUENCE_POS =", sequence_pos)

            # Recreate sequence
            feature_DNA_seq = ""
            for sub_sequence in sequence_pos: # Réordonner pour être sûr
                sub_sequence_start = sub_sequence[0]
                sub_sequence_end = sub_sequence[1]
                print("DNA_sub_sequence =", DNA_seq[sub_sequence_start : sub_sequence_end])
                feature_DNA_seq += DNA_seq[sub_sequence_start : sub_sequence_end]


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
