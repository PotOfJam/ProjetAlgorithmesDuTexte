from parse import lectureBorne
from parse import XinY
from parse import geneSequences

import os

#os.chdir(os.path.dirname(os.path.abspath("test.txt")))
chaine="join (1..2,4..5,  abc   23..56) \\"

bornes=lectureBorne(chaine)
print(bornes)

chaine="atg ccccatg taa   at g  atgcccc tg a atg cccc"
alphabet="tacg"
x=XinY(["atg","taa","tga","tag"],chaine,alphabet)
print(x)
res=geneSequences(["atg","taa","tga","tag"],chaine,alphabet)
print(res)
for bornes in res:
    print(chaine[bornes[0]:bornes[1]])

