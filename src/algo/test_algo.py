from XinY import XinY,verifie
import os

#os.chdir(os.path.dirname(os.path.abspath("test.txt")))

f = open("test/test.txt", "r")

chaine=f.read(3)

print("chaine:", chaine)

res=XinY(["atg","taa","tag","tga"],f)

print("res",res)

verifie(res,f)


f.close()

with open("test.txt", "r") as fichier:
	print (fichier.read())
	


#j-m+1 - j
#i - i+m-1
#j = j+ der
#i+m-1 = i+m-1+der
#i = i+der
