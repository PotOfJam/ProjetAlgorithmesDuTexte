# fonction reprise du cour
def dernier_occurence(x,m):

    # ici alphabet = acgt
    dern_occ=[m,m,m,m]

    for k in range(0,m):
        if x[k]=="a":
            dern_occ[0]=m-k-1
        elif x[k]=="c":
            dern_occ[1]=m-k-1
        elif x[k]=="g":
            dern_occ[2]=m-k-1
        elif x[k]=="t":
            dern_occ[3]=m-k-1

    print("dern_occ",dern_occ)

    return dern_occ

# lis 3 lettre du fichier, ignore les lettres qui ne sont pas dans l'alphabet
def dansAlphabet(pointeur,fichier):

    i=0   #nombre de lettre lu apparetenant à l'aphabet
    tot=0 #nombre de lettre lu au total

    chaine=""
    boucle=1

    while i<3:
        c=fichier.read(1)
        if (not c):
            boucle=0
            break
        elif not( c=="a" or c=="c" or c=="g" or c=="t"):
            tot=tot+1
            fichier.seek(pointeur+tot,0)
        else:
            chaine=chaine+c
            tot=tot+1
            i=i+1

    return [boucle,tot, chaine]

# renvoit les listes des occurences de atg/taa/tag/tga dans le fichier
def XinY(x,fichier):
    occ=[[],[],[],[]]

    pointeur=0
    fichier.seek(pointeur,0)

    boucle=1
    while(boucle):

        [boucle,tot,chaine]=dansAlphabet(pointeur,fichier)
        trouve=0

        if(not chaine):
            break

        elif(x[0] == chaine):
            trouve=1
            occ[0].append(pointeur)

        elif(x[1] == chaine):
            trouve=1
            occ[1].append(pointeur)

        elif(x[2] == chaine):
            trouve=1
            occ[2].append(pointeur)

        elif(x[3] == chaine):
            trouve=1
            occ[3].append(pointeur)

        if(trouve):
            pointeur=pointeur+tot
        else:
            pointeur=pointeur+1

        fichier.seek(pointeur,0)
    return occ

# affiche les mot du fichier 
def verifie(x,fichier):
    for i in x:
        print(i)
        for j in i:
            pointeur=j
            fichier.seek(pointeur,0)
            [boucle,tot,chaine]=dansAlphabet(pointeur,fichier)
            print(" ",chaine)