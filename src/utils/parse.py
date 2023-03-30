def lectureBorne(chaine):
    """
    lis les bornes au début du fichier des gènes
    Args:
        chaine (string): chaine contenant le fichier
    Returns:
        tableau de tableau de int: renvoit les différentes occurence des mots
    """
    bornes=[]

    numeros="1234567890"
    pointeur=0
    basseOUhaute=0

    uneBorne=[]
    b0=""
    b1=""

    while(1):

        char=chaine[pointeur]

        # on lit un numéros -> une borne est agrandie
        if(numeros.__contains__(char)):
            if(basseOUhaute):
                b1+=char
            else:
                b0+=char

        # on lit un point -> un autre point doit suivre
        elif(char=="." and chaine[pointeur+1]=="."):
            uneBorne.append(int(b0))

            basseOUhaute=1
            b0=""
            
            #on saute le "." prochain
            pointeur+=1 
        
        # on lit un caractère spéciale qui indique que l'on n'est plus dans les bornes
        elif(char=="\\" or char=="O"):
            break

        # on a lut un autre caractère -> on vient de lire l'entièreté de la borne
        elif(basseOUhaute ):
            uneBorne.append(int(b1))
            bornes.append(uneBorne)

            uneBorne=[]
            basseOUhaute=0
            b1=""


        pointeur+=1

    return bornes


def XinY(x,chaine,alphabet):
    """_summary_

    Args:
        x (tableau de string): mots à retrouvé dans le fichier
        chaine (string): chiane à parser
        alphabet (string): alphabet des lettres étudié, aide à determiner quels caractère sont à ignoré

    Returns:
        tableau de tableau de int: renvoit les différentes occurence des mots
    """
    tailleListe=len(x)
    print("taillListe: ",tailleListe, range(tailleListe))
    tailleMot=len(x[0])
    
    occ=[ [] for x in range(tailleListe)]

    pointeur=0
    chaine[0]

    boucle=1
    while(boucle<10 and pointeur < len(chaine)):

        [boucle,tot,mot,debut]=dansAlphabet(pointeur,chaine,alphabet,tailleMot)
        trouve=-1

        for k in range(tailleListe):
            if(x[k]==mot):
                trouve=k
                occ[k].append(debut)

        if(trouve>-1):
            pointeur=pointeur+tot
        else:
            pointeur=pointeur+1
        
        boucle+=1
        

    return occ

def dansAlphabet(pointeur,chaine,alphabet,tailleMot):
    """_summary_

    Args:
        pointeur (int): pointeur dans la chaine
        chaine (string): chaine à parser
        alphabet (string): liste de caractère à considérer pendant la lecture (les autres lettres sont ignoré)
        tailleMot (int): taille du mot à lire

    Returns:
        [boucle,tot, chaine,debut]: boucle dit si on est arrivé à la fin du fichier, tot et le nombre total de lettre lu, chaine est le mot lu, debut est un int renvoyant au début du mot
    """

    i=0   #nombre de lettre lu apparetenant à l'alphabet
    tot=0 #nombre de lettre lu au total
    debut=-1 # position du début du mot dans la chaine
    position=pointeur

    mot=""
    boucle=1

    while i<tailleMot and position < len(chaine):
        char=chaine[position]

        if(alphabet.__contains__(char)):
            mot+=char
            i+=1

            if(debut<0):
                debut=position
        
        
        tot=tot+1
        position+=1

    return [boucle,tot, mot,debut]
           
        
def geneSequences(x,chaine,alphabet):
    """_summary_

    Args:
        x (tableau de string): mots à retrouvé dans le fichier
        chaine (string): chiane à parser
        alphabet (string): alphabet des lettres étudié, aide à determiner quels caractère sont à ignoré

    Returns:
        tableau de tableau de int: les séquaneces encadré par des introns
    """

    tailleListe=len(x)
    tailleMot=len(x[0])
    
    occ=[]

    pointeur=0
    chaine[0]

    start=False
    borneInf=-1
    borneSup=-1

    while( pointeur < len(chaine)):

        [boucle,tot,mot,debut]=dansAlphabet(pointeur,chaine,alphabet,tailleMot)
        trouve=-1

        for k in range(tailleListe):
            if(x[k]==mot):
                trouve=k

        if(trouve>-1):
            pointeur=pointeur+tot

            if trouve==0 and start==False and borneInf<0:
                borneInf=pointeur
                start=True
                
            elif trouve!=0 and start==True and borneSup<0:
                borneSup=pointeur-tot
        else:
            pointeur=pointeur+1

        if borneInf>0 and borneSup>0:
            occ.append([borneInf,borneSup])
            borneSup=-1
            borneInf=-1
            start=False
        
        

    return occ

