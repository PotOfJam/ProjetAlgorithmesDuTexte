from filecmp import cmp
from os import path, chdir

def compare(filelist1, filelist2):
    chdir(path.dirname(path.abspath(__file__)))

    if len(filelist1) == len(filelist2) and \
        isinstance(filelist1, list) and isinstance(filelist2, list):
        for i in range(0, len(filelist1)):
            check = cmp(filelist1[i], filelist2[i])
            print(check)
    else:
        raise ValueError("Both arguments of compare should be lists of same length!")


f1 = 'f1.txt'
f2 = 'f2.txt'
f3 = 'f3.txt'
f4 = 'f4.txt'
l1 = [f1, f3]
l2 = [f2, f4]
l3 = [f2, f4, f3]
compare(l1, l2)
# compare(l1, l3) # raise a ValueError