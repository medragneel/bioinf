import re
import sys
def find_zinc_finger(seq):
    rgx= "C.H.[LIVMFY]C.{2}C[LIVMYA]"
    mo= re.search(rgx,seq)
    if(mo!= None):
        print(mo.span())
        return mo.span()[0]
    else:
        return -1

def prosite(sq,profile):
    rgx = profile.replace("−","")
    rgx = rgx.replace("x",".")
    rgx = rgx.replace("(","{")
    rgx = rgx.replace(")","}")
    print(f"Prosite: { profile }")
    print(f"Regex: { rgx }")
    mo = re.search(rgx,sq)
    if(mo != None):
        return mo.span()[0]
    else:
        return -1


def check_iub_code(iub):
    codes={
            "A":"A",
            "T":"T",
            "G":"G",
            "C":"C",
            "K": "[TG]",
            "M": "[AC]",
            "R": "[GA]",
            "y": "[TC]",
            "S": "[GC]",
            "W": "[AT]",
            "B": "[TGC]",
            "D": "[ATG]",
            "H": "[ATC]",
            "V": "[AGC]",
            "N": "[ATCG]"

            }
    site = iub.replace("^","")
    rgx=""
    for c in site:
        rgx += codes[c]
    return rgx

def cut_pos(enz,sq):
    pos= enz.find("^")
    #the index of the caret in the enz string
    print(pos)
    # return a regex expression
    rgx= check_iub_code(enz)
    print(rgx)

    # return an iterable of the matches
    matches = re.finditer(rgx,sq)
    locs=[]
    for m in matches:
        locs.append(m.start() + pos)

    # return the locactions of the matches
    return locs

def cut_sub_sq(locs,sq):
    res=[]
    positions= locs
    positions.insert(0,0)
    positions.append(len(sq))
    print(positions)
    for i in range(len(positions)-1):
        # cut the seq to part
        # i=0 i+1= 1 sq[i] = [sq[0]=0 : sq[1] = 7 ]
        # i=1 i+1= 2 sq[i] = [sq[1]=7 : sq[2] = 19 ]
        # i=2 i+1= 3 sq[i] = [sq[2]=19 : sq[3] = 23 ]
        res.append(sq[positions[i]:positions[i+1]])
    return res


def create_mat(nrow,ncol):
    mat = []
    for i in range(nrow):
        mat.append([])
        for j in range(ncol):
            mat[i].append(0)
    return mat

def dotplot(seq1,seq2):
    mat = create_mat(len(seq1),len(seq2))
    for i in range(len(seq1)):
        for j in range(len(mat)):
            if seq1[i] == seq2[j]:
                mat[i][j] = 1
    return mat

def print_dplot(m,s1,s2):
    sys.stdout.write("  " + s2 + "\n")
    for i in range(len(m)):
        sys.stdout.write(s1[i])
        for j in range( len(m[i] )):
            if m[i][j] >= 1:
                sys.stdout.write(' *')
            else:
                sys.stdout.write(" ")
        sys.stdout.write("\n")

s1 = "ATG"
s2 = "ATG"
mat = dotplot(s1,s2)
print_dplot(mat,s1,s2)

# just a test function 
def test():
    # sq = "HKMMLASCKHLLCLKCIVKLG"
    # print(f"seq Length: { len(sq) }")
    # print(prosite(sq=sq,profile="C−x−H−x−[LIVMFY]−C−x(2)−C−[LIVMYA]"))
    # print(check_iub_code("G^AMTV"))
    pos = cut_pos(enz="G^ATTC",sq="GTAGAAGATTCTGAGATCGATTC")
    print(cut_sub_sq(pos,"GTAGAAGATTCTGAGATCGATTC"))

# test()

