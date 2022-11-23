
def create_sub_mat(m,mm,alphabet):
    sm={}
    for c1 in alphabet:
        for c2 in alphabet:
            if (c1==c2):
                sm[c1+c2] = m
            else:
                sm[c1+c2] = mm
    return sm

sub = create_sub_mat(alphabet="ACGT",m=1,mm=0)
def read_sub_mat_from_file(filename):
    sm={}
    with open(filename ,"r") as blosum62:
        line=blosum62.readline().split("\t")
        alphabet=[ line[i][0] for i in range(0,len(line))]
        for i in range(0,len(line)):
            ln = blosum62.readline().split("\t")
            for j in range(0,len(ln)):
                k = alphabet[i] + alphabet[j]
                sm[k] = int(ln[j])
        return sm


# read_sub_mat_from_file('./blosum62.mat')

def score_pos(c1,c2,sm,g):
    if c1  == "−" or c2 == "−":
        return g
    else:
        return sm[c1+c2]
    
def score_al(s1,s2,sm,g):
    res=0
    for i in range(len(s1)):
        res += score_pos(s1[i],s2[i],sm,g)
    return res
    
def improved_score_al(s1,s2,sm,g,r):
    res = 0
    in_g1=False
    in_g2=False
    for i in range(len(s1)):
        if s1[i] == "−":
            if in_g1 : res+= r
            else:
                in_g1=True
                res += g
        elif s2[i] == "−":
            if in_g2 : res += r
            else:
                in_g2=True
                res += g
        else:
            if in_g1: in_g1 = False
            if in_g2: in_g2 = False
            res += sm[s1[i]+s2[i]]
    return res





def test_dna():
    sm=create_sub_mat(m=2,mm=-2,alphabet="AGCT")
    s1 = "−CAGTGCATG−ACATA"
    s2="TCAG−GC−TCTACAGA"
    g=-3
    print(score_al(s1,s2,sm,g))

def test_prot():
    sm= read_sub_mat_from_file("./blosum62.mat")
    s1="LGPSSGCASRIWTKSA"
    s2="TGPS−G−−S−IWSKSG"
    g=-8
    r=-2
    print(improved_score_al(s1,s2,sm,g,r))
    print(score_al(s1,s2,sm,g))
    



test_dna()
test_prot()
