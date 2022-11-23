def print_mat(mat):
    for i in range (0,len(mat)):
        print(mat[i])
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
def score_pos(c1,c2,sm,g):
    if c1  == "−" or c2 == "−":
        return g
    else:
        return sm[c1+c2]
def max3t(v1,v2,v3):
    if v1 >v2:
        if v1 > v3: return 1
        else: return 3
    else:
        if v2 > v3: return 2
        else: return 3

def sw(s1,s2,sm,g): 
    s= [[0]]
    t= [[0]]
    maxscore = 0
    
    #for rows 
    for _ in range(1,len(s2)+1):
        s[0].append(0)
        t[0].append(0)
    #for cols  
    for _ in range(1,len(s1)+1):
        s.append([0])
        t.append([0])
    #### 
    for i in range(0,len(s1)):
        for j in range(len(s2)):
            rs1 = s[i][j] + score_pos(s1[i],s2[j],sm,g)
            rs2= s[i][j+1] + g
            rs3= s[i+1][j] + g
            b = max(rs1,rs2,rs3)
            if(b <= 0 ):
                s[i+1].append(0)
                t[i+1].append(0)
            else:
                s[i+1].append(b)
                t[i+1].append(max3t(rs1,rs2,rs3))
                if b > maxscore:
                    maxscore = b


            # s[i+1].append(max( rs1,rs2,rs3 ))
            # t[i+1].append(max( rs1,rs2,rs3 ))
    return (s,t,maxscore)


def max_mat(mat):
    maxval = mat[0][0]
    maxrow=0
    maxcol=0
    for i in range(0,len(mat)):
        for j in range(0,len(mat[i])):
            if mat[i][j] > maxval:
                maxval = mat[i][j]
                maxrow = i
                maxcol = j
    return (maxrow,maxcol)
def final_align(t,s,s1,s2):
    res=["",""]
    i,j = max_mat(s)
    while t[i][j] > 0:
        if t[i][j] == 1:
            res[0] = s1[i-1] + res[0]
            res[1] = s2[j-1] + res[1]
            i -= 1
            j -= 1
        elif t[i][j] == 3:
            res[0] = "−" + res[0]
            res[1] = s2[j-1] + res[1]
            j -= 1
        else:
            res[0] = s1[i-1] + res[0]
            res[1] = "−"  + res[1]
            i -= 1
    return res


def test_sw():
    sm = read_sub_mat_from_file('./blosum62.mat')
    s1 = "ANDDR"
    s2 = "AARRD"

    res = sw(s1,s2,sm,-8)
    S= res[0]
    T= res[1]
    print(f"Score Optimal is : {res[2]}")
    print_mat(S)
    print_mat(T)
    align = final_align(T,S,s1,s2)
    print(align[0])
    print(align[1])


test_sw()
