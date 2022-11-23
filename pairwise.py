from ore_no_align import Align
from myseq import MySeq
from submat import SubMat


class Pairwise:
    def __init__(self,sm,g):
        self.sm = sm 
        self.g = g
        self.S = {}
        self.T = {}
        self.s1 = ""
        self.s2 = ""
     
    def score_pos(self,c1,c2):
        if c1  == "−" or c2 == "−":
            return self.g
        else:
            return self.sm[c1+c2]
    def score_alin (self, al):
        res = 0;
        for i in range(len(al)):
            res += self.score_pos (al[0][i], al[1][i])
        return res
    def nw(self,s1,s2):
        if (s1.seq_type != s2.seq_type): return None
        self.S = [[0]]
        self.T = [[0]]
        self.s1 = s1
        self.s2 = s2
        for j in range(1,len(s2) + 1):
            self.S[0].append(self.g * j)
            self.T[0].append(3)
        for i in range(1,len(s1) + 1):
            self.S.append([self.g * i])
            self.T.append([2])
        for i in range(0,len(s1)):
            for j in range(len(s2)):
                rs1= self.S[i][j] + self.score_pos(s1[i],s2[j])
                rs2 = self.S[i][j+1] + self.g
                rs3 = self.S[i+1][j] + self.g
                self.S[i+1].append(max(rs1,rs2,rs3))
                self.T[i+1].append(max3t(rs1,rs2,rs3))
        return self.S[len(s1)][len(s2)]
    def final_align(self):
        res=["",""]
        i = len(self.s1)
        j= len(self.s2)
        while i > 0 or j > 0:
            if self.T[i][j] == 1:
                res[0] =  self.s1[i-1] + res[0]
                res[1] = self.s2[j-1] + res[1]
            elif self.T[i][j] == 3:
                res[0] = "-" + res[0]
                res[1] = self.s2[j-1] + res[1] 
                j -= 1
            elif self.T[i][j] ==2:
                res[0] = self.s1[i-1] + res[0]
                res[1] = "-" + res[1]
                i -= 1
            else:break
        return Align(res,self.s1.seq_type)
    def sw(self,s1,s2):
        if (s1.seq_type != s2.seq_type): return None
        self.S = [[0]]
        self.T = [[0]]
        self.s1 = s1
        self.s2 = s2
        maxscore = 0 
        # column insirtion
        for _ in range(1,len(s2) + 1 ):
            self.S[0].append(0)
            self.T[0].append(0)
        #row shit
        for _ in range(1,len(s1) + 1):
            self.S.append([ 0 ])
            self.T.append([ 0 ])
        for i in range(0,len(s1)):
            for j in range(0,len(s2)): 
                rs1 = self.S[i][j] + self.score_pos(s1[i],s2[j])
                rs2 = self.S[i][j+1] + self.g
                rs3 = self.S[i+1][j] + self.g
                b= max(rs1,rs2,rs3)
                if(b <= 0):
                    self.S[i+1].append(0)
                    self.T[i+1].append(0)
                else:
                    self.S[i+1].append(b)
                    self.T[i+1].append(max3t(rs1,rs2,rs3))
                    if b > maxscore:
                        maxscore = b
        return maxscore



    def recover_align_local(self):
        res = ["",""]
        maxscore = 0
        maxrow = 0
        maxcol = 0
        for i in range(1,len(self.S)):
            for j in range(1,len(self.S[i])):
                if self.S[i][j] > maxscore:
                    maxscore = self.S[i][j]
                    maxrow = i
                    maxcol = j
        i = maxrow
        j=maxcol
        while i > 0 or j > 0:
            if self.T[i][j]==1:
                res[0] = self.s1[i-1] + res[0]
                res[1] = self.s2[j-1] +res[1]
                i -= 1
                j-= 1
            elif self.T[i][j] == 3 :
                res[0] = "-" + res[0]
                res[1] = self.s2[j-1] +res[1]
                j-=1
            elif self.T[i][j] == 2:
                res[0] = self.s1[i-1] +res[0]
                res[1] = "-" + res[1]
                i -=1
            else: break
        return Align(res,self.s1.seq_type)





def max3t (v1, v2, v3):
    if v1 > v2:
        if v1 > v3: return 1
        else: return 3
    else:
        if v2 > v3: return 2
        else: return 3


def print_mat(mat):
    for i in range (0,len(mat)):
        print(mat[i])

def test():
    seq1 = MySeq("ATGATATGATGATT")
    seq2 = MySeq("GATGAATAGATGTGT")
    sm = SubMat()
    sm.create_sub_mat(3, -1, "ACGT")
    alin = Pairwise(sm, -3)
    print("------- smith and waterman------------")
    print(alin.sw(seq1, seq2))
    print_mat(alin.S)
    print(alin.recover_align_local())
    
    print("------- needleman and wunch -------------")
    print(alin.nw(seq1,seq2))
    print_mat(alin.S)
    print(alin.final_align())
    

if __name__ == "__main__":   

    test()

