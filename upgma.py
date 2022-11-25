from numat import NumMatrix
from h_clustering import H_clustering
from myseq import MySeq
from pairwise import Pairwise
from submat import SubMat

class Upgma:
    def __init__(self,seqs,alseq):
        self.seqs = seqs
        self.alseq = alseq
        self.create_mat_dist()
    def create_mat_dist(self):
        self.mat_dist = NumMatrix(len(self.seqs),len(self.seqs))
        for i in range(len(self.seqs)):
            for j in range(i,len(self.seqs)):
                s1 = self.seqs[i]
                s2 = self.seqs[j]
                self.alseq.nw(s1,s2)
                alin = self.alseq.final_align()
                ncd = 0
                for k in range(len(alin)):
                    col = alin.col(k)
                    if (col[0] != col[1]):ncd +=1
                self.mat_dist.set_val(i,j,ncd)
    def run(self):
        ch = H_clustering(self.mat_dist)
        t = ch.execute_clustering()
        return t

def test():
    seq1 = MySeq("ATAGCGAT")    
    seq2 = MySeq("ATAGGCCT")    
    seq3 = MySeq("CTAGGCCC")
    seq4 = MySeq("CTAGGCCT")    
    sm = SubMat()    
    sm.create_sub_mat(1, -1, "ACGT")    
    alseq = Pairwise(sm, -2)    
    up  = Upgma([seq1, seq2, seq3, seq4], alseq)    
    arv = up.run()    
    arv.print_tree() 

if __name__ == "__main__":
    test()
