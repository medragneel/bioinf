from ore_no_align import Align
from myseq import MySeq
from submat import SubMat
from pairwise import Pairwise

class Msa:
    """docstring for Msa."""
    def __init__(self,seqs,alignseq):
        self.seqs = seqs
        self.alignseq = alignseq
    def add_seq_alignement(self,alignement,seq):
        res=[]
        for _ in range(len(alignement.ls())+1):
            res.append("")
        cons = MySeq(alignement.consensus(),alignement.al_type)
        self.alignseq.nw(cons,seq)
        align2 = self.alignseq.final_align()
        orig = 0
        for i in range(len(align2)):
            if align2[0,i] == '−':
                for k in range(len(alignement.ls)):
                    res[k] += "−"
            else:
                for k in range(len(alignement.ls)):
                    res[k] += alignement[k,orig] 
                orig +=1
        res[len(alignement.ls) ] = align2.ls[1]
        return Align(res,alignement.al_type)
    def align_consensus(self):
        self.alignseq.nw(self.seqs[0],self.seqs[1])
        res = self.alignseq.final_align()

        for i in range(2,len(self.seqs)):
            res = self.add_seq_alignement(res,self.seqs[i])
        return res



def printMat (mat):
    for i in range(0, len(mat)):
        print(mat[i])
def test_prot():  
    s1 = MySeq("PHWAS","protein")
    s2 = MySeq("HWASW","protein")
    s3 = MySeq("HPHWA","protein")
    sm = SubMat()
    sm.read_sub_mat_from_file("blosum62.mat", "\t")
    aseq = Pairwise(sm, -8)
    ma = Msa([s1,s2,s3], aseq)
    alinm = ma.align_consensus()
    print(alinm)
    

def test():
    s1 = MySeq("ATAGC")
    s2 = MySeq("AACC")
    s3 = MySeq("ATGAC")
    
    sm = SubMat()
    sm.create_sub_mat(1,-1,"ACGT")
    aseq = Pairwise(sm,-1)
    ma = Msa([s1,s2,s3], aseq)
    al = ma.align_consensus()
    print(al)
    
def exercise1():
    s1 = MySeq("ACATATCAT")
    s2 = MySeq("AACAGATCT")
    s3 = MySeq("AGATATTAG")
    s4 = MySeq("GCATCGATT")
    
    sm = SubMat()
    sm.create_sub_mat(1,-1,"ACGT")
    aseq = Pairwise(sm,-1)
    ma = Msa([s1,s2,s3,s4], aseq)
    al = ma.align_consensus()
    print(al)

if __name__ == "__main__": 
    test_prot()
    print()
    test()
    print()
    exercise1()
