from myseq import MySeq
from ore_no_motifs import Ore_no_motifs
from random import randint,random

class MotifFinding:
    def __init__(self,size=8,seqs=None):
        self.motif_size = size
        if seqs != None:
            self.seqs = seqs
            self.alphabet = seqs[0].alphabet()
        else:
            self.seqs = []
            self.alphabet = "ACGT" 
    def __len__(self):
        return len(self.seqs)
    def __getItem__(self,n):
        return self.seqs[n]
    def seq_size(self,i):
        return len(self.seqs[i])
    def read_file(self,file,t):
        with open(file,"r") as f:
            self.seqs = [ MySeq( s.strip().upper(),t ) for s in f.readlines() ] 
            self.alphabet = self.seqs[0].alphabet()
            print(self.alphabet)
    def create_motif_from_indexes(self,indexes):
        pseqs=[]
        for i,val in enumerate(indexes):
            pseqs.append(MySeq(self.seqs[i][val:(val + self.motif_size)]))
        return Ore_no_motifs(pseqs)
            

    def score(self,s):
        score = 0
        motif = self.create_motif_from_indexes(s)
        motif.do_count()
        mat = motif.counts
        for j in range(len(mat[0])):
            maxcol = mat[0][j]
            for i in range(1,len(mat)):
                if mat[i][j] > maxcol:
                    maxcol = mat[i][j]
            score += maxcol
        return score

    def dup_score(self,s):
        score = 1.0
        motif = self.create_motif_from_indexes(s)
        motif.create_pwm()
        mat = motif.pwm
        for j in range(len(mat[0])):
            maxcol = mat[0][j]
            for i in range(1,len(mat)):
                if mat[i][j] > maxcol:
                    maxcol = mat[i][j]
            score *= maxcol
        return score
    def heuristic_stochastic(self):
        s = [0] * len(self.seqs)
        for k in range(len(self.seqs)):
            s[k] = randint(0,self.seq_size(k)-self.motif_size)
        motif = self.create_motif_from_indexes(s)
        motif.create_pwm()
        sc = self.dup_score(s)
        best_sol = s
        improve = True
        while improve:
            for k in range(len(s)):
                s[k] = motif.most_probable_seq(self.seqs[k])
            if self.dup_score(s) > sc:
                sc = self.dup_score(s)
                best_sol = s
                motif = self.create_motif_from_indexes(s)

                motif.create_pwm()
            else: improve = False
        return best_sol
    def gibbs(self,iterations=100):
        s = []
        for _ in range(len(self.seqs)):
           s.append(randint(0,len(self.seqs[0])-self.motif_size - 1))
        best_s = list(s)
        best_score = self.dup_score(s)
        for _ in range(iterations):
            seq_idx = randint(0,len(self.seqs)-1)
            seq_sel = self.seqs[seq_idx]
            s.pop(seq_idx)
            removed = self.seqs.pop(seq_idx)
            motif = self.create_motif_from_indexes(s)
            motif.create_pwm()
            self.seqs.insert(seq_idx,removed)
            r = motif.probability_all_pos(seq_sel)
            pos = self.roulette(r)
            s.insert(seq_idx,pos)
            score = self.dup_score(s)
            if score > best_score:
                best_score = score
                best_s = list(s)
        return best_s


    def roulette(self, f):
        tot = 0.0
        for x in f: tot += (0.01 + x) 
        val = random() * tot
        acum = 0.0
        idx = 0
        while acum < val:
            acum += (f[idx] + 0.01)
            idx += 1
        return idx-1




        

         

         
         
    
def test():
    sm = MotifFinding()
    sm.read_file("example_motifs.txt","dna")
    sol = [25,20,2,55,59]
    sa = sm.score(sol)
    print(sa)
    scm = sm.dup_score(sol)
    print(scm)



def test1():
    mf = MotifFinding()
    mf.read_file("example_motifs.txt","dna")
    print("Heuristic stochastic")
    sol = mf.heuristic_stochastic()
    print ("Solution: " , sol)
    print ("Score:" , mf.score(sol))
    print ("Score mult:" , mf.dup_score(sol))
    print("Consensus:", mf.create_motif_from_indexes(sol).consensus()) 
    
    print("gibbs")
    sol2 = mf.gibbs(10000)
    print ("Score:" , mf.score(sol2))
    print ("Score mult:" , mf.dup_score(sol2))
    print("Consensus:", mf.create_motif_from_indexes(sol2).consensus())




if __name__ == "__main__":
    test()
    test1()

        
