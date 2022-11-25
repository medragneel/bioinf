from myseq import MySeq
class MotifFinder:
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
        res = [[0]*self.motif_size for _ in range(len(self.alphabet))]
        for i,val in enumerate(indexes):
            subseq = self.seqs[i][val:(val + self.motif_size)]
            for ms in range(self.motif_size):
                for k in range(len( self.alphabet )):
                    if subseq[ms] == self.alphabet[k]:
                        res[k][ms] = res[k][ms] + 1
        return res

    def score(self,s):
        score = 0
        mat = self.create_motif_from_indexes(s)
        for j in range(len(mat[0])):
            maxcol = mat[0][j]
            for i in range(1,len(mat)):
                if mat[i][j] > maxcol:
                    maxcol = mat[i][j]
            score += maxcol
        return score

    def dup_score(self,s):
        score = 1.0
        mat = self.create_motif_from_indexes(s)
        for j in range(len(mat[0])):
            maxcol = mat[0][j]
            for i in range(1,len(mat)):
                if mat[i][j] > maxcol:
                    maxcol = mat[i][j]
            score *= maxcol
        return score
         
    def next_sol(self,s):
        next_sol = [0] * len(s)
        pos = len(s) - 1
        while pos >= 0 and s[pos] == self.seq_size(pos) - self.motif_size:
             pos -= 1 
        if pos < 0:
            next_sol = None
        else:
            for i in range(pos):
                next_sol[i] = s[i]
            next_sol[pos] = s[pos] + 1
            for i in range(pos+1,len(s)):
                next_sol[i] = 0
                
        return next_sol

    def exhaustive_search(self):
        best_score = -1
        res = []
        s = [0] * len(self.seqs)
        while s != None:
            sc = self.score(s)
            if sc > best_score:
                best_score = sc
                res = s
            s = self.next_sol(s) 
        return res
    def next_vertex(self,s):
        res = []
        if len(s) < len(self.seqs):
            for i in range(len(s)):
                res.append(s[i])
            res.append(0)
            #bypass
        else:
            pos = len(s) - 1
            while pos >= 0 and s[pos] == self.seq_size(pos) - self.motif_size :
                pos -= 1
            if pos < 0: res = None
            else:
                for i in range(pos): res.append(s[i])
                res.append(s[pos] + 1)
        return res

    def bypass(self,s):
        res = []
        pos = len(s) - 1
        while pos >= 0 and s[pos] == self.seq_size(pos) - self.motif_size :
            pos -= 1
        if pos < 0: res = None
        else:
            for i in range(pos): res.append(s[i])
            res.append(s[pos] + 1)
        return res
    def branch_and_bound(self):
        best_score = -1
        best_motif = None
        size = len(self.seqs)
        s = [0] * size
        while s != None:
            if len(s) < size:
                opt_score = self.score(s) + (size-len(s)) * self.motif_size
                if opt_score < best_score: s = self.bypass(s)
                else: s = self.next_vertex(s)
            else:
                sc = self.score(s)
                if sc > best_score:
                    best_score = sc
                    best_motif = s
                s = self.next_vertex(s)
        return best_motif

    def heuristic_consensus(self):
        res = [0] * len(self.seqs)
        max_score = -1
        partial = [0,0]
        for i in range(self.seq_size(0) - self.motif_size):
            for j in range(self.seq_size(1) - self.motif_size):
                partial[0] = i 
                partial[1] = j
                sc = self.score(partial)
                if (sc > max_score):
                    max_score = sc
                    res[0] = i
                    res[1] = j
        for k in range(2,len(self.seqs)):
            partial = [0] * (k+1)
            for j in range(k):
                partial[j] = res[j]
            max_score  = -1
            for i in range(self.seq_size(k)-self.motif_size):
                partial[k] = i
                sc = self.score(partial)
                if (sc > max_score):
                    max_score = sc 
                    res[k] = i
            return res



         
         
    
def test():
    m = MotifFinder()
    m.read_file("./seqBlast.txt","DNA")
    sm = m.create_motif_from_indexes([25,20,2,55,59])
    print(m.score([25,20,2,55,59]))
    print(m.dup_score([25,20,2,55,59]))
    print(sm)


    for i in range(len(sm)):
        print("ACGT"[i],end="\t")
        for j in range(len(sm[i])):
            print(sm[i][j],end="\t")
        print()
def test2():
    seq1 = MySeq("ATAGAGCTGA","DNA")
    seq2 = MySeq("ACGTAGATGA","DNA")
    seq3 = MySeq("AAGATAGGGG","DNA")
    mf = MotifFinder(3, [seq1,seq2,seq3])
    
    print ("Exhaustive:")
    sol = mf.exhaustive_search()
    print ("Solution: " , sol)
    print ("Score: ", mf.score(sol))
    

    print ("\nBranch and Bound:")
    sol2 = mf.branch_and_bound()
    print ("Solution: " , sol2)
    print ("Score:" , mf.score(sol2))
    
    print ("\nHeuristic consensus: ")
    sol3 = mf.heuristic_consensus()
    print ("Solution: " , sol3)
    print ("Score:" , mf.score(sol3))
def test3():
    mf = MotifFinder()
    mf.read_file("./example_motifs.txt","DNA")
    print ("Branch and Bound:")
    sol = mf.branch_and_bound()
    print ("Solution: ", sol)
    print ("Score:", mf.score(sol))
    
    print ("\nHeuristic consensus: ")
    sol2 = mf.heuristic_consensus()
    print ("Solution: " , sol2)
    print ("Score:" , mf.score(sol2)) 


if __name__ == "__main__":
    test()
    test2()
    test3()

        
