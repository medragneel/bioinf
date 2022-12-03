
from myseq import MySeq

def create_mat(nrows,ncols):
    # for _ in range(0,nrows):
        # res.append([0]*ncols)
    res = [[0] * ncols for _ in range(0,nrows)]
    return res
def print_mat(mat,alphabet="ACGT"):
    for i in range(len(mat)):
        print(alphabet[i],end="\t")
        for j in range(len(mat[i])):
            print(mat[i][j],end="\t")
        print()


# mat = create_mat(2,3)
# print_mat(mat)

class Ore_no_motifs:
    def __init__(self,seqs=[],pwm=[],alphabet=None):
        if seqs:
            self.size = len( seqs[0] )
            self.seqs = seqs
            self.alphabet = seqs[0].alphabet() 
            self.do_count()
            self.create_pwm()
        else:
            self.pwm = pwm
            self.size = len(pwm[0])
            self.alphabet = alphabet
    ###### functions ##########
    def __len__(self):
        return self.size

    def do_count(self):
        self.counts = create_mat(len(self.alphabet),self.size)
        for s in self.seqs:
            for i in range(self.size):
                idx = self.alphabet.index(s[i])
                self.counts[idx][i] += 1
    def create_pwm(self):
        if self.counts == None: self.do_count()
        self.pwm = create_mat(len(self.alphabet),self.size)
        for i in range(len(self.alphabet)):
            for j in range(self.size):
                self.pwm[i][j] = float(self.counts[i][j]) / len(self.seqs)
    def consensus(self):
        res=""
        for j in range(self.size):
            max_col = self.counts[0][j]
            max_coli = 0
            for i in range(1,len(self.alphabet)):
                if self.counts[i][j] > max_col:
                    max_col = self.counts[i][j]
                    max_coli = i
            res += self.alphabet[max_coli]
        return res
    def masked_consensus(self):
        res=""
        for j in range(self.size):
            max_col = self.counts[0][j]
            max_coli = 0
            for i in range(1,len(self.alphabet)):
                if self.counts[i][j] > max_col:
                    max_col = self.counts[i][j]
                    max_coli = i
            if max_col > len(self.seqs) /2:
                res += self.alphabet[max_coli]
            else:
                res += "-"
        return res
    def probability_seq(self,seq):
        res = 1.0
        for i in range(self.size):
            idx = self.alphabet.index(seq[i])
            res *= self.pwm[idx][i]
        return res
    def probability_all_pos(self,seq): 
       res= [self.probability_seq(seq)  for _ in range(len(seq) - self.size +1)]
       return res
    def most_probable_seq(self,seq):
        maximum = -1.0
        max_index = -1
        for k in range(len(seq)-self.size):
            p = self.probability_seq(seq[k:k + self.size])
            if p > maximum:
                maximum = p
                max_index = k
        return max_index
    def create_motif(self,seqs):
        l = []
        for s in seqs:
            ind = self.most_probable_seq(s.seq)
            subseq = MySeq(s[ind:(ind +self.size)],s.get_seq_biotype())
            l.append(subseq)

        return  Ore_no_motifs(l)




def test():
    seq1 = MySeq("AAAGTT")
    seq2 = MySeq("CACGTG")
    seq3 = MySeq("TTGGGT")
    seq4 = MySeq("GACCGT")
    seq5 = MySeq("AACCAT")
    seq6 = MySeq("AACCCT")
    seq7 = MySeq("AAACCT")
    seq8 = MySeq("GAACCT")
    lseqs = [seq1, seq2, seq3, seq4, seq5, seq6, seq7, seq8]
    motifs = Ore_no_motifs(lseqs)
    print("alphabet")
    print(motifs.alphabet)

    print ("Counts matrix")
    print_mat(motifs.counts,motifs.alphabet)

    print("PWM")
    print_mat(motifs.pwm ,motifs.alphabet)

    for s in lseqs:
        print (s)
    print ("Consensus sequence")
    print(motifs.consensus())

    print ("Masked Consensus sequence")
    print(motifs.masked_consensus())

    print("Probabilty of Sequences")
    print(motifs.probability_seq("AAACCT"))
    print(motifs.probability_seq("ATACAG"))
    print(motifs.most_probable_seq("CTATAAACCTTACATC"))


    s1 = MySeq("TAAAGTTATGA")
    s2 = MySeq("ATGACACGTG")
    s3 = MySeq("TTTGGGTAT")

    newmotif = motifs.create_motif([s1,s2,s3])
    print_mat(newmotif.counts)
    print(newmotif.most_probable_seq("AAAACT"))
    




if __name__ == "__main__":
    test()
