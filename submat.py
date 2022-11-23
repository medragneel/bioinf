class SubMat():
    """docstring for SubMat."""
    def __init__(self):
        self.alphabet = ""
        self.sm = {}
    def __getitem__(self,ij):
        i,j=ij
        return self.score_pair(i,j)

     
    def score_pair(self,c1,c2):
        if c1 not in self.alphabet or c2 not in self.alphabet:
            return None
        return self.sm[c1+c2]
        
     
    def create_sub_mat(self,m,mm,alphabet):
        self.alphabet = alphabet
        for c1 in alphabet:
            for c2 in alphabet:
                if (c1==c2):
                    self.sm[c1+c2] = m
                else:
                   self.sm[c1+c2] = mm
        return None

    def read_sub_mat_from_file(self,filename,sep):
        with open(filename ,"r") as blosum62:
            line=blosum62.readline().split(sep)
            self.alphabet = ""
            for i in range(0,len(line)):
                self.alphabet += line[i][0]
            
            for i in range(0,len(line)):    
                ln = blosum62.readline().split(sep)
                for j in range(0,len(ln)):
                    k = self.alphabet[i] + self.alphabet[j]
                    self.sm[k] = int(ln[j])
            return None
    
def test1():
    sm = SubMat()
    sm.read_sub_mat_from_file("blosum62.mat", "\t")
    print(sm.alphabet)
    print(sm.score_pair("G", "M"))
    print(sm.score_pair("W", "W"))
    print(sm.score_pair("A", "S"))
    print(sm.score_pair("X", "X"))
    print(sm["G","K"])
    print(sm["T","T"])


def test2():
    sm = SubMat()
    sm.create_sub_mat(3, -1, "ACGU")
    print(sm.alphabet)
    print(sm.score_pair("A", "A"))
    print(sm.score_pair("A", "U"))
    print(sm.score_pair("T", "T"))
    print(sm["G","G"])


if __name__ == "__main__":   
    test1()
    print()
    test2()
