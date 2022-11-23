class Align():
    def __init__(self,ls,al_type="P") -> None:
        self.ls = ls
        self.al_type = al_type
        self.alphabet = ""
        self.sm = {}
    def __len__(self):
        return len(self.ls[0])
    def __getitem__(self,n):
        if type(n) is tuple and len(n) == 2:
            i,j = n
            return self.ls[i][j]
        elif type(n) is int: return self.ls[n]
        return None
    def __str__(self):
        res=""
        for sq in self.ls:
            res += "\n" + sq
        return res
    def num_seqs(self):
        return len(self.ls)
    def col(self,idx):
        return [self.ls[k][idx] for k in range(len(self.ls))]
    def consensus(self):
        cons=""
        for i in range(len(self)):
            cont={}
            for k in range(len(self.ls)):
                c = self.ls[k][i]
                if c in cont:
                    cont[c] = cont[c] + 1
                else:
                    cont[c] = 1
            maximum=0
            cmax=None
            for key in cont.keys():
                if key != "−" and cont[key] > maximum:
                    maximum= cont[key]
                    cmax= key
            cons =  cons + cmax
        return cons




def test_align():
    al= Align(ls=["ATGA−A","AA−AT−"],al_type="D")
    print(len(al))
    print (al.__len__())
    print(al.col(idx=1))
    print(al[1,1])
    print(al.consensus())
        

