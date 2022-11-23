class Boyer:
    def __init__(self,alphabet,pat) -> None:
       self.alphabet = alphabet 
       self.pat = pat
       self.preprocess()
    def preprocess(self):
        self.process_bcr()
        self.process_gsr()
    def process_bcr(self):
        self.occ = {}
        for s in self.alphabet:
            self.occ[s] = -1
        for j in range(len(self.pat)):
            # print(j)
            c = self.pat[j]
            self.occ[c] = j
        print(self.occ)
    def process_gsr(self):
        self.f = [0] * (len(self.pat) + 1)
        self.s = [0] * (len(self.pat) + 1)
        i = len(self.pat)
        j = len(self.pat) + 1
        self.f[i] = j
        while i > 0 :
            while j <= len(self.pat) and self.pat[i-1] != self.pat[j-1]:
                if self.s[j] == 0 : self.s[j] = j-i
                j = self.f[j]
                j -= 1
                i -= 1
                self.f[i] = j
                self.f[0]
                for i in range(len(self.pat)):
                    if self.s[i] ==0 : self.s[i] = j
                    if i == j: j= self.f[j]
                print(self.f)
                print(self.s)
                     
    def search_pattern(self, text):
         res = []
         i = 0
         while i <= len(text) - len(self.pat):
            j= len(self.pat)- 1

            while j>=0 and self.pat[j]==text[j+i]: 
                j -= 1 
                if (j<0):
                    res.append(i)
                    i += self.s[0]
                else:
                    c = text[j+i]            
                    i += max(self.s[j+1], j- self.occ[c])

            return res




def test_boyer():
    b= Boyer(alphabet="ATGC",pat="TATA")
    print(b.search_pattern('ATGTATATAGCTGCTACT'))

test_boyer()
    
    
        
