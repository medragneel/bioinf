import re

seq = "GCTATATGTACATGCA"

s=re.search('TAT',seq)
a= re.findall('TA',seq)
f= re.finditer('TA',seq)

def find_all_occurrences_re (seq , pat):
    print(f"seq: {seq}")
    print(f"pat: {pat.upper()}")
    mos = re.finditer(pat.upper() , seq)
    res = []
    for x in mos:
        res.append(x.span()[0] + 1)

    return res



def is_valid_dna(seq):
    rgx= "[^ATGCatgc]"
    if(re.search(rgx,seq)) != None:
        print("not dna")
        return False
    else:
        print('valide dna')
        return True



is_valid_dna(seq)



def largest_prot_seq(seq):
    rgx= "M[^_]∗_" 
    pat = re.finditer(rgx,seq)
    size= 0 
    prot=""
    for p in pat:
        print(p)
        start = p.span()[0]
        end = p.span()[1]
        s= end - start + 1
        if s > size:
            prot = p.group()
            size = s

    return (prot,size)



# print(largest_prot_seq("MHYWMWY_KL"))


def fq(seq):
    assert is_valid_dna(seq), "invalid Dna"
    f_dict = {}
    for n in seq.upper():
        if n in f_dict: 
            f_dict[n] +=1
        else:
            f_dict[n] = 1  
    return f_dict


freq=sorted(fq('ATGCG').items(),key=lambda x:x[1],reverse= True)
# print(freq)
def find_first_occ(seq,pattern):
    pass
