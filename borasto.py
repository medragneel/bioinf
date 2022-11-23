from typing import SupportsAbs


def read_db(filename):
    with open(filename,"r") as f:
        db = [line.rstrip() for line in f]
    return db 


def build_map(q,w):
    res ={}
    for i in range(len(q)-w+1):
        sub_seq = q[i:i+w]
        if sub_seq in res:
            res[sub_seq].append(i)
        else:
            res[sub_seq]= [i]
    return res

def get_matches(seq,m,w):
    mts= []
    for i in range(len(seq)-w+1):
        sub_seq = seq[i:i+w] 
        if sub_seq in m:
            l= m[sub_seq]
            for ind in l:
                mts.append((ind,i))
    return mts
    
def extends_match(seq,hit,q,w):
    stq,sts = hit[0],hit[1]
    
    matfw = 0
    k=0
    bestk=0
    while 2 * matfw >= k and stq+w+k < len(q) and sts+w+k < len(seq):
        if q[stq+w+k] == seq[sts+w+k]:
            matfw +=1
            bestk = k+1
    #move backward
    size = w + bestk
    k=0
    matbw = 0
    bestk = 0
    while 2 * matfw >= k and stq > k and sts > k:
        if q[stq-k-1] == seq[sts-k-1]:
            matbw +=1
            bestk = k+1
        k +=1
    size += bestk
    return (stq-bestk,sts-bestk,size,w+matfw+matbw)

def hit_best_score(seq,q,m,w):
    hits=get_matches(seq,m,w)
    best_score = -1.0
    best=()
    for h in hits:
        ext = extends_match(seq,h,q,w)
        score = ext[3]
        if score > best_score or (score == best_score and ext[2] < best[2]):
            best_score = score
            best = ext
    return best

def best_alignement(db,q,w):
    m = build_map(q,w)
    best_score = -1.0
    res = (0,0,0,0)
    for k in range(0,len(db)):
        best_seq = hit_best_score(db[k],q,m,w)
        if best_seq !=():
            score = best_seq[3]
            if score > best_score or (score == best_score and best_seq[2] < res[2]):
                best_score = score
                res = best_seq[0], best_seq[1],best_seq[2],best_seq[3],k
            if best_score < 0: return ()
            else: return res


db = read_db("seqBlast.txt")
query = "gacgcctcgcgctcgcgcgctgaggcaaaaaaaaaaaaaaaaaaaatcggatagctagctgagcgctcgatagcgcgttcgctgcatcgcgtatagcgctgaagctcccggcgagctgtctgtaaatcggatctcatctcgctctatcct"
r = best_alignement(db, query, 11)
print(r)

#query2 = "cgacgacgacgacgaatgatg"
#r = best_alignement(db, query2, 11)
#print(r)
