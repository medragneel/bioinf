from math import log2
import argparse
import sys
import time

def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        if 'log_time' in kw:
            name = kw.get('log_name', method.__name__.upper())
            kw['log_time'][name] = int((te - ts) * 1000)
        else:
            print('%r  %2.2f ms' %(method.__name__, (te - ts) * 1000))
        return result
    return timed

# creer une matrice de nrows (lines) et ncols (columns)
def create_mat(nrows,ncols):
    mat = [[0.0] * ncols for _ in range(0,nrows)]
    return mat
def print_mat(mat,alphabet="ACGT"):
    for i in range(len(mat)):
        print(alphabet[i],end="\t")
        for j in range(len(mat[i])):
            print(str(mat[i][j]).rstrip("0").rstrip("."),sep=" ",end="\t")
        print()


    
def mat_counts(seqs,alphabet="ACGT"):
    counts = create_mat(len(alphabet),len(seqs[0]))
    for s in seqs:
        for i in range(len(seqs[0])):
            idx = alphabet.index(s[i])
            counts[idx][i] += 1
    return counts


def fq_mat(seqs,alphabet="ACGT"):
    c = mat_counts(seqs)
    f = create_mat(len(alphabet),len(seqs[0]))
    for s in seqs:
        for i in range(len(seqs[0])):
            idx = alphabet.index(s[i])
            f[idx][i] = float("{0:.2f}".format(c[idx][i] / len(seqs)))
    return f

def pwm(seqs,alphabet="ACGT"):
    f  = fq_mat(seqs)
    p = create_mat(len(alphabet),len(seqs[0]))
    for s in seqs:
        for i in range(len(seqs[0])):
            idx = alphabet.index(s[i])
            p[idx][i] =  float("{0:.2f}".format(log2(f[idx][i] / 0.25)))
    return p 
    

def get_motif(seqs,c,alphabet="ACGT"):
    m = ""
    for j in range(len(seqs[0])):
        max_c = c[0][j]
        max_ci = 0
        for i in range(1,len(alphabet)):
            if c[i][j] > max_c:
                max_c = c[i][j]
                max_ci =  i
        m += alphabet[max_ci] 
    return m

parser=argparse.ArgumentParser(description='Bio_mat')
parser.add_argument('-f','--file',help='read from file')
parser.add_argument('-s','--seqs',nargs="*",help='read seqs from inputs')
parser.add_argument('-o','--output',nargs="?",help='output in a file')
parser.add_argument('-i','--info',action="store_true",help='info about the program')
args= parser.parse_args()

def read_seq_form_file(file):
    with open(file,"r") as f:
        lines = [line.strip("\n") for line in f.readlines() ]
    return lines

def print_file_result():
    s = read_seq_form_file(args.file)
    c = mat_counts(s)
    f = fq_mat(s)
    p = pwm(s)
    print()
    print("matrice de contage")
    print()
    print_mat(c)
    print()
    print("matrice de frequence")
    print()
    print_mat(f)
    print()
    print("pwm")
    print()
    print_mat(p)
    print()
    print("motif")
    print()
    m= get_motif(s,p)
    print(m)
def print_stdin_result():
     c = mat_counts(args.seqs)
     f = fq_mat(args.seqs)
     p = pwm(args.seqs)
     print()
     print("matrice de contage")
     print()
     print_mat(c)
     print()
     print("matrice de frequence")
     print()
     print_mat(f)
     print()
     print("pwm")
     print()
     print_mat(p)
     print()
     print("motif")
     print()
     m= get_motif(args.seqs,p)
     print(m)





seqs2=[
        "GCCGGAAGTG",
        "ACCGGAAGCA",
        "GCCGGATGTA",
        "ACCGGAAGCT",
        "ACCGGATATA",
        "CCCGGAAGTG",
        "ACAGGAAGTC",
        "GCCGGATGCA",
        "TCCGGAAGTA",
        "ACAGGAAGCG",
        "ACAGGATATG",
        "TCCGGAAACC",
        "ACAGGATATC",
        "CAAGGACGAC",
        "TCTGGACCCT",
        ]

seqs = [
        "AAAGTT",
        "CACGTG",
        "TTGGGT",
        "GACCGT",
        "AACCAT",
        "AAACCT",
        "GAACCT"
        ]
@timeit
def test():
    if args.file:
        print_file_result()
        if args.output:
             with open(args.output,"w") as sys.stdout:
                 print_file_result()

    elif args.seqs:
        print_stdin_result()
        if args.output:
             with open(args.output,"w") as sys.stdout:
                 print_stdin_result()
    elif args.info:
        print('devInfo: created by Mahdjoub Mohamed')
        print('website: https://medmh.netlify.app/')
        print('createdat: 2-dec-2022')


if __name__ == "__main__":
    test()



# GCCGGAAGTG ACCGGAAGCA GCCGGATGTA ACCGGAAGCT ACCGGATATA CCCGGAAGTG
