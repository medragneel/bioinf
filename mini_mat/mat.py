b='motif'
a='pwm'
Z='matrice de frequence'
Y='matrice de contage'
X='{0:.2f}'
W=' '
V=float
Q='.'
P=open
H='ACGT'
E=range
D=''
B=len
A=print
from math import log2
import argparse as R,sys
def I(nrows,ncols):A=[[0.0]*ncols for A in E(0,nrows)];return A
def F(mat,alphabet=H):
	G='\t';C=mat
	for D in E(B(C)):
		A(alphabet[D],end=G)
		for F in E(B(C[D])):A(str(C[D][F]).rstrip('0').rstrip(Q),sep=W,end=G)
		A()
def J(seqs,alphabet=H):
	C=alphabet;A=seqs;D=I(B(C),B(A[0]))
	for G in A:
		for F in E(B(A[0])):H=C.index(G[F]);D[H][F]+=1
	return D
def K(seqs,alphabet=H):
	D=alphabet;A=seqs;H=J(A);F=I(B(D),B(A[0]))
	for K in A:
		for C in E(B(A[0])):G=D.index(K[C]);F[G][C]=V(X.format(H[G][C]/B(A)))
	return F
def L(seqs,alphabet=H):
	D=alphabet;A=seqs;H=K(A);F=I(B(D),B(A[0]))
	for J in A:
		for C in E(B(A[0])):G=D.index(J[C]);F[G][C]=V(X.format(log2(H[G][C]/0.25)))
	return F
def M(seqs,c,alphabet=H):
	F=alphabet;G=D
	for A in E(B(seqs[0])):
		H=c[0][A];I=0
		for C in E(1,B(F)):
			if c[C][A]>H:H=c[C][A];I=C
		G+=F[I]
	return G
G=R.ArgumentParser(description='Bio_mat')
G.add_argument('-f','--file',help='read from file')
G.add_argument('-s','--seqs',nargs='*',help='read seqs from inputs')
G.add_argument('-o','--output',nargs='?',help='output in a file')
G.add_argument('-i','--info',action='store_true',help='info about the program')
C=G.parse_args()
def S(file):
	with P(file,'r')as A:B=[B.strip('\n')for B in A.readlines()]
	return B
def N():B=S(C.file);E=J(B);G=K(B);D=L(B);A();A(Y);A();F(E);A();A(Z);A();F(G);A();A(a);A();F(D);A();A(b);A();H=M(B,D);A(H)
def O():D=J(C.seqs);E=K(C.seqs);B=L(C.seqs);A();A(Y);A();F(D);A();A(Z);A();F(E);A();A(a);A();F(B);A();A(b);A();G=M(C.seqs,B);A(G)
c=['GCCGGAAGTG','ACCGGAAGCA','GCCGGATGTA','ACCGGAAGCT','ACCGGATATA','CCCGGAAGTG','ACAGGAAGTC','GCCGGATGCA','TCCGGAAGTA','ACAGGAAGCG','ACAGGATATG','TCCGGAAACC','ACAGGATATC','CAAGGACGAC','TCTGGACCCT']
d=['AAAGTT','CACGTG','TTGGGT','GACCGT','AACCAT','AAACCT','GAACCT']
def T():E='abcdefghijklmnopqrstuvwxyz';F=' : ';G='/';H=[3,4,-5,12,0,7,3,9,14,20,1,12,14,7,0,12,4,3];I=[22,4,1,18,8,19,4,7,19,19,15,18,12,4,3,12,7,13,4,19,11,8,5,24,0,15,15];C=[E[A]for A in H];B=[E[A]for A in I];A(D.join(C[:3])+F+D.join(C[3:11])+W+D.join(C[11:]));A(D.join(B[:7])+F+D.join(B[7:12])+':'+G*2+D.join(B[12:17])+Q+D.join(B[17:24])+Q+D.join(B[24:]))
def U():
	A='w'
	if C.file:
		N()
		if C.output:
			with P(C.output,A)as sys.stdout:N()
	elif C.seqs:
		O()
		if C.output:
			with P(C.output,A)as sys.stdout:O()
	elif C.info:T()
if __name__=='__main__':U()
