def read_fasta_file(filename,sc=">"):
    fileobj = open(filename,"r")
    seqs = []
    seq = ''
    for line in fileobj:
        if line.startswith(sc):
            if seq:
                seqs.append(seq)
            seq = ''
        else:
            seq += line.rstrip()
    if seq:
        seqs.append(seq)
    fileobj.close()
    return seqs

def read_pdb_file(filename):
    fileobj = open(filename,'r')
    natoms = 0
    xsum = ysum = zsum = 0
    for l in fileobj:
        if l[:6] == 'ATOM ':
            natoms += 1
            x = float(l[30:38])
            y = float(l[38:46])
            z = float(l[46:54])
            xsum += x 
            ysum += y
            zsum += z
    fileobj.close()
    if natoms == 0:
        xavg = yavg = zavg = 0
    else:
        xavg = xsum /   natoms 
        yavg = ysum /   natoms
        zavg = zsum /   natoms
        zavg = 0
    return (natoms, xavg, yavg, zavg)


def test():
    s = read_fasta_file("test.fasta")
    print(s)
    print(read_pdb_file("./Glycophorin.pdb"))

if __name__ == "__main__":
    test()
