import Bio
from Bio.Seq import Seq
from Bio import motifs

instances = []
instances.append(Seq("TATAA"))
instances.append(Seq("TATTA"))
instances.append(Seq("TTTAT"))
instances.append(Seq("TATAC"))

m = motifs.create(instances)

print(type(m))
print(m)
print(len(m))
print(m.consensus)
print(m.pwm)
print(m.counts)
print(m.pssm)

# m.weblogo("motif.png")
