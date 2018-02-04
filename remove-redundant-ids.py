#!/usr/bin/env python
import sys
from Bio import SeqIO

outfile = open(sys.argv[2], 'w')
ids = []
for rec in SeqIO.parse(open(sys.argv[1], 'rU'), "fasta"):
    if rec.id in ids:
        next
        print("I found this one more than onece",rec.id)
    else:
        outfile.write(">"+str(rec.id)+"\n"+str(rec.seq)+"\n")
        ids.append(rec.id)

