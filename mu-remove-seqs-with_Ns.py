#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
import argparse

parser = argparse.ArgumentParser(description='''Takes a fasta file and removes any sequence containing an "N"''')
parser.add_argument('-f', help = 'The input fasta file')
parser.add_argument('-o', help = 'The output fasta file')
args = parser.parse_args()

fa = open(args.f, 'rU')
o = open(args.o, 'w')
n = 0
c = 0 

fasta = SeqIO.parse(fa, "fasta")
for rec in fasta:
    if "N" in rec.seq:
        n += 1
    else:
        SeqIO.write(rec, o, "fasta" )
        c += 1
              

print("Nice work microbial scientist\nfound %d sequences with N's and %d clean sequences" %(n,c))
