#!/usr/bin/env python
import sys
from Bio import SeqIO

fa = SeqIO.parse(open(sys.argv[1], 'rU'), "fasta")
hits = open(sys.argv[2], 'rU')
out = open(sys.argv[3], 'w')

names = {}
for hit in hits:
    x = hit.strip().split("|")
    if x[0][0] == "G":
        part1 = x[0].split("_")[0]+"_"+x[0].split("_")[1]
        part2 = x[1]+"|"+x[2]
        hit_name = part1+"|"+part2
        print hit_name
        names[hit_name] = x
    elif x[0][0] == "D":
        part1 = x[0].split("_")[0]
        part2 = x[1]+"|"+x[2]
        hit_name = part1+"|"+part2
        names[hit_name] = x
    
for rec in fa:
    if rec.id in names.keys():
        out.write(">"+str(rec.id)+"\n"+str(rec.seq)+"\n")
