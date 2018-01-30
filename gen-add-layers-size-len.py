#!/usr/bin/env python

from Bio import SeqIO
import sys

outfile = open(sys.argv[2], 'w')
addlyr = open("additional_layers.txt", 'w')
addlyr.write("split"+'\t'+"seq_length"+'\t'+"count"+'\n')
for rec in SeqIO.parse(open(sys.argv[1], 'rU'), "fasta"):
    id =  rec.id.split("|")[0]
    start = rec.id.split("|")[1].split(":")[1]
    stop = rec.id.split("|")[2].split(";")[0].split(":")[1]
    cnt =  rec.id.split(";")[1].split("=")[1]
    outfile.write(">"+id+'_'+start+'_'+stop+'\n'+str(rec.seq)+'\n')
    addlyr.write(id+'_'+start+'_'+stop+'\t'+str(len(rec.seq))+'\t'+cnt+'\n')
    
outfile.close()
addlyr.close()
