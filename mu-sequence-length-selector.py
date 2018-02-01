#!/usr/bin/env python
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='''Define a lenght and a fasta file to recieve only sequences of the dength you dictate.'\n'mu-sequence-length-selector.py -l "length" -fa "fasta.fa" -o "file to hold sequences that are just the right length''')
parser.add_argument('-l' , help = ' a numerical value for the longest sequence length')
parser.add_argument('-s', help = ' a numerical value for the shortest acceptable length')
parser.add_argument('-fa' , help = 'a fasta file')
parser.add_argument('-o' , help = 'the name of the fasta output file')
args = parser.parse_args()

out = open(args.o,'w')

for rec in SeqIO.parse(open(args.fa,'rU'), "fasta"):
    if int(len(rec.seq)) <= int(args.l) and int(len(rec.seq)) >= int(args.s): 
        print rec.seq
        out.write(">" + str(rec.id) + '\n' + str(rec.seq) + '\n')
