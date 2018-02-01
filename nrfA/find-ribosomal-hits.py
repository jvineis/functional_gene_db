#!/usr/bin/env python

import argparse
from Bio import SeqIO
parser = argparse.ArgumentParser(description='''Select sequences from the ncbi-parks-ribosomal.fa that are found in the fasta file of functional gene hits (hmm) to the same genome''')
parser.add_argument('-phylo' , help = 'The single copy concatenated genes fasta for all 16K reference genomes (if they had them)')
parser.add_argument('-fa' , help = 'The fasta file of functional genes. We only want sequeuces from "phylo" with headers in this file')
parser.add_argument('-out', help = 'The output fasta file')
args = parser.parse_args()

out = open(args.out,'w')
fa_names = {}
fa = SeqIO.parse(open(args.fa, 'rU'), 'fasta')

for rec in fa:
    x = rec.id.split('_')[0]
    y = rec.id.split('_')[1]
    x_y = x+"_"+y+"_contigs"
    fa_names[x_y] = rec.id

for rec in SeqIO.parse(open(args.phylo, 'rU'), 'fasta'):
    if rec.id in fa_names.keys():
        out.write(">"+fa_names[rec.id]+'\n'+str(rec.seq)+'\n')

