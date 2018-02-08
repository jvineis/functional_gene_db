#!/usr/bin/env python 

import argparse
from Bio import SeqIO
from Bio import Entrez
Entrez.email = "jvineis@gmail.com"

parser = argparse.ArgumentParser(description='''This takes the genbank genome id and returns the entire lineage string''')
parser.add_argument('-ncbi', help = 'a single column list of genbank ids to search for tax strings')
parser.add_argument('-parks', help = 'a two column list of genomeids and NCBI UID numbers')
parser.add_argument('-out', help = 'an outfile with columns containg the id of the genome as the first column and the hit in the following columns')
args = parser.parse_args()

out = open(args.out,'w')

for line in open(args.parks, 'rU'):
    Entrez.email = "jvineis@gmail.com"
    x = line.strip().split('\t')
    handle = Entrez.efetch(db="nucleotide", id=x[2], retmode = "xml")
    record = Entrez.read(handle)
    ts = record[0]["GBSeq_taxonomy"].split(';')
    print("parks", len(ts), ts)
    
    if len(ts) == 0 :
        out.write(x[2][0]+x[2][1]+x[2][2]+x[2][3]+x[2][4]+'\t'+"na"+'\t'+"na"+'\t'+"na"+'\t'+"na"+'\t'+"na"+'\t'+"na"+'\t'+record[0]["GBSeq_organism"]+'\t'+record[0]["GBSeq_source"]+'\n')
    elif len(ts) == 1:
        out.write(x[2][0]+x[2][1]+x[2][2]+x[2][3]+x[2][4]+'\t'+ts[0]+'\t'+"na"+'\t'+"na"+'\t'+"na"+'\t'+"na"+'\t'+"na"+'\t'+record[0]["GBSeq_organism"]+'\t'+record[0]["GBSeq_source"]+'\n')
    elif len(ts) == 2:
        out.write(x[2][0]+x[2][1]+x[2][2]+x[2][3]+x[2][4]+'\t'+ts[0]+'\t'+ts[1]+'\t'+"na"+'\t'+"na"+'\t'+"na"+'\t'+"na"+'\t'+record[0]["GBSeq_organism"]+'\t'+record[0]["GBSeq_source"]+'\n')
    elif len(ts) == 3:
        out.write(x[2][0]+x[2][1]+x[2][2]+x[2][3]+x[2][4]+'\t'+ts[0]+'\t'+ts[1]+'\t'+ts[2]+'\t'+"na"+'\t'+"na"+'\t'+"na"+'\t'+record[0]["GBSeq_organism"]+'\t'+record[0]["GBSeq_source"]+'\n')
    elif len(ts) == 4:
        out.write(x[2][0]+x[2][1]+x[2][2]+x[2][3]+x[2][4]+'\t'+ts[0]+'\t'+ts[1]+'\t'+ts[2]+'\t'+ts[3]+'\t'+"na"+'\t'+"na"+'\t'+record[0]["GBSeq_organism"]+'\t'+record[0]["GBSeq_source"]+'\n')
    elif len(ts) == 5:
        out.write(x[2][0]+x[2][1]+x[2][2]+x[2][3]+x[2][4]+'\t'+ts[0]+'\t'+ts[1]+'\t'+ts[2]+'\t'+ts[3]+'\t'+ts[4]+'\t'+"na"+'\t'+record[0]["GBSeq_organism"]+'\t'+record[0]["GBSeq_source"]+'\n')
    elif len(ts) >= 6:
        out.write(x[2][0]+x[2][1]+x[2][2]+x[2][3]+x[2][4]+'\t'+ts[0]+'\t'+ts[1]+'\t'+ts[2]+'\t'+ts[3]+'\t'+ts[4]+'\t'+ts[5]+'\t'+record[0]["GBSeq_organism"]+'\t'+record[0]["GBSeq_source"]+'\n')

#for line in open(args.ncbi, 'rU'):
#    Entrez.email = "jvineis@gmail.com"
#    x = line.strip().split('\t')
#    handle = Entrez.efetch(db="taxonomy", id=x[6], rettype="xml", retmode = "text")
#    record = Entrez.read(handle)
#    for taxon in record:
#        ts = taxon["Lineage"].split(";")
#        print("ncbi", ts)
#        if len(ts) < 7:
#            print(x[0], "is not a complete record")
#            out.write(str(x[0].split('.')[0])+'\t'+"na"+ '\t' +"na"+ '\t' +"na"+ '\t' +"na"+ '\t' +"na"+ '\t' +"na"+ '\t' +"na"+ '\t' +"na"+'\n')
#        else:#            out.write(str(x[0].split('.')[0]) + '\t' + ts[1]+ '\t' + ts[2] + '\t' + ts[3]+ '\t' +ts[4]+ '\t' +ts[5]+ '\t' +ts[6]+ '\t'+x[7]+'\t'+"na"+'\n')
