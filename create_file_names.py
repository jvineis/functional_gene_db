#!/usr/bin/env python
import sys

file = open(sys.argv[1],'rU') 
outfile = open(sys.argv[2],'w')
for line in file:
    a = line.strip()
    b = line[0:2]
    c = line[2:4]
    d = line[0:5]
    print(a,b,c,d)
    outfile.write('ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux'+'/'+b+'/'+c+'/'+d+'1'+'/'+d+'1.1.fsa_nt.gz'+'\n')

