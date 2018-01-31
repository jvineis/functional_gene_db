#!/usr/bin/env python 

import sys
from Bio import SeqIO
from operator import itemgetter


addlyrs = open(sys.argv[1], 'rU')
nodes = open(sys.argv[2], 'rU')
nodelyrs = open("additional_layers_w_nodes.txt", 'w')

addlyrs_dic = {}
with addlyrs as f:
    firstline = f.readline()
    for line in f:
        x = line.strip().split('\t')
        addlyrs_dic[x[0]] = x[0:len(x)]


nodes_dic = {}

for line in nodes:
    x = line.strip().split(',')
    if round(float(x[2]),1) > float(sys.argv[3]):
        nodes_dic[x[1]] = x

nodelyrs.write("split"+'\t'+"len"+'\t'+"size"+'\t'+"pidnet"+'\t'+"alnlen"+'\t'+"mism"+'\t'+"gapopen"+'\t'+"qstart"+'\t'+"qend"+'\t'+"sstart"+'\t'+"send"+'\t'+"evalue"+'\n')

for key in addlyrs_dic.keys():
    if key in nodes_dic.keys():
        print key
        nodelyrs.write(key+'\t'+addlyrs_dic[key][1]+'\t'+addlyrs_dic[key][2]+'\t'+nodes_dic[key][2]+'\t'+nodes_dic[key][3]+'\t'+nodes_dic[key][4]+'\t'+nodes_dic[key][5]+'\t'+nodes_dic[key][6]+'\t'+nodes_dic[key][7]+'\t'+nodes_dic[key][8]+'\t'+nodes_dic[key][9]+'\t'+nodes_dic[key][10]+'\n')
    else:
        nodelyrs.write(key+'\t'+addlyrs_dic[key][1]+'\t'+addlyrs_dic[key][2]+'\t'.join(["0"]*10)+'\n')
