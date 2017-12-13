#!/usr/bin/env python

import sys

# A file with a single column of deflines with the ">" removed
deflines = open(sys.argv[1], 'rU')
# A list of the ncbi genomes and assigned taxonomy.
ncbi_archive = open(sys.argv[2], 'rU')
# The name of the output file where you would like to write the functional_gene_ids.txt file.  I suggest 
# you call it something like nrfA_ids.txt
out = open(sys.argv[3], 'w')

defline = {}
for line in deflines:
    x = line.strip().split("|")
    print x[0]
    defline[x[0]] = [x[0], x[1], x[2]]


archive_dict = {}
for line in ncbi_archive:
    x =  line.strip().split("\t")
    if "UBA" in x[4]:
        archive_dict[x[0]] = [x[2],x[5]]
    else:
        archive_dict[x[0]] = [x[0],x[7]]


#for key in archive_dict.keys():
#    print archive_dict[key]

for i in defline.keys():
    for key in archive_dict.keys():
#        print i[0] , archive_dict[key][0]
        if defline[i][0] in archive_dict[key][0]:
            out.write(defline[i][0] +'|'+defline[i][1]+'|'+defline[i][2]+'\t'+archive_dict[key][1]+'\n')

            
