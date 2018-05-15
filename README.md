# Creation of a database for important functional genes using hmm models.  Currently, I have run the following hmm models - which required a lot of work and testing with models from pfam and other sources.  I found the following models to work best.  The anvio hmm directories are included in this repository.

dsrB - ftp://ftp.jcvi.org/pub/data/TIGRFAMs/

nirS - http://fungene.cme.msu.edu/hmm_detail.spr?hmm_id=21

nrfA - http://pfam.xfam.org/family/Cytochrom_C552


### Steps to create a fasta file of functional genes from 8344 ncbi "Complete" genomes explained in get_genomes_and_builddb.txt and 7904 genomes derived from https://www.nature.com/articles/s41564-017-0012-7 explained below. note-to-selfAt the time of the creation of this README I housed the scripts etc for this git here /Users/joevineis/Documents/BOWEN/functional_gene_db and running all searches and contigsdb creation from here /workspace/jvineis/functional-gene-project


1.I downloaded the genbank file ids from the parks publication. included in the repository [41564_2017_12_MOESM3_ESM.xls](https://github.com/jvineis/functional_gene_db/blob/master/41564_2017_12_MOESM3_ESM.xls)

2.I made a list of the first column only and created two files from this file.  One contains the name of the NCBI file with the trailing zeros removed, the other contains the ftp addresses of each of the files.  Here is how i did it.  All necessary files are int the repository.

    sed 's/0000000//g' genome_ids_list.txt > fix
    mv fix genome_ids_list.txt  

    python create_file_names.py genome_ids_list.txt complete-geneome-list-paths.txt

3.Then I made a bunch of temp direcories using a file with a single column of numbers 1:40 called nums

    for i in `cat nums`; do mkdir 'temp'$i; done

4.Now I went through thigs manually to create subsets of file names and paths that I want to process in each of the temp sub directories.  Something like this.

    head -n 200 genome_ids_list.txt > temp1/genome_ids_list.txt
    sed '1:200d' genome_ids_list.txt 
    head -n 200 genome_ids_list.txt > temp2/genome_ids_list.txt
    .... etc

    head -n 200 complete-genome-list-paths.txt > temp1/cgl.txt
    sed '1:200d' complete-genome-list-paths.txt
    head -n 200 complete-genome-list-paths.txt > temp2/cgl.txt
    .... etc 


5.Then I copied the find_functional_genes.shx bash script into each of those directories.  also found in this repository.

    for i in `cat nums`; do cp find_functional_genes.shx 'temp'$i/ ; done

6.Then I ran each of the bash scripts in a separate screen session.

### The steps above can be much easier to complete using a cluster configuration - which is a luxury that not everyone can afford.  In this case, I would suggest that all genomes remain in a single directory.  

7.Now concatenate all of the genes into a single fasta file for each of the single copy genes

    cat temp*/*dsrB.fa > dsrB_parks.fa
    cat temp*/*nirS.fa > nirs_parks.fa
    cat temp*/*nrfA.fa > nrfA_parks.fa

8.I concatenated the ncbi and parks data into a single file for each of the functional genes. eg

    cat dsrB_parks.fa dsrB_ncbi.fa > dsrB_ncbi_parks.fa

9.You can do the same for anmino acid sequences by altering the find_functional_genes.shx script to collect amino acid sequeneces instead.  I have a bunch of bash scripts that are designed to collect specific items from the database.  some of which are on this main page.

### Create a Taxonomy string to link with each of your sequences: Note the script below creates a single file with genome ids for the reference database and taxonomic string.  I had to run this script twice to recover the parks and ncbi data becuse the connection to NCBI is unstable.  I merged each of the taxonomy files and called it "ncbi-parks-LINEAGE-STRINGS.txt".  This is the file that you will use when linking taxonomy to your amplicon data from the sequence hit in your functional gene database

    python extract-taxonomy-strings.py -ncbi archaea-bacterial-complete-list.txt -parks parks-complete-list.txt -out parks-nrfA.tax

## functional gene tree building, phylogeny, and linking ASVs to the full tree.  These steps help us to identify where our amplicons fall in our collection of reference genes and inspect the full diversity of reference sequences.  Its best to run them in a single bash script that looks like this

    #!/bin/bash
    # all the standard stuff
    python ~/scripts/ncbi_parks_name_fix.py amoA_A-master.fa amoA_A-master_names.faa
    python /Users/joevineis/Documents/BOWEN/functional_gene_db/mu-sequence-length-selector.py -l 3000 -s 300 -fa amoA_A-master_names.faa -o amoA_ncbi_parks-sized.faa
    python /Users/joevineis/Documents/BOWEN/functional_gene_db/mu-remove-seqs-with_Ns.py -f amoA_ncbi_parks-sized.faa -o amoA_ncbi_parks-sized-Ns.faa
    python /Users/joevineis/Documents/BOWEN/functional_gene_db/mu-dereplicate-seqs.py -fa amoA_ncbi_parks-sized-Ns.faa -o amoA_ncbi_parks-sized-Ns-derep.faa
    python /Users/joevineis/Documents/BOWEN/functional_gene_db/gen-add-layers-size-len.py amoA_ncbi_parks-sized-Ns-derep.faa amoA_ncbi_parks-sized-Ns-derep-add.faa
    famsa amoA_ncbi_parks-sized-Ns-derep-add.faa amoA_ncbi_parks-sized-Ns-derep-add-famsa.faa
    trimal -in amoA_ncbi_parks-sized-Ns-derep-add-famsa.faa -out amoA_ncbi_parks-sized-Ns-derep-add-famsa-trimal.faa -gappyout
    FastTree amoA_ncbi_parks-sized-Ns-derep-add-famsa-trimal.faa > amoA_ncbi_parks-sized-Ns-derep-add-famsa-trimal.tre

    #The stuff for the ribosomal tree
    python /Users/joevineis/Documents/BOWEN/functional_gene_db/find-ribosomal-hits.py -phylo /Users/joevineis/Documents/BOWEN/functional_gene_db/ncbi-parks-ribosomal.faa -fa amoA_ncbi_parks-sized-Ns-derep-add.faa -out ribosomal-amoA-sequences-RAW.faa
    python /Users/joevineis/Documents/BOWEN/functional_gene_db/mu-sequence-length-selector.py -fa ribosomal-amoA-sequences-RAW.faa -l 5000 -s 1600 -o ribosomal-amoA-sequences-1600min.faa
    famsa ribosomal-amoA-sequences-1600min.faa ribosomal-amoA-sequences-1600min-famsa.faa
    trimal -in ribosomal-amoA-sequences-1600min-famsa.faa -out ribosomal-amoA-sequences-1600min-famsa-trimal.faa -gappyout
    FastTree ribosomal-amoA-sequences-1600min-famsa-trimal.faa > ribosomal-amoA-sequences-1600min-famsa-trimal.tre

    rm amoA_ncbi_parks-sized.faa amoA_ncbi_parks-sized-Ns.faa amoA_ncbi_parks-sized-Ns-derep.faa amoA_ncbi_parks-sized-Ns-derep-add.faa amoA_ncbi_parks-sized-Ns-derep-add-famsa.faa


    # run the commands below on a cluster
    #makeblastdb -in nrfA_ncbi_parks-sized-Ns-derep-add.faa --dbtype 'prot' -title nrfA-blastdb -out nrfA-blastdb
    clusterize blastp -db dsrB_tigr_blastdb -out dsrB_ncbi_parks_tigr-node-blast-hits.txt -query dsrB-pooled-samples-node-representatives.faa -outfmt "6 qseqid sseqid pident length qlen sstart send qstart qend evalue" -max_target_seqs 1

    #Continue with the script below to make your amazing layers
    #python ../add-node-hits-to-additional-layers -i additional_layers.txt -blast nrfA_ncbi_parks-node-blast-hits.txt -per 60 -tax ~/scripts/databas/functional_fa_dbs/ncbi-parks-LINEAGE-STRINGS.txt

1.Remove sequences that are too divergent in legth to be the target.  Most amino acid sequences are below 350aa.  I found this after multiple iterations of building trees and alignments.  This is very helpful 

    python ../mu-sequence-length-selector.py -l 3000 -s 300 -fa nrfA_ncbi_parks.faa -o nrfA_ncbi_parks-sized.faa

we started with 2189 sequences in the fix.faa file and ended with 2095 in the sized file.

remove amino acid sequences with "X" They contained "Ns" in their nucleotide scaffold assembly.
 
    python ../mu-remove-seqs-with_Ns.py -f nirS_ncbi_parks-sized.faa -o nirS_ncbi_parks-sized-Ns.faa

2.Dereplicate the sequences

    python ~/scripts/mu-dereplicate-seqs.py -fa nrfA_ncbi_parks-sized-Ns.faa -o nrfA_ncbi_parks-sized-Ns-derep.faa

1012 sequences remain after dereplicating.

3.We want to add some additional layers that will help when visualizing this beast with Anvio.  Lets start with the size of the fragment and the number of sequences represented in the dereplicated data.  This information can be produced from headers that look like this. ">GCA_001278275|source_start:1367327|stop:1368764;size=116;" using a script called gen-add-layers-size-len.py like below.  This will produce a faa with a headers that looks like this ">GCA_001278275_1367327_1368764" and a file called "additional_layers.txt".

     python gen-add-layers-size-len.py nrfA_ncbi_parks-sized-derep.faa nrfA_ncbi_parks-sized-derep-add.faa

4.Run famsa to align the dereplicated sequences

    famsa nrfA_ncbi_parks-sized-derep-add.faa nrfA_ncbi_parks-sized-derep-add-famsa.faa

4a.Run trimal on the alignment.

    trimal -in nrfA_ncbi_parks-sized-derep-add-famsa.faa -out nrfA_ncbi_parks-sized-derep-add-famsa-trimal.faa -gappyout

5.Run FastTree on the alignment.  This step produces our tree for visualization

    FastTree nrfA_ncbi_parks-sized-derep-add-famsa-trimal.faa > nrfA_ncbi_parks-sized-derep-add-famsa-trimal.tre

6.In order to see where the nodes from our nrfA amplicon study from the thin sections map onto our tree in an unbiased way, we used blast to search our translated NODE-REPRESENTATIVES.fa.  Here are the steps to create a layer for the node hits and details of the match for visualization in anvio.  Beyond these steps will need the "additional_layers.txt" file produced in step 3 

     prodigal -i ~/Dropbox/oxygen_gradient/nrfA/swarm-d1-fastidious/swarm-d1-fastidious-node-representatives.fa -a ~/Dropbox/oxygen_gradient/nrfA/swarm-d1-fastidious/swarm-d1-fastidious-node-representatives.faa -n -p meta
     makeblastdb -in nrfA_ncbi_parks-sized-derep-add.faa -dbtype 'prot' -title nrfA_ncbi_parks-sized-derep-add -out nrfA_ncbi_parks-sized-derep-add
     blastp -db nrfA_ncbi_parks-sized-derep-add -out nrfA_ncbi_parks-node-blast-hits.txt -outfmt 10 -query ~/Dropbox/oxygen_gradient/nrfA/nrfA-MED-OUTPUT/NODE-REPRESENTATIVES.faa


9.Now we can run the step below which will create a text file called additional_layers_w_nodes.txt which can then be visualized with anvio

    python add-node-hits-to-additional-layers.py -i additional_layers.txt -blast nrfA_ncbi_parks-node-blast-hits.txt -per 60 -tax ncbi-parks-LINEAGE-STRINGS.txt
    anvi-interactive -t nrfA_ncbi_parks-sized-derep-add-famsa-trimal.tre -p profile.db -d additional_layers_w_nodes.txt --manual-mode

## This tree may be not as phylogenetically accurate as the tRNA genes used in the [Hug et.al](https://www.nature.com/articles/nmicrobiol201648) to infer phylogeny.  Lets grab those genes, create the tree, map on our ASVs and taxonomy.  

1. Lets try some Hug tRNA genes from the 16K genomes using anvi-get-sequences-for-hmm-hits like this. 

    bash x_collect_ribosomal.shx

which looks like this.
      
    #!/bin/bash

    for i in `cat samples.txt`;
    do 
       anvi-get-sequnces-for-hmm-hits -c "$i"_contigs.db --hmm-source Campbell_et_al --gene-names ../list_of_ribosomal_proteins.txt --return-best-hit --get-aa-sequences --concatenate;
    done

here is the list of ribosomal proteins used in the script above

    Ribosomal_L2
    Ribosomal_L3
    Ribosomal_L4
    Ribosomal_L5
    Ribosomal_L6
    Ribosomal_L10
    Ribosomal_L11
    Ribosomal_L12
    Ribosomal_L13
    Ribosomal_L14
    Ribosomal_L16
    Ribosomal_L22
    Ribosomal_S8
    Ribosomal_S10
    Ribosomal_S17
    Ribosomal_S19

I concatenated all of the ribosomal proteins 

    cat temp*/*ribosomal.fa > all_ncbi_ribosomal.fa
    cat temp*/*ribosomal.fa > all_parks_ribosomal.fa

To collect the ribosomal sequences associated with a particular functional gene, I used this magical script

    python find-ribosomal-hits.py -phylo ncbi-parks-ribosomal.faa -fa dsrB_ncbi_parks-nonredundant-blastfilter-sized-Ns-derep-add.faa -out ribosomal-dsrB-sequences-RAW.faa

Then I selected sequences that were greater than 1600bp.

    python ~/scripts/mu-sequence-length-selector.py -fa ribosomal-nrfA-sequences-RAW.faa -l 5000 -s 1600 -o ribosomal-nrfA-sequences-1600min.faa

You can now continue on with this output at step 3 of "nrfA Tree building" to add additional layers etc..

# I'm currently using SWARM v2.2.2 to cluster ASVs.  Here is how I use the nrfA database to assign taxonomy to each node and set up for Phyloseq visualization etc.  
# Build a Phyloseq Object using the output nodes from MED or SWARM etc and vsearch taxonomic assignment - This is amplicon stuff that is somewhat misplaced here, but I'm leaving it for now.  Here is and example of how to run things for SWARM.  MED and DADA2 can be run in a similar way but will need a little development.  Let me know if you desire this.

1.Run swarm starting with a concatenated fasta file of all of your samples
    vsearch --derep_fulllength pooled-samples.fa --sizeout --output pooled-samples-derep.fa
    swarm -d 1 -f -t 10 -z pooled-samples-derep.fa -s pooled-samples-derep-stats.txt -w pooled-samples-node-representatives.fa -o pooled-samples-node-table.txt

2.Run a bunch of scripts to make useful tables,phyloseq objects, and anvio trees for visualizing relevant metadata.  I have this set up so that swarms with a minimum abundace of 0,10,3 are all created and run.  You can run one or all of them. Just bust the scripts out of the loop if you are not interested in multiple abundance thresholds.

    #!/bin/bash

    for i in 10 3 0
    do
        mkdir min-count"$i"
        python ~/Documents/swarm-analysis/mu-swarms-to-ASVs-table.py -s pooled-samples-node-table.txt -o reduced-min"$i"-matrix-count.txt -l ../samples.txt -n pooled-samples-node-representatives.fa -min "$i"
        vsearch --usearch_global reduced-node-fasta-min"$i".fa --db ~/Documents/BOWEN/functional_gene_db/dsrB/pfam_NIR_SIR/dsrB_ncbi_parks.fa --blast6out NODES-min"$i"-HITS.txt --id 0.6
        python ~/Documents/swarm-analysis/mu-swarm-to-phyloseq-objects.py -tax_ref ~/scripts/databas/functional_fa_dbs/ncbi-parks-LINEAGE-STRINGS.txt -hits NODES-min"$i"-HITS.txt -med reduced-min"$i"-matrix-count.txt -fa reduced-node-fasta-min"$i".fa
        mv PHYLOSEQ-MATRIX-OBJECT.txt PHYLOSEQ-TAX-OBJECT.txt reduced-node-fasta-min"$i".fa reduced-min"$i"-matrix-count.txt min-count"$i"
        famsa min-count"$i"/reduced-node-fasta-min"$i".fa min-count"i"/reduced-node-fasta"$i"-famsa.fa
        trimal -in min-count"$i"/reduced-node-fasta"$i"-famsa.fa -out min-count"i"/reduced-node-fasta"$i"-famsa-trimal.fa -gappyout
        FastTree -nt min-count"$i"/reduced-node-fasta"$i"-famsa-trimal.fa > min-count"i"/reduced-node-fasta"$i"-famsa-trimal.tre
    done
 


3.Now you have everything that you need to build a Phyloseq object using R.  Make sure that you have R version 3.4, bioconductor 3.6 and phyloseq 1.22.  There are bugs in the earlier version of phyloseq.  You don't need to have any metadata but it sure does help. The sample names in the metadata file must be exactly the same as those contained in the PHYLOSEQ-MATRIX-OBJECT and you will need a header for each column.  

    library("phyloseq")
    library("ape")
    mat_nirS = read.table("PHYLOSEQ-MATRIX-OBJECT.txt", header = TRUE, sep = "\t", row.names = 1)
    tax_nirS = read.table("PHYLOSEQ-TAX-OBJECT.txt", header = TRUE, sep = ";", row.names = 1)
    meta_nirS = read.table("nrfA-map.txt", header = TRUE, sep = "\t", row.names = 1)
    tree_nirS = read.tree("reduced-node-famsa-trimal.tre")    

    mat_nirS = as.matrix(mat_nirS)
    tax_nirS = as.matrix(tax_nirS)

    OTU = otu_table(mat_nirS, taxa_are_rows = TRUE)
    TAX = tax_table(tax_nirS)
    META = sample_data(meta_nirS)

    nrfA_physeq = phyloseq(OTU,TAX,META,tree_nirS)

## Update - I am working toward trying to cluster amino acid sequences for the short reads using MED.  This involves the dangerous practice of translating the short reads using prodigal.. 



1.Start out by using prodigal to call amino acid sequences for all short reads - this takes a while.

	prodigal -i sequences-to-decompose.fa -a sequences-to-decompose.faa -n -p meta



Each of these collections can now be used as a db to search using vsearch or converted into a blast databaase etc. in the same way that you us\
e the silva database for your 16S data.  Enjoy!  If you have another maker that you would like to explore, just let me know and I'll make it h\
appen :)