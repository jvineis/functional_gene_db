# Creation of a database for important functional genes
### Steps to create a fasta file of functional genes from 6000 ncbi "Gold" genomes explained in get_genomes_and_builddb.txt and 8000 genomes derived from https://www.nature.com/articles/s41564-017-0012-7 explained below. At the time of the creation of this README, I was working from here /Users/joevineis/Documents/BOWEN/functional_gene_db.

1.I downloaded the genbank file ids from the parks publication. included in the repository [41564_2017_12_MOESM3_ESM.xls](https://github.com/jvineis/functional_gene_db/blob/master/41564_2017_12_MOESM3_ESM.xls)

2.I made a list of the first column only and created two files from this file.  One contains the name of the NCBI file with the trailing zeros removed, the other contains the ftp addresses of each of the files.  Here is how i did it.  All necessary files are int he repository

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

6.Then I ran each of the bash scripts in a separate screen session.  I know that this is laborious and tiresome but its the only way I know of how to do it that maximizes speed and efficiency.  

7.Now concatenate all of the genes into a single fasta file for each of the single copy genes

    cat temp*/*dsrB.fa > dsrB_parks.fa
    cat temp*/*nirS.fa > nirs_parks.fa
    cat temp*/*nrfA.fa > nrfA_parks.fa

8.I concatenated the ncbi and parks data into a single file for each of the functional genes. eg

    cat dsrB_parks.fa dsrB_ncbi.fa > dsrB_ncbi_parks.fa
    sed 's/_contigs//g' dsrB_ncbi_parks.fa | sed 's/ //g' > dsrB

9.Fix the deflines that are produced by anvio using the handy ncbi_parks_name_fix.py.  eg

    python ncbi_parks_name_fix.py dsrB dsrB_ncbi_parks.fa

If you are comfortable with the results, you can remove all the temporary stuff.       

### DONT DO THIS ANYMORE
10.The taxa for each of the genomes are contained in the NCBI archive tables for archeal, bacterial, and parks genomes are found in the tables called archaea-bacterial-complete-list.txt and parks-complete-list.txt. So I guess its time to create a taxonomy and defline table that looks something like this.

    Defline_from_anvio    taxonomy_from_archive_table 
    GCA_000005825|source_start:73135|stop:74284    Methanococcus voltae

   begin by making a list of the deflines 
    
    grep ">" dsrB_ncbi_parks.fa | sed 's/>//g' > dsrB_deflines.txt

   add the taxonomy to the defline using the add_taxa_to_defline.py.  I created the ncbi-parks-genome-2-taxa-list.txt simply by concatenating the lists of each of the genomes.  There are unequal numbers of columns in this dataset.  The ncbi "gold" genomes information contains more values than the parks genome id table.  
    
    add_taxa_to_defline.py nirS/nirS_deflines.txt ncbi-parks-genome-2-taxa-list.txt nirS_ids.tx

11.I ran taxize on the genus level taxonomy of the ncbi_parks_tax_genera.txt, where "d" is the comma separated list of the genus names.
     
    d =  c("Bacteroides","Colostidium","Alphaproteobacteria".....)
    work = tax_name(d, get = c("phylum","class","order","family","genus"), db = "ncbi")
    write.table(work, "~/scripts/databas/functional_fa_dbs/raw_taxit.txt")


12.I added the taxonomic string found in raw_taxit.txt produced by the R code above to the appropriate defline using mu-convert-ncbi-parks-ids-to-taxastring.py. Done like this for each functional gene.    

    python ~/scripts/mu-convert-ncbi-parks-ids-to-taxastring.py nirS_ids.txt raw_taxit.txt nirS.tax

if you see an error here its possible that you need to remove the first line (header information) of your raw_taxit.txt file

Each of these collections can now be used as a db to search using vsearch or converted into a blast databaase etc. in the same way that you use the silva database for your 16S data.  Enjoy!  If you have another maker that you would like to explore, just let me know and I'll make it happen :) 

## Update - I have found there are some questionable sequences containing Ns.  I used mu-remove-seqs-with_Ns.py to remove these sequences from the database. I'm continuing to work toward refinement and improvement of these databases. 

## Here are some improvements to nrfA that I have made.  

### Build a tree based on nrfA to expole the possibility of non-nrfA sequences.  
1.Remove sequences that are too the wrong size to be the target.  Most amino acid sequences are below 350aa.  I found this after multiple iterations of building trees and alignments.  This is very helpfu.

    python ../mu-sequence-length-selector.py -l 650 -s 200 -fa nrfA_ncbi_parks.faa -o nrfA_ncbi_parks-sized.faa

we started with 2189 sequences in the fix.faa file and ended with 1825 in the sized file. Add the Thioalialivibrio nitratireducens from [here](https://www.ncbi.nlm.nih.gov/protein/2OT4_B?report=fasta) to the fasta file.  This is based on recommendation by Walsh et.al. !!!! NOT ADDING THE T.nitratireducens.. it is not different enough from nrfA given our approach

2.Dereplicate the sequences

    vsearch --derep_fulllength nrfA_ncbi_parks-sized.faa --output nrfA_ncbi_parks-sized-derep.faa --sizeout

1012 sequences remain after dereplicating.

3.We want to add some additional layers that will help when visualizing this beast with Anvio.  Lets start with the size of the fragment and the number of sequences represented in the dereplicated data.  This information can be produced from headers that look like this. ">GCA_001278275|source_start:1367327|stop:1368764;size=116;" using a script called gen-add-layers-size-len.py like below.  This will produce a "nrfA_ncbi_parks-derep-fix.faa" with a headers that looks like this ">GCA_001278275_1367327_1368764" and a file called "additional_layers.txt".

     python gen-add-layers-size-len.py nrfA_ncbi_parks-sized-derep.faa nrfA_ncbi_parks-sized-derep-add.faa

4.Run famsa to align the dereplicated sequences

    famsa nrfA_ncbi_parks-sized-derep-add.faa nrfA_ncbi_parks-sized-derep-add-famsa.faa

4a.Run trimal on the alignment.

    trimal -in nrfA_ncbi_parks-sized-derep-add-famsa.faa -out nrfA_ncbi_parks-sized-derep-add-famsa-trimal.faa -gappyout

5.Run FastTree on the alignment

    FastTree nrfA_ncbi_parks-sized-derep-add-famsa-trimal.faa > nrfA_ncbi_parks-sized-derep-add-famsa-trimal.tre

6a.After inspecting the tree its obvious that there are some non-nrfA seqs that cluster with our Thioalkalivibrio.  I selected all of the branches in the tree and then clicked on the "split" box to see the genome id.  I copied and pasted this into a file called non_nrfA_ids.txt.  The used mu-selectseq_from_fasta.py like this.. ### After furthrther inspection.. I'm ok with the sequences clustered around the Thioalkalivibrio because they are not divergent enough given the amount of diversity that we see with nrfA.. However, we can play with this if we like and return to step 1 with the reduced sequence list.

    python ~/scripts/mu-selectseq_from_fasta.py --x nrfA_ncbi_parks.fix-sized-derep.faa --out nrfA_ncbi_parks.fix-sized-derep-clean.faa --infile non_nrfA_ids.txt --list TRUE

then I reran all of the famsa, trimal, and FastTree to create my clean tree.  Hope this looks nice!

7.In order to see where the nodes from our nrfA amplicon study from the thin sections map onto our tree in an unbiased way, we used blast to search our translated NODE-REPRESENTATIVES.fa.  Here are the steps to create a layer for the node hits and details of the match for visualization in anvio.  These steps will need the "additional_layers.txt" file produced in step 6.  

     prodigal -i ~/Dropbox/oxygen_gradient/nrfA/nrfA-MED-OUTPUT/NODE-REPRESENTATIVES.fa -a ~/Dropbox/oxygen_gradient/nrfA/nrfA-MED-OUTPUT/NODE-REPRESENTATIVES.faa -n -p meta
     makeblastdb -in nrfA_ncbi_parks-sized-derep-add.faa -dbtype 'prot' -title nrfA_ncbi_parks-sized-derep-add -out nrfA_ncbi_parks-sized-derep-add
     blastp -db nrfA_ncbi_parks-sized-derep-add -out nrfA_ncbi_parks-node-blast-hits.txt -outfmt 10 -query ~/Dropbox/oxygen_gradient/nrfA/nrfA-MED-OUTPUT/NODE-REPRESENTATIVES.faa
     python add-node-hits-to-additional-layers.py -i additional_layers.txt -blast nrfA_ncbi_parks-node-blast-hits.txt -per 60
     anvi-interactive -t nrfA_ncbi_parks.fix-derep-famsa.tre -d additional_layers_w_nodes.txt -p profile.db --manual-mode

8.We are getting a little greedy now and would like to add the taxonomy to the tree.  Here is how I'm now using Entrez in biopython to collect the taxonomy information from all records in the 16,000 genomes, not just nrfA

9.If I ever wanted to collect ribosomal RNA genes from the 16K genomes, I could use the script below.  But there is not enough congurencey among genomes to do anything useful with this information so this step 9 is a dead end for now.

     bash x_collect_rnas.shx

which looks like this.

    #!/bin/bash
    module load python/3.6.3-201710101533
    for i in `cat samples.txt`
    do
        anvi-get-sequences-for-hmm-hits -c "$i"_contigs.db --hmm-sources Ribosomal_RNAs -o "$i"_rRNA_output.fa
    done 

then I concatenate all of the ribosomal RNAs from the ncbi and parks db directories in my home directory on the MBL servers.

     cat temp*/*rRNA_output.fa > all_ncbi_ribosomal_RNAs.fa
     cat temp*/*rRNA_output.fa > all_parks_ribosomal_RNAs.fa

now select only 16S sequences (ignoring 23S for now).  Well, I would do this, but only 797 out of the 8000 parks genomes contained 16S sequences.  Can this be right?  I need to check the Parks paper again.  

10.Lets try some Hug tRNA genes from the 16K genomes using anvi-get-sequences-for-hmm-hits like this.

    bash x_collect_ribosomal.shx

which looks like this.
      
    #!/bin/bash
    module load python/3.6.3-201710101533
    for i in `cat samples.txt`;
    do 
       anvi-get-sequnces-for-hmm-hits -c "$i"_contigs.db --hmm-source Campbell_et_al --gene-names ../list_of_ribosomal_proteins.txt --return-best-hit --get-aa-sequences --concatenate;
    done

I concatenated all of the ribosomal proteins 

    cat temp*/*ribosomal.fa > all_ncbi_ribosomal.fa
    cat temp*/*ribosomal.fa > all_parks_ribosomal.fa

To collect the ribosomal sequences associated with a particular functional gene, I used this magical script

    python find-ribosomal-hits.py -phylo ncbi-parks-ribosomal.faa -fa nrfA_ncbi_parks.fix-derep.faa -out ribosomal-nrfA-sequences-RAW.fa

Then I selected sequences that were greater than 1600bp. ran famsa and FastTree

    python ~/scripts/mu-sequence-length-selector.py -fa ribosomal-nrfA-sequences-RAW.fa -l 5000 -s 1600 -o ribosomal-nrfA-sequences-1600min.fa
    famsa ribosomal-nrfA-sequences-1600min.fa ribosomal-nrfA-sequences-1600min-famsa.faa
    trimal -in ribosomal-nrfA-sequences-1600min-famsa.faa -out ribosomal-nrfA-sequences-1600min-famsa-trimal.faa -gappyout
    FastTree ribosomal-nrfA-sequences-1600-famsa.faa > ribosomal-nrfA-sequences-1600-famsa.tre

##Here are some fixes that I made to dsrB.  This one was pretty phucked up due to the presence of multiple gene families in the pfam.

Just in case you forgot. Here is how to run a blast of the protein sequences against the most recent protein ref contained at MBL. I run it on the cluster like thus

    bash x_extract-aa-and-submit-blast.shx

which looks like this

    #!/bin/bash
    module load python/3.6.3-201710101533
    for i in `cat samples.txt`
    do
        anvi-get-sequences-for-hmm-hits -c "$i"_contigs.db --hmm-sources dsrB --get-aa-sequences -o "$i".dsrB.faa
        sed 's/ /_/g' "$i".dsrB.faa > "$i".dsrB.fix.faa
        mv "$i".dsrB.fix.faa "$i".dsrB.faa
        python ../ncbi_parks_name_fix.py "$i".dsrB.faa "$i".dsrB.fix.faa
        mv "$i".dsrB.fix.faa "$i".dsrB.faa
        clusterize blastp -query "$i".dsrB.faa -db /usr/local/blastdb/refseq_protein.26 -outfmt \"6 qseqid sseqid length mismatch evalue qstart qend sstart send salltitles\" -out "$i".dsrB_blast.txt
    done

1.collect all of the hits from both the parks and ncbi dirctories and then make a single concatenated list of hits.  Then pull out just the node ids that have the "sulfite reductase" in the hit id

    cat parks*/temp*/*dsrB_blast.txt > parks-dsrb-hits.txt
    cat ncbi*/temp*/*dsrB_blast.txt > ncbi-dsrb-hits.txt
    cat parks-dsrb-hits.txt ncbi-dsrb-hits.txt > ncbi-parks-dsrB-blast-hits.txt
    grep "sulfite reductase" ncbi-parks-dsrb-allblast.txt | cut -f 1 | sort | uniq > ncbi-parks-dsrb-blast-truehits.txt
    
2.Then I use this filter to return only the sequences that are in the hit table created above.  Its crazy but sometimes anvio returns two ids that are nearly identical.  I use the script below to remove any sequnce that is redundant.  They are acutally differnt sequences but nearly identical.. I don't know exactly whats happening here but I need to talk with Meren on this one.  In 14,000 dsrB sequences, I only found one instance of this DCPD0|source_start:0|stop:231.  

    python remove-redundant-ids.py dsrB_ncbi_parks.fa dsrB_ncbi_parks-nonredundant.fa
    python filter-dsrB-using-blast.py dsrB_ncbi_parks-nonredundant.faa ncbi-parks-dsrb-blast-truehits.txt dsrB_ncbi_parks-nonredundant-blastfilter.faa

3.Now run famsa and do other awesome stuff.  But remember that famsa requires you to remove some of the characters that you like in you fasta file.  This works most of the time.

    sed 's/\|/_/g' dsrB_ncbi_parks-nonredundant-blastfilter.faa | sed 's/:/_/g' > dsrB_ncbi_parks-nonredundant-blastfilter-sed.faa

# Build a Phyloseq Object using the output nodes from MED and vsearch taxonomic assignment


1.Run MED in the way you like.  Something like this will do.

    decompose -M 20 sequences-padded.fa

2.In the output, you will find a file called NODE-REPRESENTATIVES.fasta.  We need to fix this file in order to proceed.

    sed 's/-//g' NODE-REPRESENTATIVES.fasta | sed 's/\|/_/g' | sed 's/:/_/g' > NODE-REPRESENTATIVES.fa

3.Obtain taxonomy using vsearch.  

    vsearch --usearch_global NODE-REPRESENTATIVES.fa --db nirS_ncbi_parks.fa --blast6out NODE-HITS.txt --id 0.6
	 
4.The following script will produce PHYLOSEQ-TAX-OBJECT.txt and PHYLOSEQ-MATRIX-OBJECT.txt that will be used for input into R as the otu_table and tax_table.

     mu-silva_ids_to_phyloseq.py -tax_ref nirS.tax -hits NODE-HITS.txt -med MATRIX-COUNT.txt -fa NODE-REPRESENTATIVES.fa

5.Now you need a tree file to do some amazing expoloration with phyloseq.  I use MUSCLE v3.8.31 to align my NODE-REPRESENTATIVES.fa like this.

    muscle -in NODE-REPRESENTATIVES.fa -out nirS-muscle-alignment.fa

5a.I'm now using [famsa](https://github.com/refresh-bio/FAMSA)

    e.g. famsa NODE-REPRESENTATIVES.fasta NODE-REPRESENTATIVE-famsa.fa

6.Build the tree using FastTree http://www.microbesonline.org/fasttree/

    FastTree -nt nirS-muscle-alignment.fa > nirS_fasttree.tre

7.Now you have everything that you need to build a Phyloseq object using R.  Make sure that you have R version 3.4, bioconductor 3.6 and phyloseq 1.22.  There are bugs in the earlier version of phyloseq.  You don't need to have any metadata but it sure does help. The sample names in the metadata file must be exactly the same as those contained in the PHYLOSEQ-MATRIX-OBJECT and you will need a header for each column.  

    library("phyloseq")
    library("ape")
    mat_nirS = read.table("PHYLOSEQ-MATRIX-OBJECT.txt", header = TRUE, sep = "\t", row.names = 1)
    tax_nirS = read.table("PHYLOSEQ-TAX-OBJECT.txt", header = TRUE, sep = ";", row.names = 1)
    meta_nirS = read.table("nrfA-map.txt", header = TRUE, sep = "\t", row.names = 1)
    tree_nirS = read.tree("nirS_fasttree.tre")    

    mat_nirS = as.matrix(mat_nirS)
    tax_nirS = as.matrix(tax_nirS)

    OTU = otu_table(mat_nirS, taxa_are_rows = TRUE)
    TAX = tax_table(tax_nirS)
    META = sample_data(meta_nirS)

    nrfA_physeq = phyloseq(OTU,TAX,META,tree_nirS)

## Update - I am working toward trying to cluster amino acid sequences for the short reads using MED.  This involves the dangerous practice of translating the short reads using prodigal.. 

1.Start out by using prodigal to call amino acid sequences for all short reads - this takes a while.

	prodigal -i sequences-to-decompose.fa -a sequences-to-decompose.faa -n -p meta



