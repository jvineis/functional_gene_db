# Creation of a database for important functional genes
### Steps to create a fasta file of functional genes from 6000 ncbi "Gold" genomes and 8000 genomes derived from https://www.nature.com/articles/s41564-017-0012-7. At the time of the creation of this README, I was working from here /Users/joevineis/Documents/BOWEN/functional_gene_db.

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
    sed 's/_contigs//g' dsrB_ncbi_parks.fa > dsrB
    mv dsrB dsrB_ncbi_parks.fa

9.Fix the deflines that are produced by anvio using the handy ncbi_parks_name_fix.py.  eg

    python ncbi_parks_name_fix.py nirS_ncbi_parks.fa nirS
    mv nirS nirS_ncbi_parks.fa 

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


