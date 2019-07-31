## Identify clusters in differentially expressed data

#### Kate Lee July 2019


### AIM:
Identify clusters of differentially expressed genes that are shared between species.

### INPUT:
Input has columns for each gene in analysis (reformatted using reformat_core_gene_set.Rmd):  
"contig", "start", "stop", "gene_id","orthogroup", "log2fc", "lfcSE", "svalue_1", "svalue_2", "species"

### USAGE:
perl id_DE_clusters.pl <reformated_core_gene_set_file> > log.txt

### WORKFLOW
* Iterates by contig and gene postion through differential expression output from core gene DE analysis (rows for each gene on each species). 
    * Identifies genes that could initiate or extend a cluster based on their fold change and svalues
* Iterates through DE values a second time. 
    * Initiates clusters where there is a gene that could initiate or extend
    * Allows for a given number of gaps
    * Evaluates potential cluster when the direction of differential expression changes, or the allowed number of gaps is exceeded, or the contig changes.
    * Trailing gaps are trimmed off (gaps are only allowed within the cluster 
    * Cluster must have a minimum number of genes 
    * If all requirements are met, cluster is recored
* Iterates through clusters to find any that are linked across species
    * For each ortholog in the cluster, checks its presence in clusters in other species
    * For each ortholog not in original cluster, but present in other species, checks for more linked clusters
    * (On the off chance that clusters from the same species are linked across two contigs by clusters in other species, make a third pass with any unchecked orthologs).
    * Writes out table of clusters identified
* Iterates through groups of clusters
    * identifies shared orthologs
    * writes table of cluster groups with flat file (i.e. one row for each species and ortholog) 
    * Writes python script to graphs cluster with added annotation to show end of contigs, gaps and shared orthologs using biopython.
    