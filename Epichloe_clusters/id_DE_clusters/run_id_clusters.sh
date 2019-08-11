#!usr/bin/bash

SCRIPTDIR=
DATADIR=


######### RUN PERL script to identify and graph clusters

mkdir STR_PS
cd STR_PS
perl ../id_DE_clusters.pl ../../rfmt_core_gene_sets/core_genes_STR_PS_rfmt.txt
cd ..

mkdir STR_INF
cd STR_INF
perl ../id_DE_clusters.pl ../../rfmt_core_gene_sets/core_genes_STR_INF_rfmt.txt
cd ..

mkdir INF_PS 
cd INF_PS 
perl ../id_DE_clusters.pl ../../rfmt_core_gene_sets/core_genes_INF_PS_rfmt.txt
cd ..

# re-format output to append to Daniel's Annotations
R --slave --vanilla < summarise_clusters_by_orthogroup.R


######## RUN R script to carry out permutations of cluster identification
