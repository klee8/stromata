#!usr/bin/bash

SCRIPTDIR="/media/kate/Massey_linux_onl/projects/STROMATA/stromata_analysis/Epichloe_clusters/id_DE_clusters"
RESDIR="/media/kate/Massey_linux_onl/projects/STROMATA/results/Epichloe_clusters/id_DE_clusters"

cd $RESDIR

######### RUN PERL script to identify and graph clusters

mkdir STR_PS
cd STR_PS
mkdir graphs
perl $SCRIPTDIR/id_DE_clusters.pl ../../rfmt_core_gene_sets/core_genes_STR_PS_rfmt.txt
cd graphs
for i in *.py; do python $i; done
cd ../..

mkdir STR_INF
cd STR_INF
mkdir graphs
perl $SCRIPTDIR/id_DE_clusters.pl ../../rfmt_core_gene_sets/core_genes_STR_INF_rfmt.txt
cd graphs
for i in *.py; do python $i; done
cd ../..

mkdir INF_PS
cd INF_PS
mkdir graphs
perl $SCRIPTDIR/id_DE_clusters.pl ../../rfmt_core_gene_sets/core_genes_INF_PS_rfmt.txt
cd graphs
for i in *.py; do python $i; done
cd ../..

# re-format output to append to Daniel's Annotations
R --slave --vanilla < $SCRIPTDIR/summarise_clusters_by_orthogroup.R
