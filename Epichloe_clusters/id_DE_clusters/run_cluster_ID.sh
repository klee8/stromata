# Identify clusters of genes that have at least 2-fold change

#  Cluster permutations
#  usage: Rscript --vanilla cluster_ID.R <basefilename>  <cluster_fun.R> <datadir> <resultsdir>
#  output: *_cluster_results.txt, _permuted_DE_cluster_results.txt, STR_PS_permuted_DE_clusters.pdf

#!usr/bin/bash

SCRIPTDIR="/media/kate/Massey_linux_onl/projects/stromata_analysis/Epichloe_clusters/id_DE_clusters/R"
DATADIR="/media/kate/Massey_linux_onl/projects/results/stromata/Epichloe_clusters/rfmt_core_gene_sets/"
RESDIR="/media/kate/Massey_linux_onl/projects/results/stromata/Epichloe_clusters/id_DE_clusters/"

cd $RESDIR

mkdir STR_PS INF_PS STR_INF

#####   Cluster permutations
cd STR_PS
Rscript --vanilla $SCRIPTDIR/cluster_ID.R STR_PS 10000 $SCRIPTDIR/cluster_fun.R $DATADIR $RESDIR
cd ../INF_PS
Rscript --vanilla $SCRIPTDIR/cluster_ID.R INF_PS 10000 $SCRIPTDIR/cluster_fun.R $DATADIR $RESDIR
cd ../STR_INF
Rscript --vanilla $SCRIPTDIR/cluster_ID.R STR_INF 10000 $SCRIPTDIR/cluster_fun.R $DATADIR $RESDIR
cd ..
