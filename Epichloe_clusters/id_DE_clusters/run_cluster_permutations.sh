# Randomly permute DE expression values among genes for each species
#    check the number of clusers that are found when DE is randomised

#  Cluster permutations
#  usage: Rscript --vanilla clusters_permutation.R <basefilename> <number of simulations>
#  output: *_cluster_results.txt, _permuted_DE_cluster_results.txt, STR_PS_permuted_DE_clusters.pdf

#!usr/bin/bash

SCRIPTDIR="/media/kate/Massey_linux_onl/projects/stromata_analysis/Epichloe_clusters/id_DE_clusters/R"
DATADIR="/media/kate/Massey_linux_onl/projects/results/stromata/Epichloe_clusters/rfmt_core_gene_sets/"
RESDIR="/media/kate/Massey_linux_onl/projects/results/stromata/Epichloe_clusters/id_DE_clusters/"

cd $RESDIR

mkdir STR_PS INF_PS STR_INF

#####   Cluster permutations
cd STR_PS
Rscript --vanilla $SCRIPTDIR/clusters_permutation.R STR_PS 10000 $SCRIPTDIR/cluster_fun.R $DATADIR $RESDIR
cd ../INF_PS
Rscript --vanilla $SCRIPTDIR/clusters_permutation.R INF_PS 10000 $SCRIPTDIR/cluster_fun.R $DATADIR $RESDIR
cd ../STR_INF
Rscript --vanilla $SCRIPTDIR/clusters_permutation.R STR_INF 10000 $SCRIPTDIR/cluster_fun.R $DATADIR $RESDIR
cd ..
