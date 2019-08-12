# Randomly permute DE expression values among genes for each species
#    check the number of orthologs that overlap between species
#    check the number of clusers that are found when DE is randomised

#  DE permutations
#  usage:  Rscript --vanilla DE_permutations.R <reformatted core gene file>
#  output: *_permute_DE_results.txt,  *_permuted_DE_clusters.pdf

#  Cluster permutations
#  usage: Rscript --vanilla clusters_permutation.R <basefilename> <number of simulations>
#  output:

#!usr/bin/bash


SCRIPTDIR="/media/kate/Massey_linux_onl/projects/STROMATA/stromata_analysis/Epichloe_clusters/id_DE_clusters"
RESDIR="/media/kate/Massey_linux_onl/projects/STROMATA/results/Epichloe_clusters/id_DE_clusters"

cd $RESDIR

#####   DE permutations
cd STR_PS
Rscript --vanilla $SCRIPTDIR/DE_permutations.R STR_PS 1000
cd ../INF_PS
Rscript --vanilla $SCRIPTDIR/DE_permutations.R INF_PS 1000
cd ../STR_INF
Rscript --vanilla $SCRIPTDIR/DE_permutations.R STR_INF 1000

#####   Cluster permutations
cd ../STR_PS
Rscript --vanilla $SCRIPTDIR/clusters_permutation.R STR_PS 1000
cd ../INF_PS
Rscript --vanilla $SCRIPTDIR/clusters_permutation.R INF_PS 1000
cd ../STR_INF
Rscript --vanilla $SCRIPTDIR/clusters_permutation.R STR_INF 1000
cd ..
