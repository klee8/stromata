# Randomly permute DE expression values among genes for each species
#    check the number of orthologs that overlap between species

#  DE permutations
#  usage:  Rscript --vanilla DE_permutations.R <reformatted core gene file>
#  output: *_permuted_DE_results.txt,  *_permuted_DE_clusters.pdf

#!usr/bin/bash


SCRIPTDIR="/media/kate/Massey_linux_onl/projects/stromata_analysis/Epichloe_clusters/id_DE_clusters"
DATADIR="/media/kate/Massey_linux_onl/projects/results/stromata/Epichloe_clusters/rfmt_core_gene_sets/"
RESDIR="/media/kate/Massey_linux_onl/projects/results/stromata/Epichloe_clusters/id_DE_clusters/"

cd $RESDIR

mkdir STR_PS INF_PS STR_INF

#####   DE permutations
cd STR_PS
Rscript --vanilla $SCRIPTDIR/DE_permutations.R STR_PS 10000 $SCRIPTDIR/DE_fun.R $DATADIR $RESDIR
cd ../INF_PS
Rscript --vanilla $SCRIPTDIR/DE_permutations.R INF_PS 10000 $SCRIPTDIR/DE_fun.R $DATADIR $RESDIR
cd ../STR_INF
Rscript --vanilla $SCRIPTDIR/DE_permutations.R STR_INF 10000 $SCRIPTDIR/DE_fun.R $DATADIR $RESDIR

