# Kate Lee 2019
# Bash script to reformat core gene set output for DE gene and cluster analysis


#!usr/bin/bash

DATADIR="/media/kate/Massey_linux_onl/projects/results/stromata/Epichloe_stromata_DE/3_DEseq2/"
RESDIR="/media/kate/Massey_linux_onl/projects/results/stromata/Epichloe_clusters/rfmt_core_gene_sets/"
SCRIPTDIR="/media/kate/Massey_linux_onl/projects/stromata_analysis/Epichloe_clusters/id_DE_clusters/R"

Rscript --vanilla < $SCRIPTDIR/reformat_core_gene_set.R $DATADIR $RESDIR
