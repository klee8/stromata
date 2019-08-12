#!ust/bin/bash

DATADIR="/media/kate/Massey_linux_onl/projects/STROMATA/stromata_analysis/Epichloe_stromata_DE/3_DESeq2/"
RESDIR="/media/kate/Massey_linux_onl/projects/STROMATA/results/Epichloe_clusters/rfmt_core_gene_sets/"

SCRIPTDIR="/media/kate/Massey_linux_onl/projects/STROMATA/stromata_analysis/Epichloe_clusters/rfmt_core_gene_sets/"

Rscript --vanilla < $SCRIPTDIR/reformat_core_gene_set.R $DATADIR $RESDIR
Rscript --vanilla < $SCRIPTDIR/reformat_core_gene_set.R $DATADIR $RESDIR
Rscript --vanilla < $SCRIPTDIR/reformat_core_gene_set.R $DATADIR $RESDIR
