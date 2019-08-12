# Kate Lee 2019
# format fasta files for input to SMURF
#!usr/bin/bash

# Set directories
DATADIR="media/kate/Massey_linux_onl/projects/data/gene_model_sets/"
RESDIR="/media/kate/Massey_linux_onl/projects/STROMATA/results/Epichloe_clusters/SMURF/"
SCRIPTDIR="/media/kate/Massey_linux_onl/projects/STROMATA/stromata_analysis/Epichloe_clusters/SMURF/"

mkdir $RESDIR/data 
Rscript --vanilla < $SCRIPTDIR/format_SMURF_input_files.R $DATADIR $RESDIR
