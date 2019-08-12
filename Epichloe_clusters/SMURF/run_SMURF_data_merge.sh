# Kate Lee 2019
# re-format SMURF output to merge it with core gene set
#!usr/bin/bash


# Set directories
ORTHO="/media/kate/Massey_linux_onl/projects/data/gene_model_sets/proteinortho/flat_ortho_file.txt"
RESDIR="/media/kate/Massey_linux_onl/projects/STROMATA/results/Epichloe_clusters/SMURF/"
SCRIPTDIR="/media/kate/Massey_linux_onl/projects/STROMATA/stromata_analysis/Epichloe_clusters/SMURF/"

# Reformat SMURF output
perl $SCRIPTDIR/format_SMURF_output.pl $RESDIR/E.elymi_NfE728_results/Secondary-Metabolite-Clusters.txt $RESDIR/E.elymi_NfE728_results/SMURF_ely_rfmt.txt
perl $SCRIPTDIR/format_SMURF_output.pl $RESDIR/E.festucae_E2368_results/Secondary-Metabolite-Clusters.txt $RESDIR/E.festucae_E2368_results/SMURF_fes_rfmt.txt
perl $SCRIPTDIR/format_SMURF_output.pl $RESDIR/E.typhina_E8_results/Secondary-Metabolite-Clusters.txt $RESDIR/E.typhina_E8_results/SMURF_typ_rfmt.txt

Rscript --vanilla $SCRIPTDIR/summarise_SMURF.R $ORTHO $RESDIR
