#!usr/bin/bash

# runs perl and python scripts to identify clusters and groups of clusters
# output includes
#    clusters.txt (one line for each cluster)
#    cluster_groups.txt (numbers the groups of clusters and shows how many species
#         share an ortholog(s) across clusters)
#    log_cluster_finding.txt (log output of the perl cluster identification)

# summarise clusters by ortholog uses:
# flat_ortho_file (from proteinortho output re-formated using David Winters script
# https://gist.github.com/dwinter/eeb0126e8a2243a87cca2533a6b6c633#call-orthologs
# outputs:  *_group_file_for_ann.txt (one line per ortholog to merge with annotations)


SCRIPTDIR="/media/kate/Massey_linux_onl/projects/stromata_analysis/Epichloe_clusters/id_DE_clusters"
RESDIR="/media/kate/Massey_linux_onl/projects/results/stromata/Epichloe_clusters/id_DE_clusters"

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
