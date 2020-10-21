# setup script for exonerate - to identify homologs between species
#!usr/bin/bash

# set up directories
for i in `cat genome_tag_list.txt`; do mkdir align/$i; done


# create exonerate scripts
for i in `cat genome_tag_list.txt`; do echo '''#!usr/bin/bash

# run exonerate
exonerate -q ../elymi_core-genes.fa -t TEMP.masked.fa -m protein2genome --percent 50 --showtargetgff --bestn 10 --softmasktarget --verbose 5 > TEMP.ex.out

# grab the gff3 output
# seqname source feature start end score strand frame attributes
grep 'gene_id' TEMP.ex.out > TEMP.core_genes.gff3''' > align/$i/$i.ex.sh; sed -i "s/TEMP/$i/g"  align/$i/$i.ex.sh; done


# add in the genomes
cd exonerate_genomes
for i in `cat ../genome_tag_list.txt`; do cp $i.masked.fa ../align/$i/; done
cd ..


