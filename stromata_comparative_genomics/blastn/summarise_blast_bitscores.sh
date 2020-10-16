# grab all the blast bitscores for each core gene for each genome in the dir 
#!usr/bin/bash

# if you need to start over remove any earlier output files
for i in `cat genome_tag_list.txt`; do rm $i/*bestaln.txt; done 
rm all_bitscores.txt

# get the best hits
for i in `cat genome_tag_list.txt`; do for j in `cat elymi_core_gene_names.txt`; do grep -m 1 $j $i/$i.aln.txt >> $i/$i.bestaln.txt ; done; done 

# paste the bitscores into one file
cp elymi_core_gene_names.txt all_bitscores.txt
for i in `cat genome_tag_list.txt`; do cut -f 12 $i/$i.bestaln.txt > tmp1; paste all_bitscores.txt tmp1 > tmp2; mv tmp2 all_bitscores.txt; done
echo -ne "elymi_core_genes\t" > headerline 
for i in `cat genome_tag_list.txt`; do echo -ne "$i\t" >> headerline; done; echo >> headerline
cat headerline all_bitscores.txt > temp
mv temp all_bitscores.txt
rm headerline
rm tmp1
