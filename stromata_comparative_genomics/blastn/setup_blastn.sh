# setup all the blast alignments of elymi core genes to the set of genomes for comparison

#!usr/bin/bash

# grab the EelyNfE728 core genes nucleotide sequences
rm  EelyNfE728.fna
#for i in `cat elymi_core_gene_names.txt`
for i in `cat control_gene_list.txt`
#for i in FUN_000066
#for i in FUN_000115
do
#    echo $i
    CONTIG=`grep -m 1 $i Epichloe_elymi_NfE728.gff3 | awk '{print $1}'`
    START=`grep -m 1 $i Epichloe_elymi_NfE728.gff3 | awk '{print $4}'` 
    END=`grep -m 1 $i Epichloe_elymi_NfE728.gff3 | awk '{print $5}'`
    ORIENTATION=`grep -m 1 $i Epichloe_elymi_NfE728.gff3 | awk '{print $7}'`
    # correct for reverse strand
    if [ $ORIENTATION = "-" ]; then ORI="--reverse-complement"; else ORI=""; fi
    # grab alignment from samtools
    samtools faidx ../exonerate/exonerate_genomes/EelyNfE728.masked.fa $CONTIG:$START-$END $ORI --mark-strand sign >> EelyNfE728.fna
    # add gene name to header
    sed -i "s/$CONTIG:$START-$END/$i EelyNfE728 $CONTIG:$START-$END/g"  EelyNfE728.fna
done

# set up blast dirs 
for i in `cat genome_tag_list.txt`; do mkdir $i; cp ../exonerate/align/$i/$i.masked.fa $i; done

# setup scripts to run blast and get best hits
for i in `cat genome_tag_list.txt`; do 
	echo '#!usr/bin/bash' > $i/$i.blastn.sh 
	echo "makeblastdb -in $i.masked.fa -dbtype nucl"  >> $i/$i.blastn.sh 
	echo "blastn -query ../EelyNfE728.fna -db $i.masked.fa -evalue 5 -outfmt 6 -out $i.aln.txt" >> $i/$i.blastn.sh 
done

# set screens to run the blasts
for i in `cat genome_tag_list.txt`; do cd $i; screen -dmS $i bash $i.blastn.sh; cd ..; done



