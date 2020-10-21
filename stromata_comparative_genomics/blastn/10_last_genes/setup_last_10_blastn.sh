 # setup all the blast alignments of last 10 core genes that either need a new gff start/end or to use EfeE2368 gene model

#!usr/bin/bash

rm  last_10_genes.fna

## Get the six genes from new Eelymi gene models from Dan
for i in `cat new_Eely_genes.txt`
do
    CONTIG=`grep -m 1 $i all_new_Eel_genes.gff | awk '{print $1}'`
    START=`grep -m 1 $i all_new_Eel_genes.gff | awk '{print $4}'` 
    END=`grep -m 1 $i all_new_Eel_genes.gff | awk '{print $5}'`
    ORIENTATION=`grep -m 1 $i all_new_Eel_genes.gff | awk '{print $7}'`
    # correct for reverse strand
    if [ $ORIENTATION = "-" ]; then ORI="--reverse-complement"; else ORI=""; fi
    # grab alignment from samtools
    samtools faidx ../../exonerate/exonerate_genomes/EelyNfE728.masked.fa $CONTIG:$START-$END $ORI --mark-strand sign >> last_10_genes.fna
    # add gene name to header
    sed -i "s/$CONTIG:$START-$END/$i EelyNfE728 $CONTIG:$START-$END/g"  last_10_genes.fna
done


## Get last four genes from Epichloe festucae E2368
for i in `cat new_Efe_genes.txt`
do 
#    echo $i
    CONTIG=`grep -m 1 $i Epichloe_festucae_E2368.gff3 | awk '{print $1}'`
    START=`grep -m 1 $i Epichloe_festucae_E2368.gff3 | awk '{print $4}'` 
    END=`grep -m 1 $i Epichloe_festucae_E2368.gff3 | awk '{print $5}'`
    ORIENTATION=`grep -m 1 $i Epichloe_festucae_E2368.gff3 | awk '{print $7}'`
    # correct for reverse strand
    if [ $ORIENTATION = "-" ]; then ORI="--reverse-complement"; else ORI=""; fi
    # grab alignment from samtools
    samtools faidx ../../exonerate/exonerate_genomes/EfesE2368.masked.fa $CONTIG:$START-$END $ORI --mark-strand sign >> last_10_genes.fna
    # add gene name to header
    sed -i "s/$CONTIG:$START-$END/$i EfeE2368 $CONTIG:$START-$END/g"  last_10_genes.fna
done

# add Eel homolog names back into fasta file so genes can be identified more easily later
sed -i 's/Efe2368_001458/FUN_000098 (used_query:Efe2368_001458)/g' last_10_genes.fna
sed -i 's/Efe2368_006269/FUN_000194 (used_query:Efe2368_006269)/g' last_10_genes.fna
sed -i 's/Efe2368_003497/FUN_000850 (used_query:Efe2368_003497)/g' last_10_genes.fna
sed -i 's/Efe2368_005841/FUN_007013 (used_query:Efe2368_005841)/g' last_10_genes.fna



# set up blast dirs 
#for i in `cat ../genome_tag_list.txt`; do mkdir $i; cp ../../exonerate/align/$i/$i.masked.fa $i; done

# setup scripts to run blast and get best hits
for i in `cat ../genome_tag_list.txt`; do 
	echo '#!usr/bin/bash' > ../$i/$i.blastn.sh 
	echo "makeblastdb -in $i.masked.fa -dbtype nucl"  >> ../$i/$i.blastn.sh 
	echo "blastn -query ../10_last_genes/last_10_genes.fna -db $i.masked.fa -evalue 5 -outfmt 6 -out $i.aln.txt" >> ../$i/$i.blastn.sh 
done

# set screens to run the blasts
for i in `cat ../genome_tag_list.txt`; do cd ../$i; screen -dmS $i bash $i.blastn.sh; cd ../10_last_genes; done

