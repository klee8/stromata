 #!usr/bin/bash

mkdir blastn_ann_alns



# grab the elymi exon annotations from the gff file and add an annotation line to the alignment
# this will help identify intron-exon boundaries
#for i in `cat elymi_core_gene_names.txt`
#for i in `cat control_gene_list.txt`
for i in `cat last_10_genelist.txt`
do
    echo  $i
    rm temp
    # grab line from aln.fa file
    LINE=`cat blastn_align/$i.aln.fna | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | grep 'EelyNfE728' | cut -f 2`
    # get positions of gaps in the line
    POS=`cat blastn_align/$i.aln.fna | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | grep 'EelyNfE728' | cut -f 2 | grep -b -o '-' > gaps`
    sed -i 's/:-//g' gaps
    # get sequence position from gff3 file
    GENE=`grep $i Epichloe_elymi_NfE728.gff3 | grep gene | awk '{print $4":"$5}'`
    #GENE=`grep $i Epichloe_elymi_NfE728.gff3 | grep gene | awk '{print $4":"$5}'`    # alt for last 10 genes with new models
    # get sequence orientation from gff3 file
    ORI=`grep $i Epichloe_elymi_NfE728.gff3 | grep gene | awk '{print $7}'`
    #ORI=`grep $i Epichloe_elymi_NfE728.gff3 | grep gene | awk '{print $7}'`   # alt for last 10 genes with new models
    # grab exon info from gff3 file
    rm tmp_gff
    grep $i Epichloe_elymi_NfE728.gff3 | grep exon | awk '{print $4" "$5}' >> tmp_gff
    #grep $i Epichloe_elymi_NfE728.gff3 | grep exon | awk '{print $4" "$5}' >> tmp_gff   # alt for last 10 genes with new models
    LENGTH=${#LINE}
    # use perl script to create annotation line (creates 'temp' file)
    perl insert_exon_info.pl tmp_gff gaps $LENGTH $ORI $GENE > blastn_ann_alns/$i.blastn_ann.txt
    # paste new annotation line into muscle alignment
    cat temp blastn_align/$i.aln.fna > blastn_ann_alns/$i.blastn_ann.aln.fna
#    cat temp
done



# sort fasta files into the same order
cd blastn_ann_alns

#for i in `cat ../elymi_core_gene_names.txt`
#for i in `cat ../control_gene_list.txt`
for i in `cat ../last_10_genelist.txt`
#for i in FUN_000205
do
    # linearise the files
    awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < $i.blastn_ann.aln.fna > $i.aln.lin.fna
    # sort by genome list
    for j in `cat ../genome_tags.txt`; do grep $j $i.aln.lin.fna >> $i.aln.sorted.lin.fna; done
    # re-format to fasta
    tr "\t" "\n" < $i.aln.sorted.lin.fna > $i.aln.sorted.fna
    # remove intermediate files
    rm $i.aln.lin.fna
    rm $i.aln.sorted.lin.fna
    # change the 'E' in the blastn_annotation line to 'C'
    sed -i 's/EE/CC/g' $i.aln.sorted.fna
    sed -i 's/-E-/-C-/g' $i.aln.sorted.fna
    sed -i 's/-E/-C/g' $i.aln.sorted.fna
    sed -i 's/E-/C-/g' $i.aln.sorted.fna
    sed -i 's/NE/NC/g' $i.aln.sorted.fna
    sed -i 's/EN/CN/g' $i.aln.sorted.fna
    sed -i 's/CE/CC/g' $i.aln.sorted.fna
done
cd ..

# remove intermediate files
rm blastn_ann_alns/*.blastn_ann.aln.fna

# replace space with underscore in headers (for macvector)
for i in blastn_ann_alns/*.fna; do sed -i 's/ /_/g' $i; done
for i in blastn_upstream/*.fna; do sed -i 's/ /_/g' $i; done
for i in blastn_downstream/*.fna; do sed -i 's/ /_/g' $i; done


# add fasta with blastn_upstream and blastn_downstream sequences
# pad up and blastn_downstream strings to ensure blastn_alignment is not displaced
# https://stackoverflow.com/questions/4409399/padding-characters-in-printf
padlength=99
pad="---------------------------------------------------------------------------------------------------"
#for i in `cat elymi_core_gene_names.txt`
#for i in `cat control_gene_list.txt`
for i in `cat last_10_genelist.txt`
#for i in FUN_000205
do
    # linearise the sorted blastn_alignment and the up and blastn_downstream files
    awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < blastn_ann_alns/$i.aln.sorted.fna > blastn_ann_alns/$i.lin.fna
    awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < blastn_upstream/$i.fna > blastn_upstream/$i.lin.fna
    awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < blastn_downstream/$i.fna > blastn_downstream/$i.lin.fna
    # grab the sequences for each file for each header
    for j in `cat genome_tags.txt`
    do 
        TEMP1=`grep -m 1 $j blastn_ann_alns/$i.lin.fna | awk '{print $1;}'`
        TEMP2=`grep -m 1 $j blastn_upstream/$i.lin.fna | awk '{print $2;}'`
        TEMP3=`grep -m 1 $j blastn_ann_alns/$i.lin.fna | awk '{print $2;}'`
        TEMP4=`grep -m 1 $j blastn_downstream/$i.lin.fna | awk '{print $2;}'`
        # pad additional sequences where they are shorter than expected
        PAD2=`printf '%*.*s%s%s' 0 $((padlength - ${#TEMP2} )) "$pad" "$TEMP2"`
        PAD4=`printf '%s%*.*s%s' "$TEMP4" 0 $((padlength - ${#TEMP4} )) "$pad"`
        NEWLINE="$TEMP1\t$PAD2$TEMP3$PAD4"
        echo -e $NEWLINE >> blastn_ann_alns/$i.updown.lin.fna
    done
    # change updown.lin.fna to fasta format
    tr "\t" "\n" < blastn_ann_alns/$i.updown.lin.fna > blastn_ann_alns/$i.updown.fna
done

