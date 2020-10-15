#!usr/bin/bash

mkdir ann_alns
    
for i in `cat elymi_core_gene_names.txt`
#for i in `cat long_list.txt`
#for i in FUN_000115
#for i in FUN_000205
do
    echo  $i
    rm temp
    # grab line from aln.fa file
    LINE=`cat align/$i.aln.fna | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | grep 'EelyNfE728' | cut -f 2`
    # get positions of gaps in the line
    POS=`cat align/$i.aln.fna | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | grep 'EelyNfE728' | cut -f 2 | grep -b -o '-' > gaps`
    sed -i 's/:-//g' gaps
    # get sequence position from gff3 file
    GENE=`grep $i Epichloe_elymi_NfE728.gff3 | grep gene | awk '{print $4":"$5}'`
    # get sequence orientation from gff3 file
    ORI=`grep $i Epichloe_elymi_NfE728.gff3 | grep gene | awk '{print $7}'`
    # grab exon info from gff3 file
    rm tmp_gff
    grep $i Epichloe_elymi_NfE728.gff3 | grep exon | awk '{print $4" "$5}' >> tmp_gff
    LENGTH=${#LINE}
    perl insert_exon_info.pl tmp_gff gaps $LENGTH $ORI $GENE > ann_alns/$i.ann.txt
    cat temp align/$i.aln.fna > ann_alns/$i.ann.aln.fna
done



# sort fasta files into the same order
cd ann_alns

for i in `cat ../elymi_core_gene_names.txt`
#for i in FUN_000205
do
    # linearise the files
    awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < $i.ann.aln.fna > $i.aln.lin.fna
    # sort by genome list
    for j in `cat ../genome_tags.txt`; do grep $j $i.aln.lin.fna >> $i.aln.sorted.lin.fna; done
    # re-format to fasta
    tr "\t" "\n" < $i.aln.sorted.lin.fna > $i.aln.sorted.fna
    # remove intermediate files
    rm $i.aln.lin.fna
    rm $i.aln.sorted.lin.fna
    # change the 'E' in the annotation line to 'C'
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
rm ann_alns/*.ann.aln.fna

# replace space with underscore in headers (for macvector)
for i in ann_alns/*.fna; do sed -i 's/ /_/g' $i; done
for i in upstream/*.fna; do sed -i 's/ /_/g' $i; done
for i in downstream/*.fna; do sed -i 's/ /_/g' $i; done


# add fasta with upstream and downstream sequences
# pad up and downstream strings to ensure alignment is not displaced
# https://stackoverflow.com/questions/4409399/padding-characters-in-printf
padlength=9
pad="---------"
for i in `cat elymi_core_gene_names.txt`
#for i in FUN_000205
do
    # linearise the sorted alignment and the up and downstream files
    awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < ann_alns/$i.aln.sorted.fna > ann_alns/$i.lin.fna
    awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < upstream/$i.fna > upstream/$i.lin.fna
    awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < downstream/$i.fna > downstream/$i.lin.fna
    # grab the sequences for each file for each header
    for j in `cat genome_tags.txt`
    do 
        TEMP1=`grep -m 1 $j ann_alns/$i.lin.fna | awk '{print $1;}'`
        TEMP2=`grep -m 1 $j upstream/$i.lin.fna | awk '{print $2;}'`
        TEMP3=`grep -m 1 $j ann_alns/$i.lin.fna | awk '{print $2;}'`
        TEMP4=`grep -m 1 $j downstream/$i.lin.fna | awk '{print $2;}'`
        # set output for annotation line
        #if [ $j = "EelyNFe728" ]; then TEMP2="AAAAAAAAA"; TEMP4="AAAAAAAAA"; fi
        # pad additional sequences where they are shorter than expected
        PAD2=`printf '%*.*s%s%s' 0 $((padlength - ${#TEMP2} )) "$pad" "$TEMP2"`
        PAD4=`printf '%s%*.*s%s' "$TEMP4" 0 $((padlength - ${#TEMP4} )) "$pad"`
        NEWLINE="$TEMP1\t$PAD2$TEMP3$PAD4"
#        echo -e "$TEMP1\t$PAD2\t$PAD4"
        echo -e $NEWLINE >> ann_alns/$i.updown.lin.fna
    done
    # change updown.lin.fna to fasta format
    tr "\t" "\n" < ann_alns/$i.updown.lin.fna > ann_alns/$i.updown.fna
done

rm ann_alns/*.lin.fna
rm upstream/*.lin.fna
rm downstream/*.lin.fna
