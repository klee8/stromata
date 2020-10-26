#!usr/bin/bash

mkdir blastn_upstream
mkdir blastn_downstream
mkdir blastn_query
#for i in `cat elymi_core_gene_names.txt`
#for i in `cat control_gene_list.txt`
for i in `cat last_10_genelist.txt`
do
    echo $i
    rm blastn_query/$i.positions.txt
    rm blastn_query/$i.fna
    echo "AlnType Gene Genome NumHits Contig Start End" > blastn_query/$i.positions.txt
    for j in `cat genome_tag_list.txt`
    do
	# for the E.elymiNfE728 blastn_query grab the gff3 gene positions
        if [ "$j" = "EelyNfE728" ]
        then
            CONTIG=`grep -m 1 $i Epichloe_elymi_NfE728.gff3 | awk '{print $1}'`
            START=`grep -m 1 $i Epichloe_elymi_NfE728.gff3 | awk '{print $4}'` 
            END=`grep -m 1 $i Epichloe_elymi_NfE728.gff3 | awk '{print $5}'`
            ORIENTATION=`grep -m 1 $i Epichloe_elymi_NfE728.gff3 | awk '{print $7}'`
            ALIGN="NONE"
            echo "$ALIGN $i $j $NUMHITS $CONTIG $START $END" >> blastn_query/$i.positions.txt
            # correct for reverse strand
            if [ $ORIENTATION = "-" ]; then ORI="--reverse-complement"; else ORI=""; fi
            # grab alignment from samtools
            samtools faidx ../exonerate/exonerate_genomes/$j.masked.fa $CONTIG:$START-$END $ORI --mark-strand sign >> blastn_query/$i.fna
            # get blastn_upstream and blastn_downstream sequence 
            if [ $ORIENTATION = "+" ]; then UP1=$(expr $START - 99); UP2=$(expr $START - 1); DOWN1=$(expr $END + 1); DOWN2=$(expr $END + 99);
            else UP1=$(expr $END + 1); UP2=$(expr $END + 99); DOWN1=$(expr $START - 99); DOWN2=$(expr $START - 1); fi
            #for dir in blastn_upstream blastn_downstream; do echo -e ">$LIFESTYLE B $j $CONTIG:$START-$END\n" >> $dir/$i.fna; done
            samtools faidx ../exonerate/exonerate_genomes/$j.masked.fa $CONTIG:$UP1-$UP2 $ORI --mark-strand sign >> blastn_upstream/$i.fna
            samtools faidx ../exonerate/exonerate_genomes/$j.masked.fa $CONTIG:$DOWN1-$DOWN2 $ORI --mark-strand sign >> blastn_downstream/$i.fna 
            LENGTH=`expr $END - $START`
#            echo -e "$i\t$j\t$ALIGN\t$ORIENTATION\t$CONTIG\t$UP1\t$UP2\t$START\t$END\t$DOWN1\t$DOWN2\t$LENGTH"
            # add in extra info into fasta header (including sexual/asexual/unknown, blastn_query tag and species name)
            sed -i "s/$CONTIG:$START-$END/S Q EelyNfE728 $CONTIG:$START-$END/g"  blastn_query/$i.fna
            sed -i "s/$CONTIG:$UP1-$UP2/S Q EelyNfE728 $CONTIG:$UP1-$UP2/g"  blastn_upstream/$i.fna
            sed -i "s/$CONTIG:$DOWN1-$DOWN2/S Q EelyNfE728 $CONTIG:$DOWN1-$DOWN2/g"  blastn_downstream/$i.fna
            continue
	else
            CONTIG=`grep -m 1 $i ../blastn/$j/$j.bestaln.txt | awk '{print $2}'`
            START=`grep -m 1 $i ../blastn/$j/$j.bestaln.txt | awk '{print $9}'` 
            END=`grep -m 1 $i ../blastn/$j/$j.bestaln.txt | awk '{print $10}'`
            ORIENTATION=`if [ $START -lt $END ]; then echo "+"; else echo "-"; fi`
            ALIGN="BLAST"
            echo "$ALIGN $i $j $NUMHITS $CONTIG $START $END" >> blastn_query/$i.positions.txt
            # correct for reverse strand
            if [ $ORIENTATION = "-" ]; then ORI="--reverse-complement"; BEG=$START; START=$END; END=$BEG; else ORI=""; fi
            # grab alignment from samtools
            samtools faidx ../exonerate/exonerate_genomes/$j.masked.fa $CONTIG:$START-$END $ORI --mark-strand sign >> blastn_query/$i.fna
            # get blastn_upstream and blastn_downstream sequence 
            if [ $ORIENTATION = "+" ]; then UP1=$(expr $START - 99); UP2=$(expr $START - 1); DOWN1=$(expr $END + 1); DOWN2=$(expr $END + 99);
            else UP1=$(expr $END + 1); UP2=$(expr $END + 99); DOWN1=$(expr $START - 99); DOWN2=$(expr $START - 1); fi
            #for dir in blastn_upstream blastn_downstream; do echo -e ">$LIFESTYLE B $j $CONTIG:$START-$END\n" >> $dir/$i.fna; done
            # note when the 
            if [ $UP2 -ne 0 ]; then samtools faidx ../exonerate/exonerate_genomes/$j.masked.fa $CONTIG:$UP1-$UP2 $ORI --mark-strand sign >> blastn_upstream/$i.fna; fi
            if [ $DOWN2 -ne 0 ]; then samtools faidx ../exonerate/exonerate_genomes/$j.masked.fa $CONTIG:$DOWN1-$DOWN2 $ORI --mark-strand sign >> blastn_downstream/$i.fna; fi
            LENGTH=`expr $END - $START`
#            echo -e "$i\t$j\t$ALIGN\t$ORIENTATION\t$CONTIG\t$UP1\t$UP2\t$START\t$END\t$DOWN1\t$DOWN2\t$LENGTH"
            # add in extra info into fasta header (including sexual/asexual/unknown, hit type and species name)
            LIFESTYLE=`grep -m 1 $j lifestyle.txt  | awk '{print $2}'` 
            sed -i "s/$CONTIG:$START-$END/$LIFESTYLE B $j $CONTIG:$START-$END/g"  blastn_query/$i.fna
            sed -i "s/$CONTIG:$UP1-$UP2/S Q $j $CONTIG:$UP1-$UP2/g"  blastn_upstream/$i.fna
            sed -i "s/$CONTIG:$DOWN1-$DOWN2/S Q $j $CONTIG:$DOWN1-$DOWN2/g"  blastn_downstream/$i.fna
        fi
    done
done

rm run_muscle.sh
#for i in `cat control_gene_list.txt`
for i in `cat last_10_genelist.txt`
    do echo "screen -dmS $i muscle -in blastn_query/$i.fna -out blastn_align/$i.aln.fna" >> run_muscle.sh
done

