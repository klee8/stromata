#!usr/bin/bash

mkdir upstream
mkdir downstream
mkdir query
for i in `cat elymi_core_gene_names.txt`
#for i in FUN_000430
#for i in FUN_000205
do
#    echo $i
#    rm query/$i.positions.txt
#    rm query/$i.fna
    echo "AlnType Gene Genome NumHits Contig Start End" > query/$i.positions.txt
    for j in `cat genome_tag_list.txt`
    do
	# for the E.elymiNfE728 query grab the gff3 gene positions
        if [ "$j" = "EelyNfE728" ]
        then
            CONTIG=`grep -m 1 $i Epichloe_elymi_NfE728.gff3 | awk '{print $1}'`
            START=`grep -m 1 $i Epichloe_elymi_NfE728.gff3 | awk '{print $4}'` 
            END=`grep -m 1 $i Epichloe_elymi_NfE728.gff3 | awk '{print $5}'`
            ORIENTATION=`grep -m 1 $i Epichloe_elymi_NfE728.gff3 | awk '{print $7}'`
            ALIGN="NONE"
            echo "$ALIGN $i $j $NUMHITS $CONTIG $START $END" >> query/$i.positions.txt
            # correct for reverse strand
            if [ $ORIENTATION = "-" ]; then ORI="--reverse-complement"; else ORI=""; fi
            # grab alignment from samtools
            samtools faidx ../exonerate/exonerate_genomes/$j.masked.fa $CONTIG:$START-$END $ORI --mark-strand sign >> query/$i.fna
            # get upstream and downstream sequence 
            if [ $ORIENTATION = "+" ]; then UP1=$(expr $START - 99); UP2=$(expr $START - 1); DOWN1=$(expr $END + 1); DOWN2=$(expr $END + 99);
            else UP1=$(expr $END + 1); UP2=$(expr $END + 99); DOWN1=$(expr $START - 99); DOWN2=$(expr $START - 1); fi
            samtools faidx ../exonerate/exonerate_genomes/$j.masked.fa $CONTIG:$UP1-$UP2 $ORI --mark-strand sign >> upstream/$i.fna
            samtools faidx ../exonerate/exonerate_genomes/$j.masked.fa $CONTIG:$DOWN1-$DOWN2 $ORI --mark-strand sign >> downstream/$i.fna 
            LENGTH=`expr $END - $START`
            echo -e "$i\t$j\t$ALIGN\t$ORIENTATION\t$CONTIG\t$UP1\t$UP2\t$START\t$END\t$DOWN1\t$DOWN2\t$LENGTH"
            # add in extra info into fasta header (including sexual/asexual/unknown, query tag and species name)
            sed -i "s/$CONTIG:$START-$END/S Q EelyNfE728 $CONTIG:$START-$END/g"  query/$i.fna
            sed -i "s/$CONTIG:$UP1-$UP2/S Q EelyNfE728 $CONTIG:$UP1-$UP2/g"  upstream/$i.fna
            sed -i "s/$CONTIG:$DOWN1-$DOWN2/S Q EelyNfE728 $CONTIG:$DOWN1-$DOWN2/g"  downstream/$i.fna
            continue
	fi
	
        # Grab exonerate alignments and put in a tempfile
        grep $i ../exonerate/align/$j/$j.core_genes.gff3 > temp
        NUMHITS=`wc -l temp | awk '{print $1}'`
        while IFS= read -r line
        do
            CONTIG=`echo $line | awk '{print $1}'`
            START=`echo $line | awk '{print $4}'`
            END=`echo $line | awk '{print $5}'`
	    ALIGN="Exonerate"
            ORIENTATION=`echo $line | awk '{print $7}'`
            echo "$ALIGN $i $j $NUMHITS $CONTIG $START $END" >> query/$i.positions.txt
            # correct for reverse strand
            if [ $ORIENTATION = "-" ]; then ORI="--reverse-complement"; else ORI=""; fi
            # grab alignment from samtools
            samtools faidx ../exonerate/exonerate_genomes/$j.masked.fa $CONTIG:$START-$END $ORI --mark-strand sign >> query/$i.fna
            # get upstream and downstream sequence 
            if [ $ORIENTATION = "+" ]; then  UP1=$(expr $START - 99); UP2=$(expr $START - 1); DOWN1=$(expr $END + 1); DOWN2=$(expr $END + 99); 
            else UP1=$(expr $END + 1); UP2=$(expr $END + 99); DOWN1=$(expr $START - 99); DOWN2=$(expr $START - 1); fi
            LENGTH=`expr $END - $START`
            echo -e "$i\t$j\t$ALIGN\t$ORIENTATION\t$CONTIG\t$UP1\t$UP2\t$START\t$END\t$DOWN1\t$DOWN2\t$LENGTH"
            samtools faidx ../exonerate/exonerate_genomes/$j.masked.fa $CONTIG:$UP1-$UP2 $ORI --mark-strand sign >> upstream/$i.fna 
            samtools faidx ../exonerate/exonerate_genomes/$j.masked.fa $CONTIG:$DOWN1-$DOWN2 $ORI --mark-strand sign >> downstream/$i.fna 
            # add in extra info into fasta header (including sexual/asexual/unknown, hit type and species name)
            LIFESTYLE=`grep -m 1 $j lifestyle.txt  | awk '{print $2}'` 
            sed -i "s/$CONTIG:$START-$END/$LIFESTYLE E $j $CONTIG:$START-$END/g"  query/$i.fna 
            sed -i "s/$CONTIG:$UP1-$UP2/$LIFESTYLE E $j $CONTIG:$UP1-$UP2/g"  upstream/$i.fna
            sed -i "s/$CONTIG:$DOWN1-$DOWN2/$LIFESTYLE E $j $CONTIG:$DOWN1-$DOWN2/g"  downstream/$i.fna
        done < "./temp"

	# if there are no exonerate alignments, get the top BLAST alignment hit (and any others on the same contig)
        if [ ! -s "./temp" ] 
        then
            CONTIG=`grep -m 1 $i ../blastn/$j/$j.bestaln.txt | awk '{print $2}'`
            START=`grep -m 1 $i ../blastn/$j/$j.bestaln.txt | awk '{print $9}'` 
            END=`grep -m 1 $i ../blastn/$j/$j.bestaln.txt | awk '{print $10}'`
            ORIENTATION=`if [ $START -lt $END ]; then echo "+"; else echo "-"; fi`
            ALIGN="BLAST"
            echo "$ALIGN $i $j $NUMHITS $CONTIG $START $END" >> query/$i.positions.txt
            # grab alignment from samtools
            if [ $ORIENTATION = "-" ]; then ORI="--reverse-complement"; BEG=$START; START=$END; END=$BEG; else ORI=""; fi
            samtools faidx ../exonerate/exonerate_genomes/$j.masked.fa $CONTIG:$START-$END $ORI --mark-strand sign >> query/$i.fna
            for dir in upstream downstream; do echo -e ">$LIFESTYLE B $j $CONTIG:$START-$END\n" >> $dir/$i.fna; done
            LENGTH=`expr $END - $START`
            echo -e "$i\t$j\t$ALIGN\t$ORIENTATION\t$CONTIG\t$UP1\t$UP2\t$START\t$END\t$DOWN1\t$DOWN2\t$LENGTH"
            # add in extra info into fasta header (including sexual/asexual/unknown, hit type and species name)
            LIFESTYLE=`grep -m 1 $j lifestyle.txt  | awk '{print $2}'` 
            sed -i "s/$CONTIG:$START-$END/$LIFESTYLE B $j $CONTIG:$START-$END/g"  query/$i.fna
            for dir in upstream downstream; 
               do echo -e ">$LIFESTYLE B $j $CONTIG:$NEXTSTART-$NEXTEND\n" >> $dir/$i.fna; 
               printf 'N%0.s' {1..99} >> $dir/$i.fna;  
               echo >> $dir/$i.fna;  
            done
            # if there is a second hit on the same contig in the same orientation get that too
            NUMHITS=`grep $i ../blast/$j/$j.bestaln.txt | wc -l | awk '{print $1}'`
            echo "number of blast hits = $NUMHITS"
            if [ $NUMHITS != 1 ]; 
            then
                NEXTSTART=`grep $i ../blast/$j/$j.bestaln.txt | awk '/$CONTIG/{i++}i==2' | awk '{print $2}'`
                NEXTEND= `grep $i ../blast/$j/$j.bestaln.txt | awk '/$CONTIG/{i++}i==2' | awk '{print $2}'`
                NEXTORI=`if [ $START -lt $END ]; then echo "+"; else echo "-"; fi`
                if [ $NEXTORI = $ORIENTATION ] 
                then
                    echo "second blast hit $i $j $CONTIG"
                    echo "$ALIGN $i $j $NUMHITS $CONTIG $NEXTSTART $NEXTEND" >> query/$i.positions.txt
                    samtools faidx ../exonerate/exonerate_genomes/$j.masked.fa $CONTIG:$NEXTSTART-$NEXTEND $ORI --mark-strand sign >> query/$i.fna    
                    sed -i "s/$CONTIG:$NEXTSTART-$NEXTEND/$LIFESTYLE B $j $CONTIG:$NEXTSTART-$NEXTEND/g"  query/$i.fna
                   for dir in upstream downstream; 
                       do echo -e ">$LIFESTYLE B $j $CONTIG:$NEXTSTART-$NEXTEND\n" >> $dir/$i.fna;
                       printf 'N%0.s' {1..99} >> $dir/$i.fna; 
                       echo >> $dir/$i.fna;  
                   done
                fi
            fi
        fi
    done
done



