#KATE LEE logfile for Epichloe project

###20180321
#### atom and git
set up atom and git repository for Epichloe

####reading
Kallisto does a pseudo alignment which breaks reads into kmers, takes a subset of them to map. Non-mapping kmers are dropped. Tracks kmers mapped and read of origin. Fast and accurate abundance estimates [may not be suitable for this project due to presence of plant material - I think it is likely that fungi and plants will have kmers in common]
Dupont (2015) looked at changes in grass after infection, lots of downregulation of stress response genes etc despite the phenotype being fairly close to WT
####QA
ran SolexaQA++ on Ee-STR1 fastq file
note ran the first one as undefined format, identified as Sanger, implemented sanger from here on
```{bash}
SolexaQA++ analysis /opt/stroma_RNAseq/Ee-PS1_S7_L007_R1_001.fastq.gz -d Ee_SolexaQC/Ee_PS1_S7_L007_R1_001
```
found that reads were of high quality - assumed they have previously been trimmed by David for assemblies (checked with him and they have not been qc'd after sequencing facility trimmed adapters, QC below)

####QC
ran SolexaQA++ on Ee-STR1 fastq file, trimmed with < phred 20 cutoff
```{bash}
SolexaQA++ dynamictrim /opt/stroma_RNAseq/Ee-PS1_S7_L007_R1_001.fastq.gz --phredcutoff 20 --bwa --directory E.elymi --sanger
```

####index genome with STAR
```{bash}
STAR --runMode genomeGenerate --runThreadN 4 --genomeDir genomes/Eel_728 --genomeLoad LoadAndKeep --genomeFastaFiles genomes/Eel_728/Eel_masked.fasta --outFileNamePrefix Eel_728_
```



###20180322
####overnight:
3 PS and 3 STR E.elymi files from the first run were successfully trimmed overnight with SolexaQA++
All files were analysed for initial quality with SolexaQA++

####reading
Alpine models GC bias in RNA-seq data and significantly reduces the number of false positives (especially where there are isoforms with high GC content)
Salmon implements same models as Alpine and provides accurate transcript abudance estimates either from its own light mapping or from BAM genomeFastaFiles

####notes
E6255 == NFE728 see David's email :
"Among our collaborators, there are two commonly used strain-naming conventions. U. Kentyck numbers starting with "E" and Noble Foundation ones starting "NfE". The strain used in the study is called E6255 in the U. Kentucky system and NfE728 in the Noble Foundation system."


####finish QC
run 3 INF E.elymi files through trimming with SolexaQA++


####mapping sample reads to E.elymi
STAR --runMode alignReads \
--runThreadN 8 \
--genomeDir genomes/Eel_728/ \      # NB on the forward slash for end of folder
--readFilesIn ../1_QC/E.elymi/........ \
--readFilesCommand zcat \           # input in .gz files

--sjdbGTFfile GFFfilepath \         # GTF/GFF3 file        
--sjdbGTFfeatureExon exon \         # which feature to pull out (3rd column of gff file)
--sjdbGTFtagExonParentTranscript Parent \    #NB this should be transcript_id, except we didn't put them in gff file*  
--outFileNamePrefix mappedfiles/E.elymi_whatever_the_file_prefix_is \
--outSAMtype BAM SortedByCoordinate \

* David made one canonical gene and transcript per gene, this is numbered in the 'Parent' label. No 'transcript_id' or 'gene_id' features numbered

##### should I try this parameter too?
--outSAMstrandField intronMotif \  # filters out reads with non-canonical introns (David used this)
##### example:
```{bash}
STAR --runMode alignReads --runThreadN 8 --genomeDir genomes/Eel_728/ --readFilesIn ../1_QC/E.elymi/Ee-Inf1_S14_L006_R1_001.fastq.trimmed.gz --readFilesCommand zcat --sjdbGTFfile genomes/Eel_728/Epichloe_elymi.gff3 --sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript Parent --outFileNamePrefix mappedfiles/INF1/E.elymi_Ee-Inf1_S14_L006_R1_001_
```



###20180325
####overnight:
STAR mapped first run E.elymi reads to genomes. Re-running PS1 due to missing gap in command line

####notes
moving forward with Salmon (100s of citations) and DEseq2 (still commonly used)
Had a look at trimmed data with fastqc - adapters still present, need to trim these and re-map

####meeting
meeting with Barry, Daniel and David at 3pm.
Yonathan's project is RNAseq for two mutants which have lost the ability to infect hosts.
Getting sequenced in Novagene. Wants DE analysis done between these and WT.

# Salmon
maps to transcriptome (can make fastq files from bam file)

####Rstudio
set up Rstudio added key and linux libraries to install tidyverse packages within Rstudio (intall.packages() etc.)
see https://cran.r-project.org/bin/linux/ubuntu/README.html
```
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo apt-get install libcurl4-openssl-dev
sudo apt-get install libssl-dev
sudo apt-get install libxml2-dev
```

###DONE
analyse quality of reads with SolexaQA++
trim reads < phred 20
map to Epichloe genomes with STAR


###TO DO
re-do QC with Fastqc and trimmomatic
re-map to genome
quantify abundance with Salmon
DE with DEseq2 and
see if there is a grass genome to map to to remove plant material?
remove bacterial and other genes?

###QUESTIONS
what is masked in the Epichloe genome assembly 'masked' fasta files?
