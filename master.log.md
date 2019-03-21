#KATE LEE logfile for Epichloe project

###20180321
#### atom and git
set up atom and git repository for Epichloe

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

####mapping
```{bash}
STAR --runMode alignReads \
--runThreadN 4 \
--genomeDir genomes/Eel_728 \
--outSAMtype BAM SortedByCoordinate \

```


###TO DO
trim reads < phred 15
map to Epichloe genomes with STAR

###QUESTIONS
what is masked in the Epichloe genome assembly 'masked' fasta files?
