# map reads to genome with STAR

#!/usr/bin/env bash

# make a genome index
#STAR --runMode genomeGenerate --runThreadN 4 --genomeDir genomes/Eel_728 --genomeLoad LoadAndKeep --genomeFastaFiles genomes/Eel_728/Eel_masked.fasta --outFileNamePrefix Eel_728_
#STAR --runMode genomeGenerate --runThreadN 4 --genomeDir genomes/Efe_2368 --genomeLoad LoadAndKeep --genomeFastaFiles genomes/Eel_728/Eel_masked.fasta --outFileNamePrefix Efe_2368_
#STAR --runMode genomeGenerate --runThreadN 4 --genomeDir genomes/Ety_8 --genomeLoad LoadAndKeep --genomeFastaFiles genomes/Eel_728/Eel_masked.fasta --outFileNamePrefix Ety_8_


# map reads to genome
# note transcripts in the gff3 file were identified as Parent items that end in T1, Parent items that don't end in T1 are genes
# one canonical transcript per gene was annotated in this gff3 file
# here identified transcripts in gff3 by looking at 'Parent' item
# 
#mkdir mappedfiles mappedfiles/INF1 mappedfiles/INF2 mappedfiles/INF3 mappedfiles/PS1 mappedfiles/PS2 mappedfiles/PS3 mappedfiles/STR1 mappedfiles/STR2 mappedfiles/STR3
STAR --runMode alignReads --runThreadN 8 --genomeDir genomes/Eel_728/ --readFilesIn ../1_QC/E.elymi/Ee-Inf1_S14_L006_R1_001.fastq.trimmed.gz --readFilesCommand zcat --sjdbGTFfile genomes/Eel_728/Epichloe_elymi.gff3 --sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript Parent --outFileNamePrefix mappedfiles/INF1/E.elymi_Ee-Inf1_S14_L006_R1_001_
STAR --runMode alignReads --runThreadN 8 --genomeDir genomes/Eel_728/ --readFilesIn ../1_QC/E.elymi/Ee-Inf2_S15_L006_R1_001.fastq.trimmed.gz --readFilesCommand zcat --sjdbGTFfile genomes/Eel_728/Epichloe_elymi.gff3 --sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript Parent --outFileNamePrefix mappedfiles/INF2/E.elymi_Ee-Inf2_S15_L006_R1_001_
STAR --runMode alignReads --runThreadN 8 --genomeDir genomes/Eel_728/ --readFilesIn ../1_QC/E.elymi/Ee-Inf3_S16_L006_R1_001.fastq.trimmed.gz --readFilesCommand zcat --sjdbGTFfile genomes/Eel_728/Epichloe_elymi.gff3 --sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript Parent --outFileNamePrefix mappedfiles/INF3/E.elymi_Ee-Inf3_S16_L006_R1_001_
STAR --runMode alignReads --runThreadN 8 --genomeDir genomes/Eel_728/ --readFilesIn ../1_QC/E.elymi/Ee-PS1_S7_L007_R1_001.fastq.trimmed.gz --readFilesCommand zcat --sjdbGTFfile genomes/Eel_728/Epichloe_elymi.gff3 --sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript Parent --outFileNamePrefix mappedfiles/PS1/E.elymi_Ee-PS1_S7_L007_R1_001_
STAR --runMode alignReads --runThreadN 8 --genomeDir genomes/Eel_728/ --readFilesIn ../1_QC/E.elymi/Ee-PS2_S8_L008_R1_001.fastq.trimmed.gz --readFilesCommand zcat --sjdbGTFfile genomes/Eel_728/Epichloe_elymi.gff3 --sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript Parent --outFileNamePrefix mappedfiles/PS2/E.elymi_Ee-PS2_S8_L008_R1_001_
STAR --runMode alignReads --runThreadN 8 --genomeDir genomes/Eel_728/ --readFilesIn ../1_QC/E.elymi/Ee-PS3_S13_L005_R1_001_2.fastq.trimmed.gz --readFilesCommand zcat --sjdbGTFfile genomes/Eel_728/Epichloe_elymi.gff3 --sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript Parent --outFileNamePrefix mappedfiles/PS3/E.elymi_Ee-PS3_S13_L005_R1_001_2_
STAR --runMode alignReads --runThreadN 8 --genomeDir genomes/Eel_728/ --readFilesIn ../1_QC/E.elymi/Ee-STR1_S10_L004_R1_001.fastq.trimmed.gz --readFilesCommand zcat --sjdbGTFfile genomes/Eel_728/Epichloe_elymi.gff3 --sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript Parent --outFileNamePrefix mappedfiles/STR1/E.elymi_Ee-STR1_S10_L004_R1_001_
STAR --runMode alignReads --runThreadN 8 --genomeDir genomes/Eel_728/ --readFilesIn ../1_QC/E.elymi/Ee-STR2_S11_L004_R1_001.fastq.trimmed.gz --readFilesCommand zcat --sjdbGTFfile genomes/Eel_728/Epichloe_elymi.gff3 --sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript Parent --outFileNamePrefix mappedfiles/STR2/E.elymi_Ee-STR2_S11_L004_R1_001_
STAR --runMode alignReads --runThreadN 8 --genomeDir genomes/Eel_728/ --readFilesIn ../1_QC/E.elymi/Ee-STR3_S12_L004_R1_001.fastq.trimmed.gz --readFilesCommand zcat --sjdbGTFfile genomes/Eel_728/Epichloe_elymi.gff3 --sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript Parent --outFileNamePrefix mappedfiles/STR3/E.elymi_Ee-STR3_S12_L004_R1_001_


# NOTE get a segmentation fault if you try to run them all together on the desktop
#STAR --runMode alignReads --runThreadN 10 --genomeDir genomes/Eel_728/ --readFilesIn  ../1_QC/E.elymi/Ee-Inf1_S14_L006_R1_001.fastq.trimmed.gz, ../1_QC/E.elymi/Ee-Inf2_S15_L006_R1_001.fastq.trimmed.gz, ../1_QC/E.elymi/Ee-Inf3_S16_L006_R1_001.fastq.trimmed.gz, ../1_QC/E.elymi/Ee-PS1_S7_L007_R1_001.fastq.trimmed.gz, ../1_QC/E.elymi/Ee-PS2_S8_L008_R1_001.fastq.trimmed.gz, ../1_QC/E.elymi/Ee-PS3_S13_L005_R1_001_2.fastq.trimmed.gz, ../1_QC/E.elymi/Ee-STR1_S10_L004_R1_001.fastq.trimmed.gz, ../1_QC/E.elymi/Ee-STR2_S11_L004_R1_001.fastq.trimmed.gz, ../1_QC/E.elymi/Ee-STR3_S12_L004_R1_001.fastq.trimmed.gz --readFilesCommand zcat --sjdbGTFfile genomes/Eel_728/Epichloe_elymi.gff3 --sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript Parent --outFileNamePrefix mappedfiles/E.elymi_all_RNA_first_run_
