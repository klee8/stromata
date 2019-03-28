# QC reads using SolexaQA++ v3.1.7.1
# usage info at http://solexaqa.sourceforge.net/
# trimmed reads < phred 20
#################################################

#!usr/bin/bash

# trimed with a cut off of < pred 20
mkdir SolexaQA/ SolexaQA/E.elymi
SolexaQA++ dynamictrim /opt/stroma_RNAseq/Ee-Inf1_S14_L006_R1_001.fastq.gz --phredcutoff 20 --bwa --directory SolexaQA/E.elymi --sanger
SolexaQA++ dynamictrim /opt/stroma_RNAseq/Ee-Inf2_S15_L006_R1_001.fastq.gz --phredcutoff 20 --bwa --directory SolexaQA/E.elymi --sanger
SolexaQA++ dynamictrim /opt/stroma_RNAseq/Ee-Inf3_S16_L006_R1_001.fastq.gz --phredcutoff 20 --bwa --directory SolexaQA/E.elymi --sanger
SolexaQA++ dynamictrim /opt/stroma_RNAseq/Ee-PS1_S7_L007_R1_001.fastq.gz --phredcutoff 20 --bwa --directory SolexaQA/E.elymi --sanger
SolexaQA++ dynamictrim /opt/stroma_RNAseq/Ee-PS2_S8_L008_R1_001.fastq.gz --phredcutoff 20 --bwa --directory SolexaQA/E.elymi --sanger
SolexaQA++ dynamictrim /opt/stroma_RNAseq/Ee-PS3_S13_L005_R1_001_2.fastq.gz --phredcutoff 20 --bwa --directory SolexaQA/E.elymi --sanger
SolexaQA++ dynamictrim /opt/stroma_RNAseq/Ee-STR1_S10_L004_R1_001.fastq.gz --phredcutoff 20 --bwa --directory SolexaQA/E.elymi --sanger
SolexaQA++ dynamictrim /opt/stroma_RNAseq/Ee-STR2_S11_L004_R1_001.fastq.gz --phredcutoff 20 --bwa --directory SolexaQA/E.elymi --sanger
SolexaQA++ dynamictrim /opt/stroma_RNAseq/Ee-STR3_S12_L004_R1_001.fastq.gz --phredcutoff 20 --bwa --directory SolexaQA/E.elymi --sanger
