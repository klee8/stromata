# QC reads using SolexaQA++ v3.1.7.1
# usage info at http://solexaqa.sourceforge.net/
#################################################

#!usr/bin/bash


# ran basic analysis to see read quality
SolexaQA++ analysis /opt/stroma_RNAseq/Ee-PS1_S7_L007_R1_001.fastq.gz -d Ee_SolexaQC/Ee_PS1_S7_L007_R1_001


# File looks like it has already been trimmed, skip this step and go straight to mapping
