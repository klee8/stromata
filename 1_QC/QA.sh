# Analyse read quality using SolexaQA++ v3.1.7.1
# usage info at http://solexaqa.sourceforge.net/
#################################################
#!/usr/bin/env bash

# make directories for QA output
mkdir SolexaQA/E.elymi/Ee-Inf1_S14_L006_R1_001
mkdir SolexaQA/E.elymi/Ee-Inf2_S15_L006_R1_001
mkdir SolexaQA/E.elymi/Ee-Inf3_S16_L006_R1_001
mkdir SolexaQA/E.elymi/Ee-PS1_S7_L007_R1_001
mkdir SolexaQA/E.elymi/Ee-PS2_S8_L008_R1_001
mkdir SolexaQA/E.elymi/Ee-PS3_S13_L005_R1_001_2
mkdir SolexaQA/E.elymi/Ee-STR1_S10_L004_R1_001
mkdir SolexaQA/E.elymi/Ee-STR2_S11_L004_R1_001
mkdir SolexaQA/E.elymi/Ee-STR3_S12_L004_R1_001

##### ran basic analysis to see read quality
# inflorescence samples
SolexaQA++ analysis /opt/stroma_RNAseq/Ee-Inf1_S14_L006_R1_001.fastq.gz --sanger -d SolexaQA/E.elymi/Ee-Inf1_S14_L006_R1_001
SolexaQA++ analysis /opt/stroma_RNAseq/Ee-Inf2_S15_L006_R1_001.fastq.gz --sanger -d SolexaQA/E.elymi/Ee-Inf2_S15_L006_R1_001
SolexaQA++ analysis /opt/stroma_RNAseq/Ee-Inf3_S16_L006_R1_001.fastq.gz --sanger -d SolexaQA/E.elymi/Ee-Inf3_S16_L006_R1_001
# pseudostem samples
SolexaQA++ analysis /opt/stroma_RNAseq/Ee-PS1_S7_L007_R1_001.fastq.gz --sanger -d SolexaQA/E.elymi/Ee-PS1_S7_L007_R1_001
SolexaQA++ analysis /opt/stroma_RNAseq/Ee-PS2_S8_L008_R1_001.fastq.gz --sanger -d SolexaQA/E.elymi/Ee-PS2_S8_L008_R1_001
SolexaQA++ analysis /opt/stroma_RNAseq/Ee-PS3_S13_L005_R1_001_2.fastq.gz --sanger -d SolexaQA/E.elymi/Ee-PS3_S13_L005_R1_001_2
# stroma samples
SolexaQA++ analysis /opt/stroma_RNAseq/Ee-STR1_S10_L004_R1_001.fastq.gz --sanger -d SolexaQA/E.elymi/Ee-STR1_S10_L004_R1_001
SolexaQA++ analysis /opt/stroma_RNAseq/Ee-STR2_S11_L004_R1_001.fastq.gz --sanger -d SolexaQA/E.elymi/Ee-STR2_S11_L004_R1_001
SolexaQA++ analysis /opt/stroma_RNAseq/Ee-STR3_S12_L004_R1_001.fastq.gz --sanger -d SolexaQA/E.elymi/Ee-STR3_S12_L004_R1_001
