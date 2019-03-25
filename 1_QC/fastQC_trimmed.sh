# fastQC of reads prior to trimming

#!/usr/bin/env bash
mkdir fastQC_trimmed
mkdir fastQC_trimmed/E.elymi
source /home/kate/.bash_profile
fastqc trimmed/E.elymi/Ee-Inf1_S14_L006_R1_001.trim.fastq --outdir=fastQC_trimmed/E.elymi
fastqc trimmed/E.elymi/Ee-Inf2_S15_L006_R1_001.trim.fastq --outdir=fastQC_trimmed/E.elymi
fastqc trimmed/E.elymi/Ee-Inf3_S16_L006_R1_001.trim.fastq --outdir=fastQC_trimmed/E.elymi
fastqc trimmed/E.elymi/Ee-PS1_S7_L007_R1_001.trim.fastq --outdir=fastQC_trimmed/E.elymi
fastqc trimmed/E.elymi/Ee-PS2_S8_L008_R1_001.trim.fastq --outdir=fastQC_trimmed/E.elymi
fastqc trimmed/E.elymi/Ee-PS3_S13_L005_R1_001_2.trim.fastq --outdir=fastQC_trimmed/E.elymi
fastqc trimmed/E.elymi/Ee-STR1_S10_L004_R1_001.trim.fastq --outdir=fastQC_trimmed/E.elymi
fastqc trimmed/E.elymi/Ee-STR2_S11_L004_R1_001.trim.fastq --outdir=fastQC_trimmed/E.elymi
fastqc trimmed/E.elymi/Ee-STR3_S12_L004_R1_001.trim.fastq --outdir=fastQC_trimmed/E.elymi
