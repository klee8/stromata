# trim reads including adapters

#!/usr/bin/env bash


# Trim data

# - trimlog records readname and lengths before and after trimming
#ILLUMINACLIP:TruSeq2-PE.fa:2:30:10   <fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>
#SLIDINGWINDOW:5:20   <windowSize>:<requiredQuality>
#LEADING:5         removes low quality bases from the start of a read
#TRAILING:5        removes low quality bases from the end of a read
#MINLEN:40         removes reads shorter than this

mkdir trimmed trimmed/E.elymi trimmed/E.elymi/log
java -jar /home/kate/bin/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 6 -trimlog trimmed/E.elymi/log/Ee-INF1_trim.log /opt/stroma_RNAseq/Ee-Inf1_S14_L006_R1_001.fastq.gz trimmed/E.elymi/Ee-Inf1_S14_L006_R1_001.trim.fastq ILLUMINACLIP:/home/kate/bin/Trimmomatic-0.38/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:40
java -jar /home/kate/bin/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 6 -trimlog trimmed/E.elymi/log/Ee-INF2_trim.log /opt/stroma_RNAseq/Ee-Inf2_S15_L006_R1_001.fastq.gz trimmed/E.elymi/Ee-Inf2_S15_L006_R1_001.trim.fastq ILLUMINACLIP:/home/kate/bin/Trimmomatic-0.38/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:40
java -jar /home/kate/bin/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 6 -trimlog trimmed/E.elymi/log/Ee-INF3_trim.log /opt/stroma_RNAseq/Ee-Inf3_S16_L006_R1_001.fastq.gz trimmed/E.elymi/Ee-Inf3_S16_L006_R1_001.trim.fastq ILLUMINACLIP:/home/kate/bin/Trimmomatic-0.38/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:40
java -jar /home/kate/bin/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 6 -trimlog trimmed/E.elymi/log/Ee-PS1_trim.log /opt/stroma_RNAseq/Ee-PS1_S7_L007_R1_001.fastq.gz trimmed/E.elymi/Ee-PS1_S7_L007_R1_001.trim.fastq ILLUMINACLIP:/home/kate/bin/Trimmomatic-0.38/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:40
java -jar /home/kate/bin/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 6 -trimlog trimmed/E.elymi/log/Ee-PS2_trim.log /opt/stroma_RNAseq/Ee-PS2_S8_L008_R1_001.fastq.gz trimmed/E.elymi/Ee-PS2_S8_L008_R1_001.trim.fastq ILLUMINACLIP:/home/kate/bin/Trimmomatic-0.38/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:40
java -jar /home/kate/bin/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 6 -trimlog trimmed/E.elymi/log/Ee-PS3_trim.log /opt/stroma_RNAseq/Ee-PS3_S13_L005_R1_001_2.fastq.gz trimmed/E.elymi/Ee-PS3_S13_L005_R1_001_2.trim.fastq ILLUMINACLIP:/home/kate/bin/Trimmomatic-0.38/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:40
java -jar /home/kate/bin/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 6 -trimlog trimmed/E.elymi/log/Ee-STR1_trim.log /opt/stroma_RNAseq/Ee-STR1_S10_L004_R1_001.fastq.gz trimmed/E.elymi/Ee-STR1_S10_L004_R1_001.trim.fastq ILLUMINACLIP:/home/kate/bin/Trimmomatic-0.38/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:40
java -jar /home/kate/bin/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 6 -trimlog trimmed/E.elymi/log/Ee-STR2_trim.log /opt/stroma_RNAseq/Ee-STR2_S11_L004_R1_001.fastq.gz trimmed/E.elymi/Ee-STR2_S11_L004_R1_001.trim.fastq ILLUMINACLIP:/home/kate/bin/Trimmomatic-0.38/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:40
java -jar /home/kate/bin/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 6 -trimlog trimmed/E.elymi/log/Ee-STR3_trim.log /opt/stroma_RNAseq/Ee-STR3_S12_L004_R1_001.fastq.gz trimmed/E.elymi/Ee-STR3_S12_L004_R1_001.trim.fastq ILLUMINACLIP:/home/kate/bin/Trimmomatic-0.38/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:40
