# trim reads including adapters

#!/usr/bin/env bash


# Trim data

# - trimlog records readname and lengths before and after trimming
#ILLUMINACLIP:TruSeq2-PE.fa:2:40:15   <fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>
#SLIDINGWINDOW:4:20   <windowSize>:<requiredQuality>
#LEADING:2         removes low quality bases from the start of a read
#TRAILING:2        removes low quality bases from the end of a read
#MINLEN:40         removes reads shorter than this

mkdir trimmed
mkdir trimmed/E.elymi
java -jar /home/kate/bin/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 4 -trimlog Ee-INF1_trim.log -basein /opt/stroma_RNAseq/Ee-Inf1_S14_L006_R1_001.fastq.gz -baseout trimmed/E.elymi/Ee-Inf1_S14_L006_R1_001.trim.fastq ILLUMINACLIP:TruSeq2-PE.fa:2:40:15 SLIDINGWINDOW:4:20 LEADING:2 TRAILING:2 MINLEN:40
java -jar /home/kate/bin/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 4 -trimlog Ee-INF2_trim.log -basein /opt/stroma_RNAseq/Ee-Inf2_S15_L006_R1_001.fastq.gz -baseout trimmed/E.elymi/Ee-Inf2_S15_L006_R1_001.trim.fastq ILLUMINACLIP:TruSeq2-PE.fa:2:40:15 SLIDINGWINDOW:4:20 LEADING:2 TRAILING:2 MINLEN:40
java -jar /home/kate/bin/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 4 -trimlog Ee-INF3_trim.log -basein /opt/stroma_RNAseq/Ee-Inf3_S16_L006_R1_001.fastq.gz -baseout trimmed/E.elymi/Ee-Inf3_S16_L006_R1_001.trim.fastq ILLUMINACLIP:TruSeq2-PE.fa:2:40:15 SLIDINGWINDOW:4:20 LEADING:2 TRAILING:2 MINLEN:40
java -jar /home/kate/bin/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 4 -trimlog Ee-PS1_trim.log -basein /opt/stroma_RNAseq/Ee-PS1_S7_L007_R1_001.fastq.gz -baseout trimmed/E.elymi/Ee-PS1_S7_L007_R1_001.trim.fastq ILLUMINACLIP:TruSeq2-PE.fa:2:40:15 SLIDINGWINDOW:4:20 LEADING:2 TRAILING:2 MINLEN:40
java -jar /home/kate/bin/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 4 -trimlog Ee-PS2_trim.log -basein /opt/stroma_RNAseq/Ee-PS2_S8_L008_R1_001.fastq.gz -baseout trimmed/E.elymi/Ee-PS2_S8_L008_R1_001.trim.fastq ILLUMINACLIP:TruSeq2-PE.fa:2:40:15 SLIDINGWINDOW:4:20 LEADING:2 TRAILING:2 MINLEN:40
java -jar /home/kate/bin/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 4 -trimlog Ee-PS3_trim.log -basein /opt/stroma_RNAseq/Ee-PS3_S13_L005_R1_001_2.fastq.gz -baseout trimmed/E.elymi/Ee-PS3_S13_L005_R1_001_2.trim.fastq ILLUMINACLIP:TruSeq2-PE.fa:2:40:15 SLIDINGWINDOW:4:20 LEADING:2 TRAILING:2 MINLEN:40
java -jar /home/kate/bin/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 4 -trimlog Ee-STR1_trim.log -basein /opt/stroma_RNAseq/Ee-STR1_S10_L004_R1_001.fastq.gz -baseout trimmed/E.elymi/Ee-STR1_S10_L004_R1_001.trim.fastq ILLUMINACLIP:TruSeq2-PE.fa:2:40:15 SLIDINGWINDOW:4:20 LEADING:2 TRAILING:2 MINLEN:40
java -jar /home/kate/bin/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 4 -trimlog Ee-STR2_trim.log -basein /opt/stroma_RNAseq/Ee-STR2_S11_L004_R1_001.fastq.gz -baseout trimmed/E.elymi/Ee-STR2_S11_L004_R1_001.trim.fastq ILLUMINACLIP:TruSeq2-PE.fa:2:40:15 SLIDINGWINDOW:4:20 LEADING:2 TRAILING:2 MINLEN:40
java -jar /home/kate/bin/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 4 -trimlog Ee-STR3_trim.log -basein /opt/stroma_RNAseq/Ee-STR3_S12_L004_R1_001.fastq.gz -baseout trimmed/E.elymi/Ee-STR3_S12_L004_R1_001.trim.fastq ILLUMINACLIP:TruSeq2-PE.fa:2:40:15 SLIDINGWINDOW:4:20 LEADING:2 TRAILING:2 MINLEN:40
