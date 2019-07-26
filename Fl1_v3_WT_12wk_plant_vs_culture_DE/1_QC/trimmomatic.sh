# trim reads including adapters

#!/usr/bin/env bash
source /home/kate/.bash_profile

# Trim data

# - trimlog records readname and lengths before and after trimming
#ILLUMINACLIP:TruSeq2-PE.fa:2:30:10   <fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>
#SLIDINGWINDOW:5:20   <windowSize>:<requiredQuality>
#LEADING:5         removes low quality bases from the start of a read
#TRAILING:5        removes low quality bases from the end of a read
#MINLEN:40         removes reads shorter than this
# note: stopped outputting the trimlog, I don't use it and it takes up a LOT of space.
# note: stdin was overwriting itself in the fastqc_raw/$species output directory - added $run as an extra folder in the path
mkdir trimmed

mkdir trimmed/E.festucae_Fl1 trimmed/E.festucae_Fl1/logs
java -jar /home/kate/bin/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 4 ../0_raw_data/run_1/Lane1/SRR7640294.fq.gz trimmed/E.festucae_Fl1/SRR7640294.trim.fq.gz ILLUMINACLIP:/home/kate/bin/Trimmomatic-0.38/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:40
java -jar /home/kate/bin/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 4 ../0_raw_data/run_1/Lane2/SRR7640295.fq.gz trimmed/E.festucae_Fl1/SRR7640295.trim.fq.gz ILLUMINACLIP:/home/kate/bin/Trimmomatic-0.38/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:40
java -jar /home/kate/bin/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 4 ../0_raw_data/run_1/Lane1/SRR7640296.fq.gz trimmed/E.festucae_Fl1/SRR7640296.trim.fq.gz ILLUMINACLIP:/home/kate/bin/Trimmomatic-0.38/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:40
java -jar /home/kate/bin/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 4 ../0_raw_data/run_1/Lane2/SRR7640297.fq.gz trimmed/E.festucae_Fl1/SRR7640297.trim.fq.gz ILLUMINACLIP:/home/kate/bin/Trimmomatic-0.38/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:40
java -jar /home/kate/bin/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 4 ../../Fl1_v3_WT_3day_plant_vs_culture_DE/0_raw_data/run_1/Lane1/SRR7640306.fq.gz trimmed/E.festucae_Fl1/SRR7640306.trim.fq.gz ILLUMINACLIP:/home/kate/bin/Trimmomatic-0.38/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:40
java -jar /home/kate/bin/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 4 ../../Fl1_v3_WT_3day_plant_vs_culture_DE/0_raw_data/run_1/Lane1/SRR7640308.fq.gz trimmed/E.festucae_Fl1/SRR7640308.trim.fq.gz ILLUMINACLIP:/home/kate/bin/Trimmomatic-0.38/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:40
java -jar /home/kate/bin/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 4 ../../Fl1_v3_WT_3day_plant_vs_culture_DE/0_raw_data/run_1/Lane1/SRR7640309.fq.gz trimmed/E.festucae_Fl1/SRR7640309.trim.fq.gz ILLUMINACLIP:/home/kate/bin/Trimmomatic-0.38/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:40
