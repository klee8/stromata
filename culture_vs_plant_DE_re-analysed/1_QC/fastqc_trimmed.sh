# fastQC of reads after trimming

#!/usr/bin/env bash
source /home/kate/.bash_profile
mkdir fastQC_trimmed
mkdir fastQC_trimmed/E.festucae_Fl1_v3
fastqc trimmed/E.festucae_Fl1_v3/SRR7640294.trim.fastq.gz --outdir=fastQC_trimmed/E.festucae_Fl1_v3
fastqc trimmed/E.festucae_Fl1_v3/SRR7640295.trim.fastq.gz --outdir=fastQC_trimmed/E.festucae_Fl1_v3
fastqc trimmed/E.festucae_Fl1_v3/SRR7640296.trim.fastq.gz --outdir=fastQC_trimmed/E.festucae_Fl1_v3
fastqc trimmed/E.festucae_Fl1_v3/SRR7640297.trim.fastq.gz --outdir=fastQC_trimmed/E.festucae_Fl1_v3
fastqc trimmed/E.festucae_Fl1_v3/SRR7640306.trim.fastq.gz --outdir=fastQC_trimmed/E.festucae_Fl1_v3
fastqc trimmed/E.festucae_Fl1_v3/SRR7640308.trim.fastq.gz --outdir=fastQC_trimmed/E.festucae_Fl1_v3
fastqc trimmed/E.festucae_Fl1_v3/SRR7640309.trim.fastq.gz --outdir=fastQC_trimmed/E.festucae_Fl1_v3
