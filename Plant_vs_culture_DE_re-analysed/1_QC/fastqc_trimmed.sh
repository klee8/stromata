# fastQC of reads after trimming

#!/usr/bin/env bash
source /home/kate/.bash_profile
mkdir fastQC_trimmed
mkdir fastQC_trimmed/E.festucae_Fl1

fastqc trimmed/E.festucae_Fl1/SRR7640294.trim.fq.gz --outdir=fastQC_trimmed/E.festucae_Fl1
fastqc trimmed/E.festucae_Fl1/SRR7640295.trim.fq.gz --outdir=fastQC_trimmed/E.festucae_Fl1
fastqc trimmed/E.festucae_Fl1/SRR7640296.trim.fq.gz --outdir=fastQC_trimmed/E.festucae_Fl1
fastqc trimmed/E.festucae_Fl1/SRR7640297.trim.fq.gz --outdir=fastQC_trimmed/E.festucae_Fl1
fastqc trimmed/E.festucae_Fl1/SRR7640306.trim.fq.gz --outdir=fastQC_trimmed/E.festucae_Fl1
fastqc trimmed/E.festucae_Fl1/SRR7640308.trim.fq.gz --outdir=fastQC_trimmed/E.festucae_Fl1
fastqc trimmed/E.festucae_Fl1/SRR7640309.trim.fq.gz --outdir=fastQC_trimmed/E.festucae_Fl1
