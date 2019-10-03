# fastQC of reads prior to trimming

#!/usr/bin/env bash
source /home/kate/.bash_profile
mkdir fastQC_raw
mkdir fastQC_raw/E.festucae_Fl1_v3
zcat 0_raw_data/run_1/Lane1/SRR7640294.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.festucae_Fl1_v3/SRR7640294
zcat 0_raw_data/run_1/Lane2/SRR7640295.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.festucae_Fl1_v3/SRR7640295
zcat 0_raw_data/run_1/Lane1/SRR7640296.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.festucae_Fl1_v3/SRR7640296
zcat 0_raw_data/run_1/Lane2/SRR7640297.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.festucae_Fl1_v3/SRR7640297
zcat 0_raw_data/run_1/Lane1/SRR7640306.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.festucae_Fl1_v3/SRR7640306
zcat 0_raw_data/run_1/Lane1/SRR7640308.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.festucae_Fl1_v3/SRR7640308
zcat 0_raw_data/run_1/Lane1/SRR7640309.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.festucae_Fl1_v3/SRR7640309
