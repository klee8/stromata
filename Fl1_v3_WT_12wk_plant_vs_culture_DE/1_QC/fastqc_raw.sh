# fastQC of reads prior to trimming

#!/usr/bin/env bash
source /home/kate/.bash_profile
mkdir fastQC_raw
mkdir fastQC_raw/E.festucae_Fl1

zcat ../0_raw_data/run_1/Lane1/SRR7640294.fq.gz  | fastqc stdin --outdir=fastQC_raw/E.festucae_Fl1
zcat ../0_raw_data/run_1/Lane2/SRR7640295.fq.gz  | fastqc stdin --outdir=fastQC_raw/E.festucae_Fl1
zcat ../0_raw_data/run_1/Lane1/SRR7640296.fq.gz  | fastqc stdin --outdir=fastQC_raw/E.festucae_Fl1
zcat ../0_raw_data/run_1/Lane2/SRR7640297.fq.gz  | fastqc stdin --outdir=fastQC_raw/E.festucae_Fl1
zcat ../../Fl1_v3_WT_3day_plant_vs_culture_DE/0_raw_data/run_1/Lane1/SRR7640306.fq.gz  | fastqc stdin --outdir=fastQC_raw/E.festucae_Fl1
zcat ../../Fl1_v3_WT_3day_plant_vs_culture_DE/0_raw_data/run_1/Lane1/SRR7640308.fq.gz  | fastqc stdin --outdir=fastQC_raw/E.festucae_Fl1
zcat ../../Fl1_v3_WT_3day_plant_vs_culture_DE/0_raw_data/run_1/Lane1/SRR7640309.fq.gz  | fastqc stdin --outdir=fastQC_raw/E.festucae_Fl1
