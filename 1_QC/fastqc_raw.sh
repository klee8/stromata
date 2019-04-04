# fastQC of reads prior to trimming

#!/usr/bin/env bash
source /home/kate/.bash_profile
mkdir fastQC_raw
mkdir fastQC_raw/E.elymi_NfE728
zcat /opt/stroma_RNAseq/re_runs/Lane1/Ee-Inf1_S1_L001_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.elymi_NfE728
zcat /opt/stroma_RNAseq/re_runs/Lane2/Ee-Inf1_S1_L002_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.elymi_NfE728
zcat /opt/stroma_RNAseq/re_runs/Lane1/Ee-Inf2_S2_L001_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.elymi_NfE728
zcat /opt/stroma_RNAseq/re_runs/Lane2/Ee-Inf2_S2_L002_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.elymi_NfE728
zcat /opt/stroma_RNAseq/first_run/Lane6/Ee-Inf3_S16_L006_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.elymi_NfE728
zcat /opt/stroma_RNAseq/re_runs/Lane1/Ee-Inf3_S3_L001_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.elymi_NfE728
zcat /opt/stroma_RNAseq/re_runs/Lane2/Ee-Inf3_S3_L002_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.elymi_NfE728
zcat /opt/stroma_RNAseq/re_runs/Lane6/Ee-PS1_S3_L006_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.elymi_NfE728
zcat /opt/stroma_RNAseq/re_runs/Lane7/Ee-PS1_S3_L007_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.elymi_NfE728
zcat /opt/stroma_RNAseq/first_run/Lane7/Ee-PS1_S7_L007_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.elymi_NfE728
zcat /opt/stroma_RNAseq/re_runs/Lane6/Ee-PS2_S4_L006_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.elymi_NfE728
zcat /opt/stroma_RNAseq/re_runs/Lane7/Ee-PS2_S4_L007_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.elymi_NfE728
zcat /opt/stroma_RNAseq/first_run/Lane8/Ee-PS2_S8_L008_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.elymi_NfE728
zcat /opt/stroma_RNAseq/first_run/Lane5/Ee-PS3_S13_L005_R1_001_2.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.elymi_NfE728
zcat /opt/stroma_RNAseq/re_runs/Lane6/Ee-PS3_S5_L006_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.elymi_NfE728
zcat /opt/stroma_RNAseq/re_runs/Lane7/Ee-PS3_S5_L007_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.elymi_NfE728
zcat /opt/stroma_RNAseq/first_run/Lane4/Ee-STR1_S10_L004_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.elymi_NfE728
zcat /opt/stroma_RNAseq/first_run/Lane4/Ee-STR2_S11_L004_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.elymi_NfE728
zcat /opt/stroma_RNAseq/first_run/Lane4/Ee-STR3_S12_L004_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.elymi_NfE728
mkdir fastQC_raw/E.elymi_NfE728_NfE728
zcat /opt/stroma_RNAseq/first_run/Lane6/Ee-Inf1_S14_L006_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.elymi_NfE728_NfE728
zcat /opt/stroma_RNAseq/first_run/Lane6/Ee-Inf2_S15_L006_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.elymi_NfE728_NfE728
mkdir fastQC_raw/E.festucae_E2368
zcat /opt/stroma_RNAseq/re_runs/Lane5/Ef-Inf1_S1_L005_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.festucae_E2368
zcat /opt/stroma_RNAseq/first_run/Lane3/Ef-Inf1_S7_L003_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.festucae_E2368
zcat /opt/stroma_RNAseq/first_run/Lane3/Ef-Inf2_S8_L003_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.festucae_E2368
zcat /opt/stroma_RNAseq/re_runs/Lane5/Ef-Inf3_S2_L005_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.festucae_E2368
zcat /opt/stroma_RNAseq/first_run/Lane3/Ef-Inf3_S9_L003_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.festucae_E2368
zcat /opt/stroma_RNAseq/first_run/Lane4/Ef-PS1_S4_L004_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.festucae_E2368
zcat /opt/stroma_RNAseq/first_run/Lane5/Ef-PS2_S5_L005_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.festucae_E2368
zcat /opt/stroma_RNAseq/first_run/Lane6/Ef-PS3_S6_L006_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.festucae_E2368
zcat /opt/stroma_RNAseq/first_run/Lane2/Ef-STR1_S4_L002_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.festucae_E2368
zcat /opt/stroma_RNAseq/first_run/Lane2/Ef-STR2_S5_L002_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.festucae_E2368
zcat /opt/stroma_RNAseq/first_run/Lane2/Ef-STR3_S6_L002_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.festucae_E2368
mkdir fastQC_raw/E.typhina_E8
zcat /opt/stroma_RNAseq/first_run/Lane1/Et-PS1_S1_L001_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.typhina_E8
zcat /opt/stroma_RNAseq/first_run/Lane2/Et-PS2_S2_L002_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.typhina_E8
zcat /opt/stroma_RNAseq/first_run/Lane3/Et-PS3_S3_L003_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.typhina_E8
zcat /opt/stroma_RNAseq/first_run/Lane1/Et-STR1_S1_L001_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.typhina_E8
zcat /opt/stroma_RNAseq/first_run/Lane1/Et-STR2_S2_L001_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.typhina_E8
zcat /opt/stroma_RNAseq/first_run/Lane1/Et-STR3_S3_L001_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/E.typhina_E8
mkdir fastQC_raw/UN
zcat /opt/stroma_RNAseq/re_runs/Lane1/Undetermined_S0_L001_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/UN
zcat /opt/stroma_RNAseq/re_runs/Lane2/Undetermined_S0_L002_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/UN
zcat /opt/stroma_RNAseq/re_runs/Lane5/Undetermined_S0_L005_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/UN
zcat /opt/stroma_RNAseq/re_runs/Lane6/Undetermined_S0_L006_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/UN
zcat /opt/stroma_RNAseq/re_runs/Lane7/Undetermined_S0_L007_R1_001.fastq.gz  | fastqc stdin --outdir=fastQC_raw/UN
