# Salmon mapping and read count quantification
# make sure Salmon software is in the PATH

#!/usr/bin/env bash
source /home/kate/.bash_profile

#
mkdir salmon_quant

# Salmon usage
# build an indices of the transcriptomes
# -k 31 suggested for reads over 75bp
cp transcriptome/L.perenne/best_candidates.eclipsed_orfs_removed.cds.cdhit95.minusefm3.fas transcriptome/L.perenne/L.perenne.ESTs.fas
salmon index -t transcriptome/L.perenne/L.perenne.ESTs.fas -i transcriptome/L.perenne/L.perenne.ESTs.index -k 31


# grab headers from transriptome files for use in mapping transcripts to genes later
grep '>' transcriptome/L.perenne/L.perenne.ESTs.fas > transcriptome/L.perenne/temp.txt
cut -d ' ' -f 1,2 transcriptome/L.perenne/temp.txt > transcriptome/L.perenne/EST.headers.txt
sed -i s'/>//g' transcriptome/L.perenne/EST.headers.txt
rm transcriptome/L.perenne/temp.txt



##########  QUANT OPTION (quasi-mapping and counting)
# -l U for single-end (unstranded) librarytype
# -r for SE fq file
# --validateMappings makes mapping more sensitive
# --seqBias corrects for random hexamer bias in pcr priming
# --gcBias
# --posBias enables modeling of a position-specific fragment start distribution (bias from position in read)
# --writeUnmappedNames   # this takes up loads of space - skip it

mkdir salmon_quant/E.elymi_NfE728
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.elymi_NfE728/Ee-Inf1_S14_L006_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.elymi_NfE728/Ee-Inf1_S14_L006_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.elymi_NfE728/Ee-Inf1_S1_L001_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.elymi_NfE728/Ee-Inf1_S1_L001_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.elymi_NfE728/Ee-Inf1_S1_L002_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.elymi_NfE728/Ee-Inf1_S1_L002_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.elymi_NfE728/Ee-Inf2_S15_L006_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.elymi_NfE728/Ee-Inf2_S15_L006_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.elymi_NfE728/Ee-Inf2_S2_L001_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.elymi_NfE728/Ee-Inf2_S2_L001_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.elymi_NfE728/Ee-Inf2_S2_L002_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.elymi_NfE728/Ee-Inf2_S2_L002_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.elymi_NfE728/Ee-Inf3_S16_L006_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.elymi_NfE728/Ee-Inf3_S16_L006_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.elymi_NfE728/Ee-Inf3_S3_L001_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.elymi_NfE728/Ee-Inf3_S3_L001_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.elymi_NfE728/Ee-Inf3_S3_L002_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.elymi_NfE728/Ee-Inf3_S3_L002_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.elymi_NfE728/Ee-PS1_S3_L006_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.elymi_NfE728/Ee-PS1_S3_L006_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.elymi_NfE728/Ee-PS1_S3_L007_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.elymi_NfE728/Ee-PS1_S3_L007_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.elymi_NfE728/Ee-PS1_S7_L007_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.elymi_NfE728/Ee-PS1_S7_L007_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.elymi_NfE728/Ee-PS2_S4_L006_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.elymi_NfE728/Ee-PS2_S4_L006_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.elymi_NfE728/Ee-PS2_S4_L007_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.elymi_NfE728/Ee-PS2_S4_L007_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.elymi_NfE728/Ee-PS2_S8_L008_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.elymi_NfE728/Ee-PS2_S8_L008_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.elymi_NfE728/Ee-PS3_S13_L005_R1_001_2.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.elymi_NfE728/Ee-PS3_S13_L005_R1_001_2
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.elymi_NfE728/Ee-PS3_S5_L006_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.elymi_NfE728/Ee-PS3_S5_L006_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.elymi_NfE728/Ee-PS3_S5_L007_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.elymi_NfE728/Ee-PS3_S5_L007_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.elymi_NfE728/Ee-STR1_S10_L004_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.elymi_NfE728/Ee-STR1_S10_L004_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.elymi_NfE728/Ee-STR2_S11_L004_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.elymi_NfE728/Ee-STR2_S11_L004_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.elymi_NfE728/Ee-STR3_S12_L004_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.elymi_NfE728/Ee-STR3_S12_L004_R1_001
mkdir salmon_quant/E.festucae_E2368
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.festucae_E2368/Ef-Inf1_S1_L005_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.festucae_E2368/Ef-Inf1_S1_L005_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.festucae_E2368/Ef-Inf1_S7_L003_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.festucae_E2368/Ef-Inf1_S7_L003_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.festucae_E2368/Ef-Inf2_S8_L003_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.festucae_E2368/Ef-Inf2_S8_L003_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.festucae_E2368/Ef-Inf3_S2_L005_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.festucae_E2368/Ef-Inf3_S2_L005_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.festucae_E2368/Ef-Inf3_S9_L003_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.festucae_E2368/Ef-Inf3_S9_L003_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.festucae_E2368/Ef-PS1_S4_L004_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.festucae_E2368/Ef-PS1_S4_L004_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.festucae_E2368/Ef-PS2_S5_L005_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.festucae_E2368/Ef-PS2_S5_L005_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.festucae_E2368/Ef-PS3_S6_L006_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.festucae_E2368/Ef-PS3_S6_L006_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.festucae_E2368/Ef-STR1_S4_L002_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.festucae_E2368/Ef-STR1_S4_L002_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.festucae_E2368/Ef-STR2_S5_L002_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.festucae_E2368/Ef-STR2_S5_L002_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.festucae_E2368/Ef-STR3_S6_L002_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.festucae_E2368/Ef-STR3_S6_L002_R1_001
mkdir salmon_quant/E.typhina_E8
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.typhina_E8/Et-PS1_S1_L001_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.typhina_E8/Et-PS1_S1_L001_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.typhina_E8/Et-PS2_S2_L002_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.typhina_E8/Et-PS2_S2_L002_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.typhina_E8/Et-PS3_S3_L003_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.typhina_E8/Et-PS3_S3_L003_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.typhina_E8/Et-STR1_S1_L001_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.typhina_E8/Et-STR1_S1_L001_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.typhina_E8/Et-STR2_S2_L001_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.typhina_E8/Et-STR2_S2_L001_R1_001
salmon quant -p 8 -i transcriptome/L.perenne/L.perenne.ESTs.index -l U -r <(zcat ../../stromata_DE/1_QC/trimmed/E.typhina_E8/Et-STR3_S3_L001_R1_001.trim.fastq.gz) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.typhina_E8/Et-STR3_S3_L001_R1_001
