# Salmon mapping and read count quantification

#!/usr/bin/env bash

# make sure Salmon software is in the PATH
source /home/kate/.bash_profile

#
mkdir salmon_quant

# Salmon usage
# build an indices of the transcriptomes
# -k 31 suggested for reads over 75bp
salmon index -t transcriptome/E.elymi_NfE728/Epichloe_elymi_NfE728.transcripts.fa -i transcriptome/E.elymi/E.elymi_NfE728.trans_index -k 31
salmon index -t transcriptome/E.festucae_E2368/Epichloe_festucae_E2368.transcripts.fa -i transcriptome/E.festucae/E.festucae_E2368.trans_index -k 31
salmon index -t transcriptome/E.typhina_E8/Epichloe_typhina_E8.transcripts.fa -i transcriptome/E.typhina_E8/E.typhina_E8.trans_index -k 31

##########  QUANT OPTION (quasi-mapping and counting)
# -l U for single-end (unstranded) librarytype
# -r for SE fq file
# --validateMappings makes mapping more sensitive
# --seqBias corrects for random hexamer bias in pcr priming
# --gcBias
# --posBias enables modeling of a position-specific fragment start distribution (bias from position in read)
# --writeUnmappedNames

mkdir salmon_quant/E.elymi_NfE728
salmon quant -p 8 -i transcriptome/E.elymi_NfE728/E.elymi_NfE728.trans_index -l U -r <(zcat ../1_QC/trimmed/E.elymi_NfE728/Ee-Inf1_S1_L001_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.elymi_NfE728/Ee-Inf1_S1_L001_R1_001.quant
salmon quant -p 8 -i transcriptome/E.elymi_NfE728/E.elymi_NfE728.trans_index -l U -r <(zcat ../1_QC/trimmed/E.elymi_NfE728/Ee-Inf1_S1_L002_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.elymi_NfE728/Ee-Inf1_S1_L002_R1_001.quant
salmon quant -p 8 -i transcriptome/E.elymi_NfE728/E.elymi_NfE728.trans_index -l U -r <(zcat ../1_QC/trimmed/E.elymi_NfE728/Ee-Inf2_S2_L001_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.elymi_NfE728/Ee-Inf2_S2_L001_R1_001.quant
salmon quant -p 8 -i transcriptome/E.elymi_NfE728/E.elymi_NfE728.trans_index -l U -r <(zcat ../1_QC/trimmed/E.elymi_NfE728/Ee-Inf2_S2_L002_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.elymi_NfE728/Ee-Inf2_S2_L002_R1_001.quant
salmon quant -p 8 -i transcriptome/E.elymi_NfE728/E.elymi_NfE728.trans_index -l U -r <(zcat ../1_QC/trimmed/E.elymi_NfE728/Ee-Inf3_S16_L006_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.elymi_NfE728/Ee-Inf3_S16_L006_R1_001.quant
salmon quant -p 8 -i transcriptome/E.elymi_NfE728/E.elymi_NfE728.trans_index -l U -r <(zcat ../1_QC/trimmed/E.elymi_NfE728/Ee-Inf3_S3_L001_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.elymi_NfE728/Ee-Inf3_S3_L001_R1_001.quant
salmon quant -p 8 -i transcriptome/E.elymi_NfE728/E.elymi_NfE728.trans_index -l U -r <(zcat ../1_QC/trimmed/E.elymi_NfE728/Ee-Inf3_S3_L002_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.elymi_NfE728/Ee-Inf3_S3_L002_R1_001.quant
salmon quant -p 8 -i transcriptome/E.elymi_NfE728/E.elymi_NfE728.trans_index -l U -r <(zcat ../1_QC/trimmed/E.elymi_NfE728/Ee-PS1_S3_L006_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.elymi_NfE728/Ee-PS1_S3_L006_R1_001.quant
salmon quant -p 8 -i transcriptome/E.elymi_NfE728/E.elymi_NfE728.trans_index -l U -r <(zcat ../1_QC/trimmed/E.elymi_NfE728/Ee-PS1_S3_L007_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.elymi_NfE728/Ee-PS1_S3_L007_R1_001.quant
salmon quant -p 8 -i transcriptome/E.elymi_NfE728/E.elymi_NfE728.trans_index -l U -r <(zcat ../1_QC/trimmed/E.elymi_NfE728/Ee-PS1_S7_L007_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.elymi_NfE728/Ee-PS1_S7_L007_R1_001.quant
salmon quant -p 8 -i transcriptome/E.elymi_NfE728/E.elymi_NfE728.trans_index -l U -r <(zcat ../1_QC/trimmed/E.elymi_NfE728/Ee-PS2_S4_L006_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.elymi_NfE728/Ee-PS2_S4_L006_R1_001.quant
salmon quant -p 8 -i transcriptome/E.elymi_NfE728/E.elymi_NfE728.trans_index -l U -r <(zcat ../1_QC/trimmed/E.elymi_NfE728/Ee-PS2_S4_L007_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.elymi_NfE728/Ee-PS2_S4_L007_R1_001.quant
salmon quant -p 8 -i transcriptome/E.elymi_NfE728/E.elymi_NfE728.trans_index -l U -r <(zcat ../1_QC/trimmed/E.elymi_NfE728/Ee-PS2_S8_L008_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.elymi_NfE728/Ee-PS2_S8_L008_R1_001.quant
salmon quant -p 8 -i transcriptome/E.elymi_NfE728/E.elymi_NfE728.trans_index -l U -r <(zcat ../1_QC/trimmed/E.elymi_NfE728/Ee-PS3_S13_L005_R1_001_2.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.elymi_NfE728/Ee-PS3_S13_L005_R1_001_2.quant
salmon quant -p 8 -i transcriptome/E.elymi_NfE728/E.elymi_NfE728.trans_index -l U -r <(zcat ../1_QC/trimmed/E.elymi_NfE728/Ee-PS3_S5_L006_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.elymi_NfE728/Ee-PS3_S5_L006_R1_001.quant
salmon quant -p 8 -i transcriptome/E.elymi_NfE728/E.elymi_NfE728.trans_index -l U -r <(zcat ../1_QC/trimmed/E.elymi_NfE728/Ee-PS3_S5_L007_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.elymi_NfE728/Ee-PS3_S5_L007_R1_001.quant
salmon quant -p 8 -i transcriptome/E.elymi_NfE728/E.elymi_NfE728.trans_index -l U -r <(zcat ../1_QC/trimmed/E.elymi_NfE728/Ee-STR1_S10_L004_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.elymi_NfE728/Ee-STR1_S10_L004_R1_001.quant
salmon quant -p 8 -i transcriptome/E.elymi_NfE728/E.elymi_NfE728.trans_index -l U -r <(zcat ../1_QC/trimmed/E.elymi_NfE728/Ee-STR2_S11_L004_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.elymi_NfE728/Ee-STR2_S11_L004_R1_001.quant
salmon quant -p 8 -i transcriptome/E.elymi_NfE728/E.elymi_NfE728.trans_index -l U -r <(zcat ../1_QC/trimmed/E.elymi_NfE728/Ee-STR3_S12_L004_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.elymi_NfE728/Ee-STR3_S12_L004_R1_001.quant
mkdir salmon_quant/E.elymi_NfE728_NfE728
salmon quant -p 8 -i transcriptome/E.elymi_NfE728_NfE728/E.elymi_NfE728_NfE728.trans_index -l U -r <(zcat ../1_QC/trimmed/E.elymi_NfE728_NfE728/Ee-Inf1_S14_L006_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.elymi_NfE728_NfE728/Ee-Inf1_S14_L006_R1_001.quant
salmon quant -p 8 -i transcriptome/E.elymi_NfE728_NfE728/E.elymi_NfE728_NfE728.trans_index -l U -r <(zcat ../1_QC/trimmed/E.elymi_NfE728_NfE728/Ee-Inf2_S15_L006_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.elymi_NfE728_NfE728/Ee-Inf2_S15_L006_R1_001.quant
mkdir salmon_quant/E.festucae_E2368
salmon quant -p 8 -i transcriptome/E.festucae_E2368/E.festucae_E2368.trans_index -l U -r <(zcat ../1_QC/trimmed/E.festucae_E2368/Ef-Inf1_S1_L005_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.festucae_E2368/Ef-Inf1_S1_L005_R1_001.quant
salmon quant -p 8 -i transcriptome/E.festucae_E2368/E.festucae_E2368.trans_index -l U -r <(zcat ../1_QC/trimmed/E.festucae_E2368/Ef-Inf1_S7_L003_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.festucae_E2368/Ef-Inf1_S7_L003_R1_001.quant
salmon quant -p 8 -i transcriptome/E.festucae_E2368/E.festucae_E2368.trans_index -l U -r <(zcat ../1_QC/trimmed/E.festucae_E2368/Ef-Inf2_S8_L003_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.festucae_E2368/Ef-Inf2_S8_L003_R1_001.quant
salmon quant -p 8 -i transcriptome/E.festucae_E2368/E.festucae_E2368.trans_index -l U -r <(zcat ../1_QC/trimmed/E.festucae_E2368/Ef-Inf3_S2_L005_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.festucae_E2368/Ef-Inf3_S2_L005_R1_001.quant
salmon quant -p 8 -i transcriptome/E.festucae_E2368/E.festucae_E2368.trans_index -l U -r <(zcat ../1_QC/trimmed/E.festucae_E2368/Ef-Inf3_S9_L003_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.festucae_E2368/Ef-Inf3_S9_L003_R1_001.quant
salmon quant -p 8 -i transcriptome/E.festucae_E2368/E.festucae_E2368.trans_index -l U -r <(zcat ../1_QC/trimmed/E.festucae_E2368/Ef-PS1_S4_L004_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.festucae_E2368/Ef-PS1_S4_L004_R1_001.quant
salmon quant -p 8 -i transcriptome/E.festucae_E2368/E.festucae_E2368.trans_index -l U -r <(zcat ../1_QC/trimmed/E.festucae_E2368/Ef-PS2_S5_L005_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.festucae_E2368/Ef-PS2_S5_L005_R1_001.quant
salmon quant -p 8 -i transcriptome/E.festucae_E2368/E.festucae_E2368.trans_index -l U -r <(zcat ../1_QC/trimmed/E.festucae_E2368/Ef-PS3_S6_L006_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.festucae_E2368/Ef-PS3_S6_L006_R1_001.quant
salmon quant -p 8 -i transcriptome/E.festucae_E2368/E.festucae_E2368.trans_index -l U -r <(zcat ../1_QC/trimmed/E.festucae_E2368/Ef-STR1_S4_L002_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.festucae_E2368/Ef-STR1_S4_L002_R1_001.quant
salmon quant -p 8 -i transcriptome/E.festucae_E2368/E.festucae_E2368.trans_index -l U -r <(zcat ../1_QC/trimmed/E.festucae_E2368/Ef-STR2_S5_L002_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.festucae_E2368/Ef-STR2_S5_L002_R1_001.quant
salmon quant -p 8 -i transcriptome/E.festucae_E2368/E.festucae_E2368.trans_index -l U -r <(zcat ../1_QC/trimmed/E.festucae_E2368/Ef-STR3_S6_L002_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.festucae_E2368/Ef-STR3_S6_L002_R1_001.quant
mkdir salmon_quant/E.typhina_E8
salmon quant -p 8 -i transcriptome/E.typhina_E8/E.typhina_E8.trans_index -l U -r <(zcat ../1_QC/trimmed/E.typhina_E8/Et-PS1_S1_L001_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.typhina_E8/Et-PS1_S1_L001_R1_001.quant
salmon quant -p 8 -i transcriptome/E.typhina_E8/E.typhina_E8.trans_index -l U -r <(zcat ../1_QC/trimmed/E.typhina_E8/Et-PS2_S2_L002_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.typhina_E8/Et-PS2_S2_L002_R1_001.quant
salmon quant -p 8 -i transcriptome/E.typhina_E8/E.typhina_E8.trans_index -l U -r <(zcat ../1_QC/trimmed/E.typhina_E8/Et-PS3_S3_L003_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.typhina_E8/Et-PS3_S3_L003_R1_001.quant
salmon quant -p 8 -i transcriptome/E.typhina_E8/E.typhina_E8.trans_index -l U -r <(zcat ../1_QC/trimmed/E.typhina_E8/Et-STR1_S1_L001_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.typhina_E8/Et-STR1_S1_L001_R1_001.quant
salmon quant -p 8 -i transcriptome/E.typhina_E8/E.typhina_E8.trans_index -l U -r <(zcat ../1_QC/trimmed/E.typhina_E8/Et-STR2_S2_L001_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.typhina_E8/Et-STR2_S2_L001_R1_001.quant
salmon quant -p 8 -i transcriptome/E.typhina_E8/E.typhina_E8.trans_index -l U -r <(zcat ../1_QC/trimmed/E.typhina_E8/Et-STR3_S3_L001_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/E.typhina_E8/Et-STR3_S3_L001_R1_001.quant
mkdir salmon_quant/UN
salmon quant -p 8 -i transcriptome/UN/UN.trans_index -l U -r <(zcat ../1_QC/trimmed/UN/Undetermined_S0_L001_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/UN/Undetermined_S0_L001_R1_001.quant
salmon quant -p 8 -i transcriptome/UN/UN.trans_index -l U -r <(zcat ../1_QC/trimmed/UN/Undetermined_S0_L002_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/UN/Undetermined_S0_L002_R1_001.quant
salmon quant -p 8 -i transcriptome/UN/UN.trans_index -l U -r <(zcat ../1_QC/trimmed/UN/Undetermined_S0_L005_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/UN/Undetermined_S0_L005_R1_001.quant
salmon quant -p 8 -i transcriptome/UN/UN.trans_index -l U -r <(zcat ../1_QC/trimmed/UN/Undetermined_S0_L006_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/UN/Undetermined_S0_L006_R1_001.quant
salmon quant -p 8 -i transcriptome/UN/UN.trans_index -l U -r <(zcat ../1_QC/trimmed/UN/Undetermined_S0_L007_R1_001.trim.fastq) --validateMappings --seqBias --gcBias --posBias --writeUnmappedNames -o salmon_quant/UN/Undetermined_S0_L007_R1_001.quant
