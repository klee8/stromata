# Trial run of Salmon with SolexaQA mapped and bedtofastq for fastq files

#!/usr/bin/env bash

# make sure Salmon software is in the $PATH
source /home/kate/.bash_profile

#
mkdir salmon_quant

# Salmon usage
# build an index of the transcriptome
# -k 31 suggested for reads over 75bp
#salmon index -t transcriptome/Epichloe_elymi.transcripts.fa -i transcriptome/Epichloe_elymi_trans_index -k 31


##########  QUANT OPTION (quasi-mapping and counting)
# -l U for single-end (unstranded) librarytype
# -r for SE fq file
# --validateMappings makes mapping more sensitive
# --gcBias
# --useVBOpt  variational Bayesian EM algorithm (generally slightly more accurate and can be faster)
# --biasCorrect corrects for random hexamer bias in pcr priming
# --useFSPD enables modeling of a position-specific fragment start distribution (bias from position in read) --useVBOpt --biasCorrect --useFSPD
salmon quant -p 8 -i transcriptome/Epichloe_elymi_trans_index -l U -r <(gzcat ../1_QC/trimmed/E.elymi/Ee-Inf1_S14_L006_R1_001.trim.fastq) --validateMappings --gcBias --useVBOpt --biasCorrect --useFSPD -o salmon_quant/Epichloe_elymi_INF1_quant
salmon quant -p 8 -i transcriptome/Epichloe_elymi_trans_index -l U -r <(gzcat ../1_QC/trimmed/E.elymi/Ee-Inf2_S15_L006_R1_001.trim.fastq.gz) --validateMappings --gcBias --useVBOpt --biasCorrect --useFSPD -o salmon_quant/Epichloe_elymi_INF2_quant
salmon quant -p 8 -i transcriptome/Epichloe_elymi_trans_index -l U -r <(gzcat ../1_QC/trimmed/E.elymi/Ee-Inf3_S16_L006_R1_001.trim.fastq.gz) --validateMappings --gcBias --useVBOpt --biasCorrect --useFSPD -o salmon_quant/Epichloe_elymi_PS1_quant
salmon quant -p 8 -i transcriptome/Epichloe_elymi_trans_index -l U -r <(gzcat ../1_QC/trimmed/E.elymi/Ee-PS1_S7_L007_R1_001.trim.fastq.gz) --validateMappings --gcBias --useVBOpt --biasCorrect --useFSPD -o salmon_quant/Epichloe_elymi_PS1_quant
salmon quant -p 8 -i transcriptome/Epichloe_elymi_trans_index -l U -r <(gzcat ../1_QC/trimmed/E.elymi/Ee-PS2_S8_L008_R1_001.trim.fastq.gz) --validateMappings --gcBias --useVBOpt --biasCorrect --useFSPD -o salmon_quant/Epichloe_elymi_PS2_quant
salmon quant -p 8 -i transcriptome/Epichloe_elymi_trans_index -l U -r <(gzcat ../1_QC/trimmed/E.elymi/Ee-PS3_S13_L005_R1_001_2.trim.fastq.gz) --validateMappings --gcBias --useVBOpt --biasCorrect --useFSPD -o salmon_quant/Epichloe_elymi_PS3_quant
salmon quant -p 8 -i transcriptome/Epichloe_elymi_trans_index -l U -r <(gzcat ../1_QC/trimmed/E.elymi/Ee-STR1_S10_L004_R1_001.trim.fastq.gz) --validateMappings --gcBias --useVBOpt --biasCorrect --useFSPD -o salmon_quant/Epichloe_elymi_STR1_quant
salmon quant -p 8 -i transcriptome/Epichloe_elymi_trans_index -l U -r <(gzcat ../1_QC/trimmed/E.elymi/Ee-STR2_S11_L004_R1_001.trim.fastq.gz) --validateMappings --gcBias --useVBOpt --biasCorrect --useFSPD -o salmon_quant/Epichloe_elymi_STR2_quant
salmon quant -p 8 -i transcriptome/Epichloe_elymi_trans_index -l U -r <(gzcat ../1_QC/trimmed/E.elymi/Ee-STR3_S12_L004_R1_001.trim.fastq.gz) --validateMappings --gcBias --useVBOpt --biasCorrect --useFSPD -o salmon_quant/Epichloe_elymi_STR3_quant
