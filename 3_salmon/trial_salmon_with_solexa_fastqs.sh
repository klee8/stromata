# Trial run of Salmon with SolexaQA mapped and bedtofastq for fastq files

#!/usr/bin/env bash

# make sure Salmon software is in the $PATH
source /home/kate/.bash_profile

#
mkdir trial_salmon

# Salmon usage
# build an index of the transcriptome
# -k 31 suggested for reads over 75bp
salmon index -t transcriptome/Epichloe_elymi.transcripts.fa -i transcriptome/Epichloe_elymi_trans_index -k 31

# -l U for single-end (unstranded) librarytype
# -r for SE fq file
# --validateMappings makes mapping more sensitive
# --gcBias
salmon quant -i transcriptome/Epichloe_elymi_trans_index -l U -r ../2_star_mapping/SolexaQA_mappedfiles/INF1/E.elymi_Ee-Inf1_S14_L006_R1_001.fastq --validateMappings --gcBias -o trial_salmon/Epichloe_elymi_INF1_quant
salmon quant -i transcriptome/Epichloe_elymi_trans_index -l U -r ../2_star_mapping/SolexaQA_mappedfiles/INF2/E.elymi_Ee-Inf2_S15_L006_R1_001.fastq --validateMappings --gcBias -o trial_salmon/Epichloe_elymi_INF2_quant
salmon quant -i transcriptome/Epichloe_elymi_trans_index -l U -r ../2_star_mapping/SolexaQA_mappedfiles/INF3/E.elymi_Ee-Inf3_S16_L006_R1_001.fastq --validateMappings --gcBias -o trial_salmon/Epichloe_elymi_INF3_quant
salmon quant -i transcriptome/Epichloe_elymi_trans_index -l U -r ../2_star_mapping/SolexaQA_mappedfiles/PS1/E.elymi_Ee-PS1_S7_L007_R1_001.fastq --validateMappings --gcBias -o trial_salmon/Epichloe_elymi_PS1_quant
salmon quant -i transcriptome/Epichloe_elymi_trans_index -l U -r ../2_star_mapping/SolexaQA_mappedfiles/PS2/E.elymi_Ee-PS2_S8_L008_R1_001.fastq --validateMappings --gcBias -o trial_salmon/Epichloe_elymi_PS2_quant
salmon quant -i transcriptome/Epichloe_elymi_trans_index -l U -r ../2_star_mapping/SolexaQA_mappedfiles/PS3/E.elymi_Ee-PS3_S13_L005_R1_001_2.fastq --validateMappings --gcBias -o trial_salmon/Epichloe_elymi_PS3_quant
salmon quant -i transcriptome/Epichloe_elymi_trans_index -l U -r ../2_star_mapping/SolexaQA_mappedfiles/STR1/E.elymi_Ee-STR1_S10_L004_R1_001.fastq --validateMappings --gcBias -o trial_salmon/Epichloe_elymi_STR1_quant
salmon quant -i transcriptome/Epichloe_elymi_trans_index -l U -r ../2_star_mapping/SolexaQA_mappedfiles/STR2/E.elymi_Ee-STR2_S11_L004_R1_001.fastq --validateMappings --gcBias -o trial_salmon/Epichloe_elymi_STR2_quant
salmon quant -i transcriptome/Epichloe_elymi_trans_index -l U -r ../2_star_mapping/SolexaQA_mappedfiles/STR3/E.elymi_Ee-STR3_S12_L004_R1_001.fastq --validateMappings --gcBias -o trial_salmon/Epichloe_elymi_STR3_quant




#../2_star_mapping/SolexaQA_mappedfiles/INF1/E.elymi_Ee-Inf1_S14_L006_R1_001.fastq
#../2_star_mapping/SolexaQA_mappedfiles/INF2/E.elymi_Ee-Inf2_S15_L006_R1_001.fastq
#../2_star_mapping/SolexaQA_mappedfiles/INF3/E.elymi_Ee-Inf3_S16_L006_R1_001.fastq
#../2_star_mapping/SolexaQA_mappedfiles/PS1/E.elymi_Ee-PS1_S7_L007_R1_001.fastq
#../2_star_mapping/SolexaQA_mappedfiles/PS2/E.elymi_Ee-PS2_S8_L008_R1_001.fastq
#../2_star_mapping/SolexaQA_mappedfiles/PS3/E.elymi_Ee-PS3_S13_L005_R1_001_2.fastq
#../2_star_mapping/SolexaQA_mappedfiles/STR1/E.elymi_Ee-STR1_S10_L004_R1_001.fastq
#../2_star_mapping/SolexaQA_mappedfiles/STR2/E.elymi_Ee-STR2_S11_L004_R1_001.fastq
#../2_star_mapping/SolexaQA_mappedfiles/STR3/E.elymi_Ee-STR3_S12_L004_R1_001.fastq
