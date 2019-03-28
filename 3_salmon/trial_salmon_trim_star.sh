# Trial run of Salmon with SolexaQA mapped and bedtofastq for fastq files

#!/usr/bin/env bash

# make sure Salmon software is in the $PATH
source /home/kate/.bash_profile

#
mkdir trial_salmon_trim_star

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

#salmon quant -p 8 -i transcriptome/Epichloe_elymi_trans_index -l U -r ../2_star_mapping/mappedfiles/INF1/E.elymi_Ee-Inf1_S14_L006_R1_001.fastq --validateMappings --gcBias -o trial_salmon_trim_star/Epichloe_elymi_INF1_quant
#salmon quant -p 8 -i transcriptome/Epichloe_elymi_trans_index -l U -r ../2_star_mapping/mappedfiles/INF2/E.elymi_Ee-Inf2_S15_L006_R1_001.fastq --validateMappings --gcBias -o trial_salmon_trim_star/Epichloe_elymi_INF2_quant
#salmon quant -p 8 -i transcriptome/Epichloe_elymi_trans_index -l U -r ../2_star_mapping/mappedfiles/PS1/E.elymi_Ee-PS1_S7_L007_R1_001.fastq --validateMappings --gcBias -o trial_salmon_trim_star/Epichloe_elymi_PS1_quant
#salmon quant -p 8 -i transcriptome/Epichloe_elymi_trans_index -l U -r ../2_star_mapping/mappedfiles/PS2/E.elymi_Ee-PS2_S8_L008_R1_001.fastq --validateMappings --gcBias -o trial_salmon_trim_star/Epichloe_elymi_PS2_quant




#########  Alignment OPTION (takes in transcriptome fasta and bam files and does counts)
#salmon quant -p 8  -t transcriptome/Epichloe_elymi.transcripts.fa -l U -a ../2_star_mapping/mappedfiles/INF1/E.elymi_Ee-Inf1_S14_L006_R1_001_Aligned.toTranscriptome.out.bam --gcBias -o trial_salmon_trim_star/Epichloe_elymi_INF1_aln
#salmon quant -p 8  -t transcriptome/Epichloe_elymi.transcripts.fa -l U -a ../2_star_mapping/mappedfiles/INF2/E.elymi_Ee-Inf2_S15_L006_R1_001_Aligned.toTranscriptome.out.bam --gcBias -o trial_salmon_trim_star/Epichloe_elymi_INF2_aln
#salmon quant -p 8  -t transcriptome/Epichloe_elymi.transcripts.fa -l U -a ../2_star_mapping/mappedfiles/PS1/E.elymi_Ee-PS1_S7_L007_R1_001_Aligned.toTranscriptome.out.bam --gcBias -o trial_salmon_trim_star/Epichloe_elymi_PS1_aln
salmon quant -p 8  -t transcriptome/Epichloe_elymi.transcripts.fa -l U -a ../2_star_mapping/mappedfiles/PS2/E.elymi_Ee-PS2_S8_L008_R1_001_Aligned.toTranscriptome.out.bam --gcBias -o trial_salmon_trim_star/Epichloe_elymi_PS2_aln
