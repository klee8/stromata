# convert SolexaQA++ qcd and mapped files to fastq for use in salmon with mapping option to play with data while waiting for
# fastqc trimmed and mapped data to come through
#!/usr/bin/env bash

# ensure that bedtools is in the $PATH
source /home/kate/.bash_profile

#bedtools bamtofastq -i INF1/E.elymi_Ee-Inf1_S14_L006_R1_001_Aligned.out.sam -fq INF1/E.elymi_Ee-Inf1_S14_L006_R1_001
bedtools bamtofastq -i INF2/E.elymi_Ee-Inf2_S15_L006_R1_001_Aligned.out.sam -fq INF2/E.elymi_Ee-Inf2_S15_L006_R1_001.fastq
bedtools bamtofastq -i INF3/E.elymi_Ee-Inf3_S16_L006_R1_001_Aligned.out.sam -fq INF3/E.elymi_Ee-Inf3_S16_L006_R1_001.fastq
bedtools bamtofastq -i PS1/E.elymi_Ee-PS1_S7_L007_R1_001_Aligned.out.sam -fq PS1/E.elymi_Ee-PS1_S7_L007_R1_001.fastq
bedtools bamtofastq -i PS2/E.elymi_Ee-PS2_S8_L008_R1_001_Aligned.out.sam -fq PS2/E.elymi_Ee-PS2_S8_L008_R1_001.fastq
bedtools bamtofastq -i PS3/E.elymi_Ee-PS3_S13_L005_R1_001_2_Aligned.out.sam -fq PS3/E.elymi_Ee-PS3_S13_L005_R1_001_2.fastq
bedtools bamtofastq -i STR1/E.elymi_Ee-STR1_S10_L004_R1_001_Aligned.out.sam -fq STR1/E.elymi_Ee-STR1_S10_L004_R1_001.fastq
bedtools bamtofastq -i STR2/E.elymi_Ee-STR2_S11_L004_R1_001_Aligned.out.sam -fq STR2/E.elymi_Ee-STR2_S11_L004_R1_001.fastq
bedtools bamtofastq -i STR3/E.elymi_Ee-STR3_S12_L004_R1_001_Aligned.out.sam -fq STR3/E.elymi_Ee-STR3_S12_L004_R1_001.fastq
