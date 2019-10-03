# Salmon mapping and read count quantification
# make sure Salmon software is in the PATH

#!/usr/bin/env bash
source /home/kate/.bash_profile

#
mkdir salmon_quant

# Salmon usage
# build an indices of the transcriptomes
# -k 31 suggested for reads over 75bp
salmon index -t transcriptome/E.elymi_NfE728/Epichloe_elymi_NfE728.transcripts.fa -i transcriptome/E.elymi/E.elymi_NfE728.trans_index -k 31



# grab headers from transriptome files for use in mapping transcripts to genes later
for i in transcriptome/*; do grep '>' $i/*.fa > $i/transcriptome.headers.txt; sed -i s'/>//g' $i/transcriptome.headers.txt; done

##########  QUANT OPTION (quasi-mapping and counting)
# -l U for single-end (unstranded) librarytype
# -r for SE fq file
# --validateMappings makes mapping more sensitive
# --seqBias corrects for random hexamer bias in pcr priming
# --gcBias
# --posBias enables modeling of a position-specific fragment start distribution (bias from position in read)
# --writeUnmappedNames   # this takes up loads of space - skip it

mkdir salmon_quant/E.festucae_Fl1_v3
salmon quant -p 8 -i transcriptome/E.festucae_Fl1_v3/E.festucae_Fl1_v3.trans_index -l U -r <(zcat ../1_QC/trimmed/E.festucae_Fl1_v3/SRR7640294.trim.fastq) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.festucae_Fl1_v3/SRR7640294
salmon quant -p 8 -i transcriptome/E.festucae_Fl1_v3/E.festucae_Fl1_v3.trans_index -l U -r <(zcat ../1_QC/trimmed/E.festucae_Fl1_v3/SRR7640295.trim.fastq) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.festucae_Fl1_v3/SRR7640295
salmon quant -p 8 -i transcriptome/E.festucae_Fl1_v3/E.festucae_Fl1_v3.trans_index -l U -r <(zcat ../1_QC/trimmed/E.festucae_Fl1_v3/SRR7640296.trim.fastq) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.festucae_Fl1_v3/SRR7640296
salmon quant -p 8 -i transcriptome/E.festucae_Fl1_v3/E.festucae_Fl1_v3.trans_index -l U -r <(zcat ../1_QC/trimmed/E.festucae_Fl1_v3/SRR7640297.trim.fastq) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.festucae_Fl1_v3/SRR7640297
salmon quant -p 8 -i transcriptome/E.festucae_Fl1_v3/E.festucae_Fl1_v3.trans_index -l U -r <(zcat ../1_QC/trimmed/E.festucae_Fl1_v3/SRR7640306.trim.fastq) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.festucae_Fl1_v3/SRR7640306
salmon quant -p 8 -i transcriptome/E.festucae_Fl1_v3/E.festucae_Fl1_v3.trans_index -l U -r <(zcat ../1_QC/trimmed/E.festucae_Fl1_v3/SRR7640308.trim.fastq) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.festucae_Fl1_v3/SRR7640308
salmon quant -p 8 -i transcriptome/E.festucae_Fl1_v3/E.festucae_Fl1_v3.trans_index -l U -r <(zcat ../1_QC/trimmed/E.festucae_Fl1_v3/SRR7640309.trim.fastq) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/E.festucae_Fl1_v3/SRR7640309
