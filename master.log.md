#KATE LEE logfile for Epichloe project

###20190321 Thursday
#### atom and git
set up atom and git repository for Epichloe

####reading
Kallisto does a pseudo alignment which breaks reads into kmers, takes a subset of them to map. Non-mapping kmers are dropped. Tracks kmers mapped and read of origin. Fast and accurate abundance estimates [may not be suitable for this project due to presence of plant material - I think it is likely that fungi and plants will have kmers in common]
Dupont (2015) looked at changes in grass after infection, lots of downregulation of stress response genes etc despite the phenotype being fairly close to WT
####QA
ran SolexaQA++ on Ee-STR1 fastq file
note ran the first one as undefined format, identified as Sanger, implemented sanger from here on
```{bash}
SolexaQA++ analysis /opt/stroma_RNAseq/Ee-PS1_S7_L007_R1_001.fastq.gz -d Ee_SolexaQC/Ee_PS1_S7_L007_R1_001
```
found that reads were of high quality - assumed they have previously been trimmed by David for assemblies (checked with him and they have not been qc'd after sequencing facility trimmed adapters, QC below)

####QC
ran SolexaQA++ on Ee-STR1 fastq file, trimmed with < phred 20 cutoff
```{bash}
SolexaQA++ dynamictrim /opt/stroma_RNAseq/Ee-PS1_S7_L007_R1_001.fastq.gz --phredcutoff 20 --bwa --directory E.elymi --sanger
```

####index genome with STAR
```{bash}
STAR --runMode genomeGenerate --runThreadN 4 --genomeDir genomes/Eel_728 --genomeLoad LoadAndKeep --genomeFastaFiles genomes/Eel_728/Eel_masked.fasta --outFileNamePrefix Eel_728_
```



###20190322  Friday
####overnight:
3 PS and 3 STR E.elymi files from the first run were successfully trimmed overnight with SolexaQA++
All files were analysed for initial quality with SolexaQA++

####reading
Alpine models GC bias in RNA-seq data and significantly reduces the number of false positives (especially where there are isoforms with high GC content)
Salmon implements same models as Alpine and provides accurate transcript abudance estimates either from its own light mapping or from BAM genomeFastaFiles

####notes
E6255 == NFE728 see David's email :
"Among our collaborators, there are two commonly used strain-naming conventions. U. Kentyck numbers starting with "E" and Noble Foundation ones starting "NfE". The strain used in the study is called E6255 in the U. Kentucky system and NfE728 in the Noble Foundation system."


####finish QC
run 3 INF E.elymi files through trimming with SolexaQA++


####mapping sample reads to E.elymi
STAR --runMode alignReads \
--runThreadN 8 \
--genomeDir genomes/Eel_728/ \      # NB on the forward slash for end of folder
--readFilesIn ../1_QC/E.elymi/........ \
###--readFilesCommand zcat \           # input in .gz files took this out as trimmomatic output uncompressed files

--sjdbGTFfile GFFfilepath \         # GTF/GFF3 file        
--sjdbGTFfeatureExon exon \         # which feature to pull out (3rd column of gff file)
--sjdbGTFtagExonParentTranscript Parent \    #NB this should be transcript_id, except we didn't put them in gff file*  
--outFileNamePrefix mappedfiles/E.elymi_whatever_the_file_prefix_is \
--outSAMtype BAM SortedByCoordinate \

* David made one canonical gene and transcript per gene, this is numbered in the 'Parent' label. No 'transcript_id' or 'gene_id' features numbered

##### should I try this parameter too?
--outSAMstrandField intronMotif \  # filters out reads with non-canonical introns (David used this)
##### example:
```{bash}
STAR --runMode alignReads --runThreadN 8 --genomeDir genomes/Eel_728/ --readFilesIn ../1_QC/E.elymi/Ee-Inf1_S14_L006_R1_001.fastq.trimmed.gz --readFilesCommand zcat --sjdbGTFfile genomes/Eel_728/Epichloe_elymi.gff3 --sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript Parent --outFileNamePrefix mappedfiles/INF1/E.elymi_Ee-Inf1_S14_L006_R1_001_
```



###20190325 Monday
####overnight:
STAR mapped first run E.elymi reads to genomes. Re-running PS1 due to missing gap in command line

####notes
moving forward with Salmon (100s of citations) and DEseq2 (still commonly used)
Had a look at trimmed data with fastqc - adapters still present, need to trim these and re-map

####meeting
meeting with Barry, Daniel and David at 3pm.
Yonathan's project is RNAseq for two mutants which have lost the ability to infect hosts.
Getting sequenced in Novagene. Wants DE analysis done between these and WT.

# Salmon
maps to transcriptome (can make fastq files from bam file)

####Rstudio
set up Rstudio added key and linux libraries to install tidyverse packages within Rstudio (intall.packages() etc.)
see https://cran.r-project.org/bin/linux/ubuntu/README.html
```
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo apt-get install libcurl4-openssl-dev
sudo apt-get install libssl-dev
sudo apt-get install libxml2-dev
```

###20190326 Tuesday
computer sluggish - Murray's assembly took up all the memory!

####fastQC
Ran FastQC on all raw files. (adapters still present)
####trimmomatic
ran trimmomatic on all raw fastq files, removed TrueSeq3 SE adapters, and < phred 20 over window of 5, chucked anything less than 40bp
```
java -jar /home/kate/bin/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 6 -trimlog trimmed/E.elymi/log/Ee-INF1_trim.log /opt/stroma_RNAseq/Ee-Inf1_S14_L006_R1_001.fastq.gz trimmed/E.elymi/Ee-Inf1_S14_L006_R1_001.trim.fastq ILLUMINACLIP:/home/kate/bin/Trimmomatic-0.38/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:40
```
re-ran fastQC on trimmed files

###20190327 Wednesday

####meeting
meeting with Yonathan about his project
he finishes in July, data comes in in 2-3 weeks, needs fast turnaround on Epichloe data to investigate further in the lab
he has FL1 and mutants in a Lolium perenne background
two histone mutants DsetB (K36) and DclrD (K39)
we have the FL1 gene set and a 45K gene set for Lolium. Torben Asp in Denmark has a 70K gene set for Lolium.
Data: 12 samples on 2 lanes of NOVAseq PE 150bp data. 90Gb raw data
Samples are mock (cutting), wt, mutantsX2. Three replicates of each
Samples were cut and innoculated (except mock) and harvested 3 days later. 600 samples in total, pooled into groups of 50 in order to get enough Epichloe in the samples (harvested just the innoculated area of the plant).

Note, the mutants can infect the plants if inocculated with wt, suggesting that the fungus is excreting something that allows it to infect the plant. Plants were innoculated with the combination of mutant/wt died if the mutant grew vigourously.

####mapping adapter-free reads with STAR
changed some parameters: including --quantMode and removed zcat (trimmomatic output was not zipped)
# note transcripts in the gff3 file were identified as Parent items that end in T1, Parent items that don't end in T1 are genes
# one canonical transcript per gene was annotated in this gff3 file
# here identified transcripts in gff3 by looking at 'Parent' item
# --readFilesCommand zcat  read in zipped files
# --sjdbGTFfile annotation GFF3 file
# --sjdbGTFfeatureExon exon feature to read
# --sjdbGTFtagExonParentTranscript Parent    the identifier for transcripts in th GFF3 file
# --quantMode TranscriptomeSAM output SAM/BAM alignments to transcriptome into a separate file

```
STAR --runMode alignReads --runThreadN 4 --genomeDir genomes/Eel_728/ --readFilesIn ../1_QC/trimmed/E.elymi/Ee-Inf1_S14_L006_R1_001.trim.fastq --sjdbGTFfile genomes/Eel_728/Epichloe_elymi.gff3 --sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript Parent --outFileNamePrefix mappedfiles/INF1/E.elymi_Ee-Inf1_S14_L006_R1_001_  --quantMode TranscriptomeSAM

```

####bamtofastq of SolexaQCd reads to play with in salmon
made a trial set of fastq reads to use Salmon in mapping mode

####salmon
ran initial test runs of Salmon - it is really fast
was running in mapping mode, lots of reads (fragments) got thrown out, most likely due to presence of adapters especially for inflorescene Samples
see output overview in 3_Salmon/trial_salmon/salmon_SolexaQCd_fastq_mapping.ods

NOTE: read lior pachters take-down of Salmon/Sailfish etc. and the rebuttal.
take home message - kallisto, sailfish and salmon have very similar output, gc content bias is addressed in Salmon

####updated R in ubuntu
change the /etc/apt/source.list (note needed to use bionic to match hosts in list)
```
sudo emacs /etc/apt/sources.list
```
add line:
deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/
update
```
sudo apt-get update
sudo apt-get install r-base
sudo apt-get install r-base-dev
```

#### installed DEseq2
see 3_salmon/trial_salmon/E.elymi_DEseq2_trial

####tomorrow
start with looking back at fastQC output before and after


###20190328 Thursday

####overnight
STAR mapped adapter trimmed reads for INF1 INF2 PS1, continuing with PS2 etc... on 8 processors

####trial_salmon
lots of reads got thrown out for SolexaQA trimmed data. Try trimmomatic trimmed + STAR mapped with -quant mode and also with mapped mode in salmon (get fastq from bam)
(note there was an error in the bam to fastq for PS1 only [W::sam_read1] Parse error at line 10377917)
-quant mode returns fewer mapped reads (except for PS2???????)

Tried running trimmomatic fastq output directly with salmon
-quant mode of Trimmomatic output gives similar number to fastq reads from trimmomatic+STAR bam files
-i.e. broadly speaking a large overlap between STAR and Salmon mapping, but <10% unique to each mapper

Ran basic quant options for salmon for all trimmomatic trimmed files
```
salmon quant -p 8 -i transcriptome/Epichloe_elymi_trans_index -l U -r ../1_QC/trimmed/E.elymi/Ee-Inf2_S15_L006_R1_001.trim.fastq --validateMappings --gcBias -o trial_trim_files/Epichloe_elymi_INF2_quant
```

###20190329 Friday

made setup scripts for all files fastq, trimmomatic and salmon quant steps
lab meeting 11-12.30
postgraduate drinks and boardgames @ 4



###20190401 Monday
re did config files for automated setup
back up all raw data files
started trimmomatic and fastq of everything
completed MPI assessment requirements

DESeq2 - identify clusters



###DONE
analyse quality of reads with SolexaQA++
trim reads < phred 20
map to Epichloe genomes with STAR
re-do QC with Fastqc and trimmomatic (adapters removed, trimmed phred < 20 over window size 5)

###TO DO
quantify abundance with Salmon
DE with DEseq2

###ARCHIVED ideas
see if there is a grass genome to map to to remove plant material?  Ryegrass (lolium) have 45K gene set (70K set in Denmark)
remove bacterial and other genes?

###QUESTIONS
what is masked in the Epichloe genome assembly 'masked' fasta files?
what are the unidentified samples in the illumina run?
what is the estimated fragment length in salmon?
How does STAR differ from bwa?
###Other Questions:
How do you assess genome re-shuffling?

####POSSIBLE FUTURE PROJECTs:
Indonesian data set
  - look at admixtures
Epichloe genomics (whole genomes available, hiC and expression data)
  - see how whole genomes re-shuffle
    - suspect it may happen at repetitive elements
  - look at gene expression in parents and hybrids. Any deviation from the expected value may be under selection?
    - personally I think there are too many variables in this, random noise from disruption of regulatory systems etc. would feature heavily. You would need a lot of replicates to support any possible selection.
    - need to know more about how hybridisation works.
