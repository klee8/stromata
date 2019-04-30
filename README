# Epichloe stroma RNA-Seq project analysis, Kate Lee, 2019 #

AIM:
  Find differentially expressed transcripts in three different treatments of three different Epichloe-host symbiotic pairs
  Find core set of genes that are always up/down regulated in stroma
  Find changes in the inflorescence vs psudostem (and see which of these match stroma changes)
  Identify strain-specific clusters of secondary metabolites (especially those that are down-regulated in stroma)
  Compare psudostem expression with previous studies


DATA:
  Epichloe samples
    E.typhina_E8: Epichloe typhina strain E8 (obligate stromata former), from Rikers. Causes choke on ryegrass.
      454 sequencing available for this strain (from Chris) see 2013 comparative geneome paper.
    E.festucae_E2368: Epichloe festucae strain E2368, created by Chris Shardel by crossing E189 X E434 and backcrossing
      the F1 generation to E189. Fl1 + lolium is the lab's main Epichloe-grass model.
    E.elymi_NfE728: Epichloe elymi strain NFe728, from Carolyn Young's lab (intermediate stroma former)

  David Winter's assemblies for each Epichloe strain
    masked fasta file (genome)
    gff3 files
    annotations text files
    protein fasta file
    trascript fasta files (transcriptome)

  Epichloe illumina runs listed in stroma_files.csv


PIPELINE:
  0_raw_data (mostly held in /opt or /media due to space issues)
  1_QC      ran fastQC on all files
            trimmomatic used to remove illumina adapters (TruSeq2 SE) and trim < phred 20 (SLIDINGWINDOW:5:20)
            re-ran fastQC on trimmed files to check them

  2_Salmon  random order bam file needed for use in Salmon (can randomise with samtools collate (orders by read name i.e. not in genome position order))
            mapped bam files are not in positional order, use as they are.
            Quant-map

  3_DESeq2  Take in quant counts and find differentially expressed genes with DEseq2
            Output a list of genes expressed log2foldchange > 2 and pvalue < 0.01, with their positions (from gff3)

  4_clusters  use occultercut to identify AT rich (i.e. most likely repeat regions) to look for clusters that have been split up by repeat regions
            bwa alignment of genomes against more complete genome (David has pacbio data) to see if some segments are near enough to have genes from the same cluster

  5_GOI     get fasta file of genes of interest
            get protein sequence of genes of interest (from David's annotations.txt)
            check interproscan for domains of genes of interest - can put these through GO terms




ARCHIVED
QC        initially tried SolexQC, but found this missed illumina adapters which were still on the threads
STAR      initially ran STAR alignment against genome (found some new junctions). New splicing is not of interest so abandoned this.
          Ran STAR against transcriptome (Salmon can also do this, but it seems better to do this with a
          full read mapping algorithm and extract Epichloe reads first to minimise possibility of getting contaminant reads)
          Salmon needs reads mapped to transcriptome. convert bam to fastq and try quant mapping with fastq output
          Found there was not much difference (< 10%) between mapping with STAR and quant-mapping with salmon, and just mapping quant mapping with salmon
