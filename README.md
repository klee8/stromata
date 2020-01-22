## Epichloe stromata RNA-Seq project analysis 
### Kate Lee, 2019 

### AIMS:
- Find differentially expressed transcripts in three different treatments of three different Epichloe-host symbiotic pairs
- Find core set of genes that are always up/down regulated in stroma
- Find changes in the inflorescence vs psudostem (and see which of these match stroma changes)
- Identify strain-specific clusters of secondary metabolites (especially those that are down-regulated in stroma)
- Compare psudostem expression with previous studies


### DATA:
- Epichloe samples
    - E.typhina_E8: Epichloe typhina strain E8 (obligate stromata former), from Rikers. Causes choke on ryegrass. (454 sequencing for this strain (from Chris Shardel 2013 comparative geneome paper.)
    - E.festucae_E2368: Epichloe festucae strain E2368, created by Chris Shardel by crossing E189 X E434 and backcrossing the F1 generation to E189. Fl1 + lolium is the lab's main Epichloe-grass model.
    - E.elymi_NfE728: Epichloe elymi strain NFe728, from Carolyn Young's lab (intermediate stroma former)

- David Winter's assemblies for each Epichloe strain
    - masked fasta file (genome)
    - gff3 files
    - annotations text files
    - protein fasta file
    - trascript fasta files (transcriptome)
    - (Epichloe illumina runs listed in stroma_files.csv)


### OVERVIEW:
#### Epichloe DE
- 0_raw_data (on SRA database, Bioproject PRJNA554133)
- 1_QC
    - ran fastQC on all files
    - trimmomatic used to remove illumina adapters (TruSeq2 SE) and trim < phred 20 (SLIDINGWINDOW:5:20)
    - re-ran fastQC on trimmed files to check them
- 2_Salmon
    - Quant-map against Epichloe gene sets from David Winter.
- 3_DESeq2  
    - take in quant counts and find differentially expressed genes with DEseq2
    - output a list of genes expressed log2foldchange > 1 and svalue < 0.0005, with their positions (from gff3)
    - add annotations from pannzer, SignalP, EffectorP 
    - add annotations from David Winters funnanotate pipeline for Epichloe.
    - identify core genes found in more than one species

#### Epichloe clusters 
Use DE output to identify clusters of genes that are up or down regulated. Illustrate with graphs and tabulate results.
    
    
#### Ryegrass DE
Same as above for Epichloe DE but uses the Lolium perenne gene set from Dupont et al., (2015).

#### Plant vs culture DE
Same DE pipeline as above, but using the plant and culture data from Dupont et al., (2015).

