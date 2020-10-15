
# Comparative genomics analysis

### AIM: look for presence of core genes from stromata DE analysis in a range of sexual and asexual Epichloe genomes

### DATA: data from Schardl, Young, Cox and Scott labs (sexual, unknown and asexual species)

S EelyNfE728
S EfesE2368
S EtypE8
S EamaE57
S EbacE1031
S EbraE4804
S EbroE502
S EbroE7561
S EelyE56
S EelyE757
S EglyE277
S EmolE3601
S EsylE7368
S EtypE5819
S EtypORE04
U EamaNFe708
U EamaE4668
U EbroAL0434
U Eely732
U EfesE894
A EaotE899
A EbroNFe1
A EfesAR37
A EfesAR48
A EfesAR5
A EfesAR1_2
A EfesAR6
A EfesAR3002
A EfesAR3018
A EfesAR3060
A EganE7080
A EineE818
A EsylE832
A EtypE4646
A EtypNFe76


#### METHOD
+ Used exonerate exonerate version 2.4.0 to identify genes in each genome that matched the Epichloe elymi NFe728 with relaxed parameters (--percent 50). 
+ Also did a BLAST 2.6.0+ search for the same protein in each genome (evalue 5).
+ Used positional information from exonerate gff3 file to get sequences from the genome fasta files. Used blast output positions where exonerate was not available (Added all exonerate hits for each elymi core protein to a fasta file and where there was none, added the top blast hit.)
+ Aligned genome sequences using MUSCLE v3.8.31
+ Added sex/unknown/asex column and species column to headers of the aligned fasta file
+ Took Epichloe elymi aligned sequence and gff3 information to create a sudo sequence with exon positions


#### ARCHIVED
Tried to annotate everything in funannotate, however not all genomes were suitable for annotation



