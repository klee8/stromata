# title: "formating_for_SMURF"
# author: "Kate"
# date: "5 June 2019"

#install.packages("seqinr")
library("seqinr")
#install.packages("tidyverse")
library("tidyverse")


### Set defaults
datadir <- "/media/kate/Massey_linux_onl/projects/data/gene_model_sets/"
resdir <- "/media/kate/Massey_linux_onl/projects/STROMATA/results/Epichloe_clusters/SMURF/"

### Read in command line arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) > 0) {
  datadir = args[1]
  if (length(args) >1)
  resdir = args[2]
}


###  read in fasta file and get headers
elymi_fa <- read.fasta(paste(datadir, "Epichloe_elymi_NfE728/Epichloe_elymi_NfE728.proteins.fa", sep = ""))
elymi_fa_df <- data.frame(Fragments=names(elymi_fa), Seqs=unlist(getSequence(elymi_fa, as.string=T)))
elymi_fa_df$gene_id <- gsub("-T1", "", elymi_fa_df$Fragments)

festucae_fa <- read.fasta(paste(datadir, "Epichloe_festucae_E2368/Epichloe_festucae_E2368.proteins.fa", sep = ""))
festucae_fa_df <- data.frame(Fragments=names(festucae_fa ), Seqs=unlist(getSequence(festucae_fa , as.string=T)))
festucae_fa_df$gene_id <- gsub("-T1", "", festucae_fa_df$Fragments)

typhina_fa <- read.fasta(paste(datadir, "Epichloe_typhina_E8/Epichloe_typhina_E8.proteins.fa", sep = ""))
typhina_fa_df <- data.frame(Fragments=names(typhina_fa ), Seqs=unlist(getSequence(typhina_fa , as.string=T)))
typhina_fa_df$gene_id <- gsub("-T1", "", typhina_fa_df$Fragments)

###   read in gff3 file and get mRNA start and end positions (5'-3')
elymi_gff3 <- read.delim("data/Epichloe_elymi_NfE728.gff3", header = FALSE, sep = "\t", comment.char = '#')
elymi_gff3 <- elymi_gff3[elymi_gff3$V3 == "mRNA",]
elymi_gff3$trans_name <- gsub(";.*", "", elymi_gff3$V9 )
elymi_gff3$trans_name <- gsub("ID=", "", elymi_gff3$trans_name )
elymi_gff3$product <- gsub(".*product=", "", elymi_gff3$V9 )
elymi_gff3$start <- ifelse(elymi_gff3$V7=="+", elymi_gff3$V4, elymi_gff3$V5)
elymi_gff3$end <- ifelse(elymi_gff3$V7=="+", elymi_gff3$V5, elymi_gff3$V4)
elymi_gff3 <- elymi_gff3[ , c("trans_name", "V1", "start", "end", "product")]
colnames(elymi_gff3) <-  c("trans_name", "contig", "start", "end", "product")

write.table(elymi_gff3, paste(resdir, "data/E.elymi_ann_for_SMURF.txt", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t")

festucae_gff3 <- read.delim("data/Epichloe_festucae_E2368.gff3", header = FALSE, sep = "\t", comment.char = '#')
festucae_gff3 <- festucae_gff3[festucae_gff3$V3 == "mRNA",]
festucae_gff3$trans_name <- gsub(";.*", "", festucae_gff3$V9 )
festucae_gff3$trans_name <- gsub("ID=", "", festucae_gff3$trans_name )
festucae_gff3$product <- gsub(".*product=", "", festucae_gff3$V9 )
festucae_gff3$start <- ifelse(festucae_gff3$V7=="+", festucae_gff3$V4, festucae_gff3$V5)
festucae_gff3$end <- ifelse(festucae_gff3$V7=="+", festucae_gff3$V5, festucae_gff3$V4)
festucae_gff3 <- festucae_gff3[ , c("trans_name", "V1", "start", "end", "product")]
colnames(festucae_gff3) <-  c("trans_name", "contig", "start", "end", "product")

write.table(festucae_gff3, paste(resdir, "data/E.festucae_ann_for_SMURF.txt", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t")

typhina_gff3 <- read.delim("data/Epichloe_typhina_E8.gff3", header = FALSE, sep = "\t", comment.char = '#')
typhina_gff3 <- typhina_gff3[typhina_gff3$V3 == "mRNA",]
typhina_gff3$trans_name <- gsub(";.*", "", typhina_gff3$V9 )
typhina_gff3$trans_name <- gsub("ID=", "", typhina_gff3$trans_name )
typhina_gff3$product <- gsub(".*product=", "", typhina_gff3$V9 )
typhina_gff3$start <- ifelse(typhina_gff3$V7=="+", typhina_gff3$V4, typhina_gff3$V5)
typhina_gff3$end <- ifelse(typhina_gff3$V7=="+", typhina_gff3$V5, typhina_gff3$V4)
typhina_gff3 <- typhina_gff3[ , c("trans_name", "V1", "start", "end", "product")]
colnames(typhina_gff3) <-  c("trans_name", "contig", "start", "end", "product")

write.table(typhina_gff3, paste(resdir, "data/E.typhina_ann_for_SMURF.txt", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t")

#str(elymi_fa_df)
#str(elymi_gff3)
#str(festucae_fa_df)
#str(festucae_gff3)
#str(typhina_fa_df)
#str(typhina_gff3)
