## title: "reformat_core_gene_set.Rmd"
## author: "Kate Lee"
## date: "03/07/2019"

#!/usr/bin/env Rscript



library(tidyverse)

### set default arguments
datadir <- "/media/kate/Massey_linux_onl/projects/results/stromata/Epichloe_stromata_DE/3_DEseq2/"
resdir <- "/media/kate/Massey_linux_onl/projects/results/stromata/Epichloe_clusters/rfmt_core_gene_sets/"

### Read in command line arguments

args = commandArgs(trailingOnly=TRUE)
if (length(args) == 2) {
  datadir = args[1]
  resdir = args[2]
}

  
str_ps <- read.delim(paste(datadir, "STR_PS/core_gene_set_STR_PS_ann.txt", sep = ""), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
str_ps_ele <- str_ps %>% select(elymi_Contig, elymi_Start, elymi_Stop, elymi_Strand, elymi_gene_id, orthogroup,  elymi_apeglm_log2FC, elymi_apeglm_lfcSE, elymi_apeglm_1_svalue, elymi_apeglm_2_svalue) %>% mutate(species = "elymi")
  str_ps_ele <- str_ps_ele[!is.na(str_ps_ele$elymi_gene_id),]
colnames(str_ps_ele) <- c("contig", "start", "stop", "strand", "gene_id","orthogroup", "log2fc", "lfcSE", "svalue_1", "svalue_2", "species")
str_ps_fes <- str_ps %>% select(festucae_Contig, festucae_Start, festucae_Stop, festucae_Strand, festucae_gene_id, orthogroup, festucae_apeglm_log2FC, festucae_apeglm_lfcSE, festucae_apeglm_1_svalue, festucae_apeglm_2_svalue) %>% mutate(species = "festucae")
str_ps_fes <- str_ps_fes[!is.na(str_ps_fes$festucae_gene_id),]
colnames(str_ps_fes) <- c("contig", "start", "stop", "strand", "gene_id","orthogroup", "log2fc", "lfcSE", "svalue_1", "svalue_2", "species")
str_ps_typ <- str_ps %>% select(typhina_Contig, typhina_Start, typhina_Stop, typhina_Strand, typhina_gene_id, orthogroup, typhina_apeglm_log2FC, typhina_apeglm_lfcSE, typhina_apeglm_1_svalue, typhina_apeglm_2_svalue) %>% mutate(species = "typhina")
str_ps_typ <- str_ps_typ[!is.na(str_ps_typ$typhina_gene_id),]
    colnames(str_ps_typ) <- c("contig", "start", "stop", "strand", "gene_id","orthogroup", "log2fc", "lfcSE", "svalue_1", "svalue_2", "species")
str_ps_rfmt <- rbind(str_ps_ele, str_ps_fes, str_ps_typ)

write.table(str_ps_rfmt, paste(resdir, "core_genes_STR_PS_rfmt.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)

str_inf <- read.delim(paste(datadir, "STR_INF/core_gene_set_STR_INF_ann.txt", sep = ""), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
str_inf_ele <- str_inf %>% select(elymi_Contig, elymi_Start, elymi_Stop, elymi_Strand, elymi_gene_id, orthogroup,  elymi_apeglm_log2FC, elymi_apeglm_lfcSE, elymi_apeglm_1_svalue, elymi_apeglm_2_svalue)  %>% mutate(species = "elymi")
str_inf_ele <- str_inf_ele[!is.na(str_inf_ele$elymi_gene_id),]
colnames(str_inf_ele) <- c("contig", "start", "stop", "strand", "gene_id","orthogroup", "log2fc", "lfcSE", "svalue_1", "svalue_2", "species")
str_inf_fes <- str_inf %>% select(festucae_Contig, festucae_Start, festucae_Stop, festucae_Strand, festucae_gene_id, orthogroup, festucae_apeglm_log2FC, festucae_apeglm_lfcSE, festucae_apeglm_1_svalue, festucae_apeglm_2_svalue) %>% mutate(species = "festucae")
str_inf_fes <- str_inf_fes[!is.na(str_inf_fes$festucae_gene_id),]
colnames(str_inf_fes) <- c("contig", "start", "stop", "strand", "gene_id","orthogroup", "log2fc", "lfcSE", "svalue_1", "svalue_2", "species")
str_inf_rfmt <- rbind(str_inf_ele, str_inf_fes)

write.table(str_inf_rfmt, paste(resdir, "core_genes_STR_INF_rfmt.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)

inf_ps <- read.delim(paste(datadir, "INF_PS/core_gene_set_INF_PS_ann.txt", sep = ""), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
inf_ps_ele <- inf_ps %>% select(elymi_Contig, elymi_Start, elymi_Stop, elymi_Strand, elymi_gene_id, orthogroup,  elymi_apeglm_log2FC, elymi_apeglm_lfcSE, elymi_apeglm_1_svalue, elymi_apeglm_2_svalue) %>% mutate(species = "elymi")
inf_ps_ele <- inf_ps_ele[!is.na(inf_ps_ele$elymi_gene_id),]
colnames(inf_ps_ele) <- c("contig", "start", "stop", "strand", "gene_id","orthogroup", "log2fc", "lfcSE", "svalue_1", "svalue_2", "species")
inf_ps_fes <- inf_ps %>% select(festucae_Contig, festucae_Start, festucae_Stop, festucae_Strand, festucae_gene_id, orthogroup, festucae_apeglm_log2FC, festucae_apeglm_lfcSE, festucae_apeglm_1_svalue, festucae_apeglm_2_svalue) %>% mutate(species = "festucae")
inf_ps_fes <- inf_ps_fes[!is.na(inf_ps_fes$festucae_gene_id),]
colnames(inf_ps_fes) <- c("contig", "start", "stop", "strand", "gene_id","orthogroup", "log2fc", "lfcSE", "svalue_1", "svalue_2", "species")
inf_ps_rfmt <- rbind(inf_ps_ele, inf_ps_fes)

write.table(inf_ps_rfmt, paste(resdir, "core_genes_INF_PS_rfmt.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
