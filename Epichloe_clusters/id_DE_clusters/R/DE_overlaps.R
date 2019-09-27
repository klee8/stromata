#AIM:
#  - take total number of expressed genes and total number of significantly differentially expressed genes for each species
#  - randomise the DE genes multiple times (e.g. 1 million)
#  - check to see if there are overlaps of DE genes between species

#Adapted from https://mac-theobio.github.io/QMEE/permutation_examples.html


#!/usr/bin/env Rscript
#install.packages("tidyverse")
#install.packages("ggpubr")
#install.packages("data.table")
#install.packages("doParallel")
library(doParallel)
library(tidyverse)
library(ggpubr)
library(data.table)

##### for testing only
#path_DE_fun <- "/media/kate/Massey_linux_onl/projects/stromata_analysis/Epichloe_clusters/id_DE_clusters/DE_fun.R"
#datadir <- "/media/kate/Massey_linux_onl/projects/results/stromata/Epichloe_clusters/rfmt_core_gene_sets/"
#resdir <- "/media/kate/Massey_linux_onl/projects/results/stromata/Epichloe_clusters/id_DE_clusters/"
#args <- c("STR_PS", 10, path_DE_fun, datadir, resdir)

## read in command line arguments
  args = commandArgs(trailingOnly=TRUE)

if (length(args)<1) {
  stop("At least one argument must be supplied: <base_filename>, <number of simulations [default = 10,000]> <path/DE_fun.R> <datadir> <resultsdir>.n", call.=FALSE)
} else {
  basefilename <- args[1]
  nsim <- as.numeric(args[2])
  DE_fun <- args[3]
  datadir <- args[4]
  resdir <- args[5]
  infile <- paste(datadir, "core_genes_", args[1], "_rfmt.txt", sep = "")
  raw <- paste("DE_perm_results/", args[1], "_raw_counts.txt", sep = "")
  outfile <- paste("DE_perm_results/", args[1], "_permuted_DE_results.txt", sep = "")
  summaryfile <- paste("DE_perm_results/", args[1], "_permuted_DE_summary.txt", sep = "")
  graphout <- paste("DE_perm_results/", args[1], "_permuted_DE_results.pdf", sep = "")

}

source(DE_fun)
dir.create("DE_perm_results")
print(getwd())
print(DE_fun)
print(paste("infile: ", infile))
print(paste("outfile: ", outfile))
print(paste("graph file: ", graphout))
print(paste("nsim: ", nsim))
print(paste(resdir, basefilename, sep = "" ))
setwd(paste(resdir, basefilename, sep = "" ))

#####      Read in the data from the core gene set analysis (this shows a line for every ortholog described in the gene set)

data_rfmt <- read.delim(infile, header = TRUE, sep = "\t")

# remove lines with no fold change reading
data <- data_rfmt[!is.na(data_rfmt$log2fc),]

# remove any duplicate rows
data <- data %>% distinct

# add column which sets DE at 1 for log2FC >= 1 and svalue <= 0.005, set DE at 2 for log2FC >= 2 and svalue <= 0.005
data <- DE_categories(data)

# flag direction of fold change
data <- data %>% mutate(dir_FC = ifelse((log2fc < 0), "-", "+"))

# subset data
data <- data[, c("species", "orthogroup", "gene_id", "dir_FC", "DE")]

# total genes and total number of DE genes per species
raw_counts <- species_raw_DE_counts(data)

# find ortholog overlaps
max_spp <- (length(unique(data$species)))
overlaps <- DE_overlaps(data)

# OBSERVED DATA RESULTS
# count up total number of each kind of overlap
obs <- total_overlaps(overlaps, max_spp)
names(obs) <-  c("core_up", "FC2_up", "FC4_up", "core_down", "FC2_down", "FC4_down")
