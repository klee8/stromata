#AIM:
#  - take total number of expressed genes and total number of significantly differentially expressed genes for each species
#  - randomise the DE genes multiple times (e.g. 1 million)
#  - check to see if there are overlaps of DE genes between species

#Adapted from https://mac-theobio.github.io/QMEE/permutation_examples.html


#!/usr/bin/env Rscript
#install.packages("tidyverse")
#install.packages("ggpubr")
#install.packages("data.table")
library(tidyverse)
library(ggpubr)
library(data.table)

## for testing only
#setwd("/media/kate/Massey_linux_onl/projects/results/stromata/Epichloe_clusters/id_DE_clusters/STR_PS")
#args <- c("STR_PS", 100)

## read in command line arguments
args = commandArgs(trailingOnly=TRUE)
dir.create("DE_perm_results")
if (length(args)<1) {
  stop("At least one argument must be supplied: <base_filename>, <number of simulations [default = 10,000]>.n", call.=FALSE)
} else {
  infile <- paste("../../rfmt_core_gene_sets/core_genes_", args[1], "_rfmt.txt", sep = "")
  raw <- paste("DE_perm_results/", args[1], "_raw_counts.txt", sep = "")
  outfile <- paste("DE_perm_results/", args[1], "_permuted_DE_results.txt", sep = "")
  summaryfile <- paste("DE_perm_results/", args[1], "_permuted_DE_summary.txt", sep = "")
  graphout <- paste("DE_perm_results/", args[1], "_permuted_DE_results.pdf", sep = "")
  nsim <- ifelse( length(args) ==2, as.numeric(args[2]), 10000)
}
print(getwd())
print(paste("infile: ", infile))
print(paste("outfile: ", outfile))
print(paste("graph file: ", graphout))
print(paste("nsim: ", nsim))

#####      Read in the data from the core gene set analysis (this shows a line for every ortholog described in the gene set)

data_rfmt <- read.delim(infile, header = TRUE, sep = "\t")

# remove lines with no fold change reading
data <- data_rfmt[!is.na(data_rfmt$log2fc),]

# remove any duplicate rows
data <- data %>% distinct

# set DE at 1 for log2FC >= 1 and svalue <= 0.005, set DE at 2 for log2FC >= 2 and svalue <= 0.005
data <- data %>% mutate(DE = ifelse((log2fc >= 2) &(svalue_1 <= 0.005), 2,
                                        ifelse((log2fc >= 1) &(svalue_1 <= 0.005), 1, 0)),
                            DE2 = ifelse((log2fc >= 1) &(svalue_1 <= 0.005), 1, 0),
                            DE4 = ifelse((log2fc >= 2) &(svalue_1 <= 0.005), 1, 0))

# flag direction of fold change
data <- data %>% mutate(dir_FC = ifelse((log2fc < 0), "-", "+"))

# subset data
data <- data[, c("species", "orthogroup", "gene_id", "dir_FC", "DE", "DE2", "DE4")]

# total genes and total number of DE genes per species
raw_counts <- data %>% group_by(species) %>% summarise(total = n(), total_DE2 = sum(DE2), total_DE4 = sum(DE4))

write.table(raw_counts, raw, quote = FALSE, row.names = FALSE)

# find rtholog overlaps 
DE_overlaps <- function(data) {
  res <- data %>%
    group_by(orthogroup) %>%
    mutate(number_spp = length(unique(species)),
           FC2_overlaps = ifelse((length(unique(species)) > 1) & (sum(DE2) == length(unique(species))), 1, 0),
           FC4_overlaps = ifelse((length(unique(species)) > 1) & (sum(DE4) == length(unique(species))), 1, 0),
           core_overlaps = ifelse((length(unique(species)) > 1) & (sum(DE2) == (length(unique(species)))) & (sum(DE4) >= 1), 1, 0))
  return(res)
}

data <- DE_overlaps(data)

# count up total number of each kind of overlap
total_overlaps <- function(res){
  FC4_up <- nrow(filter(res, (FC4_overlaps == 1) & (dir_FC == "+")))
  FC4_down <- nrow(filter(res, (FC4_overlaps == 1) & (dir_FC == "-")))
  FC2_up <- nrow(filter(res, (FC2_overlaps == 1) & (dir_FC == "+")))
  FC2_down <- nrow(filter(res, (FC2_overlaps == 1) & (dir_FC == "-")))
  core_up <- nrow(filter(res, (core_overlaps == 1) & (dir_FC == "+")))
  core_down <- nrow(filter(res, (core_overlaps == 1) & (dir_FC == "-")))
  tmp <- c(FC4_up, FC4_down, FC2_up, FC2_down, core_up, core_down)
  return(tmp)
}
# observed data results
obs <- total_overlaps(data)
names(obs) <-  c("FC2_up", "FC2_down", "FC4_up", "FC4_down", "core_up", "core_down")
obs

############  PERMUTE DATA

set.seed(101) ## for reproducibility

# set up permutation results outfile
headers <-  c("FC2_up", "FC2_down", "FC4_up", "FC4_down", "core_up", "core_down")
fwrite(as.list(headers), file = outfile, sep = "\t")

# permutation test
DE_permutation <- function(data) {
  ## standard approach: scramble response value
  bdat <- data.frame()
  for (spp in unique(data$species)) {
    perm <- sample(nrow(data[data$species == spp,]))
    tmp <- transform(data[data$species == spp,], DE = DE[perm])
    tmp <- transform(tmp, DE2 = ifelse((DE == 1) | (DE == 2), 1, 0), DE4 = ifelse((DE == 2), 1, 0))
    bdat <- rbind(bdat, tmp)
  }
  ## check the number of ortholog overlaps and store them
  res <- DE_overlaps(bdat)
  tmp <- as.vector(total_overlaps(res))
  fwrite(as.list(tmp), file = outfile, sep = "\t", append = TRUE)
  return(tmp)
}

# run permutation function with foreach parellisation 
n.cores <- detectCores()
registerDoParallel(n.cores)
permutations <- foreach(k = 1:nsim, .combine = "rbind") %dopar% DE_permutation(data)
permutations <- as.data.frame(permutations)
colnames(permutations) <-  c("FC2_up", "FC2_down", "FC4_up", "FC4_down", "core_up", "core_down")
rownames(permutations) <- NULL
head(permutations)



pvalue_FC2_up <- ((nrow(permutations[permutations$FC2_up >= obs["FC2_up"], ])+1)/(nsim + 1))
pvalue_FC2_down <- ((nrow(permutations[permutations$FC2_down <= obs["FC2_down"], ])+1)/(nsim + 1))
pvalue_FC4_up <- ((nrow(permutations[permutations$FC4_up >= obs["FC4_up"], ])+1)/(nsim + 1))
pvalue_FC4_down <- ((nrow(permutations[permutations$FC4_down <= obs["FC4_down"], ])+1)/(nsim + 1))
pvalue_core_up <- ((nrow(permutations[permutations$core_up >= obs["core_up"], ])+1)/(nsim + 1))
pvalue_core_down <- ((nrow(permutations[permutations$core_down <= obs["core_down"], ])+1)/(nsim + 1))

summary <- data.frame()
for (i in colnames(permutations)){
    overlap <- as.factor(i)
    range <- as.integer(range(permutations[, c(i)]))
    pvalue <- as.numeric(get(paste("pvalue_", i, sep = "")))
    row <- c(overlap, range[1], range[2], pvalue)
    summary <- rbind(summary, row)
}
colnames(summary) <- c("overlap", "range_start","range_end", "pvalue")
summary$overlap <- colnames(permutations)
write.table(summary, summaryfile, quote = FALSE, row.names = FALSE)

results <- permutations
#####  GRAPH PERMUTATIONS

up_2FC <- ggplot(results) +
            geom_bar(aes(x = FC2_up)) +
            geom_vline(xintercept = obs["FC2_up"], color = "red", size=1)

down_2FC <- ggplot(results) +
              geom_bar(aes(x = FC2_down)) +
              geom_vline(xintercept = obs["FC2_down"], color = "red", size=1)

up_4FC <- ggplot(results) +
            geom_bar(aes(x = FC4_up)) +
            geom_vline(xintercept = obs["FC4_up"], color = "red", size=1)

down_4FC <- ggplot(results) +
              geom_bar(aes(x = FC4_down)) +
              geom_vline(xintercept = obs["FC4_down"], color = "red", size=1)

up_core <- ggplot(results) +
              geom_bar(aes(x = core_up)) +
              geom_vline(xintercept = obs["core_up"], color = "red", size=1)

down_core <- ggplot(results) +
              geom_bar(aes(x = core_down)) +
              geom_vline(xintercept = obs["core_down"], color = "red", size=1)

# see combined fig example http://www.sthda.com/english/articles/32-r-graphics-essentials/126-combine-multiple-ggplots-in-one-graph/
figure <- ggarrange(up_2FC, down_2FC, up_4FC, down_4FC, up_core, down_core,
                    labels = c("A", "B", "C", "D", "E", "F"),
                    ncol = 2, nrow = 3)
ggsave(graphout, plot = figure)


