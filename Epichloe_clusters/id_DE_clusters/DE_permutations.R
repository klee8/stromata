#AIM:
#  - take total number of expressed genes and total number of significantly differentially expressed genes for each species
#  - randomise the DE genes multiple times (e.g. 1 million)
#  - check to see if there are overlaps of DE genes between species

#Adapted from https://mac-theobio.github.io/QMEE/permutation_examples.html


#!/usr/bin/env Rscript
#install.packages("tidyverse")
#install.packages("ggpubr")
library(tidyverse)
library(ggpubr)


## for testing only
#setwd("/media/kate/Massey_linux_onl/projects/results/stromata/Epichloe_clusters/id_DE_clusters/STR_PS")
#args <- c("STR_PS", 10)

## read in command line arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
  stop("At least one argument must be supplied: <base_filename>, <number of simulations [default = 10,000]>.n", call.=FALSE)
} else {
  infile <- paste("../../rfmt_core_gene_sets/core_genes_", args[1], "_rfmt.txt", sep = "")
  raw <- paste(args[1], "_raw_counts.txt", sep = "")
  outfile <- paste(args[1], "_permuted_DE_results.txt", sep = "")
  summaryfile <- paste(args[1], "_permuted_DE_summary.txt", sep = "")
  graphout <- paste(args[1], "_permuted_DE_results.pdf", sep = "")
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

# total number of ortholog overlaps in observed data
data <- data %>%
  group_by(orthogroup) %>%
  mutate(number_spp = length(unique(species)),
         FC2_overlaps = ifelse((length(unique(species)) > 1) & (sum(DE2) == length(unique(species))), 1, 0),
         FC4_overlaps = ifelse((length(unique(species)) > 1) & (sum(DE4) == length(unique(species))), 1, 0),
         core_overlaps = ifelse((length(unique(species)) > 1) & (sum(DE2) == (length(unique(species)))) & (sum(DE4) >= 1), 1, 0))


############  PERMUTE DATA

set.seed(101) ## for reproducibility
## set aside space for permutation results
FC4_up <- numeric(nsim)
FC4_down <- numeric(nsim)
FC2_up <- numeric(nsim)
FC2_down <- numeric(nsim)
core_up <- numeric(nsim)
core_down <- numeric(nsim)
permutations <- data.frame()

for (i in 1:nsim) {
  ## standard approach: scramble response value
  bdat <- data.frame()
  for (spp in unique(data_rfmt$species)) {
    perm <- sample(nrow(data[data$species == spp,]))
    tmp <- transform(data[data$species == spp,], DE = DE[perm])
    tmp <- transform(tmp, DE2 = ifelse((DE == 1) | (DE == 2), 1, 0), DE4 = ifelse((DE == 2), 1, 0))
    bdat <- rbind(bdat, tmp)
  }

  ## check the number of ortholog overlaps and store them
  res <- bdat %>% group_by(orthogroup) %>%
    mutate(number_spp = length(unique(species)),
           FC2_overlaps = ifelse((length(unique(species)) > 1) & (sum(DE2) == length(unique(species))), 1, 0),
           FC4_overlaps = ifelse((length(unique(species)) > 1) & (sum(DE4) == length(unique(species))), 1, 0),
           core_overlaps = ifelse((length(unique(species)) > 1) & (sum(DE2) == (length(unique(species)))) & (sum(DE4) >= 1), 1, 0))
  FC4_up[i] <- nrow(filter(res, (FC4_overlaps == 1) & (dir_FC == "+")))
  FC4_down[i] <- nrow(filter(res, (FC4_overlaps == 1) & (dir_FC == "-")))
  FC2_up[i] <- nrow(filter(res, (FC2_overlaps == 1) & (dir_FC == "+")))
  FC2_down[i] <- nrow(filter(res, (FC2_overlaps == 1) & (dir_FC == "-")))
  core_up[i] <- nrow(filter(res, (core_overlaps == 1) & (dir_FC == "+")))
  core_down[i] <- nrow(filter(res, (core_overlaps == 1) & (dir_FC == "-")))
  tmp <- c(FC4_up[i], FC4_down[i], FC2_up[i], FC2_down[i], core_up[i], core_down[i])
  permutations <- rbind(permutations, tmp)
}
colnames(permutations) <- c("FC2_up", "FC2_down", "FC4_up", "FC4_down", "core_up", "core_down")

write.table(permutations, outfile, quote = FALSE, row.names = FALSE)

# observed permutation results
obs_FC2_up <- nrow(filter(data, (FC2_overlaps == 1) & (dir_FC == "+")))
obs_FC2_down <- nrow(filter(data, (FC2_overlaps == 1) & (dir_FC == "-")))
obs_FC4_up <- nrow(filter(data, (FC4_overlaps == 1) & (dir_FC == "+")))
obs_FC4_down <- nrow(filter(data, (FC4_overlaps == 1) & (dir_FC == "-")))
obs_core_up <- nrow(filter(data, (core_overlaps == 1) & (dir_FC == "+")))
obs_core_down <- nrow(filter(data, (core_overlaps == 1) & (dir_FC == "-")))

pvalue_FC2_up <- ((nrow(permutations[permutations$FC2_up >= obs_FC2_up, ])+1)/(nsim + 1))
pvalue_FC2_down <- ((nrow(permutations[permutations$FC2_down <= obs_FC2_down, ])+1)/(nsim + 1))
pvalue_FC4_up <- ((nrow(permutations[permutations$FC4_up >= obs_FC4_up, ])+1)/(nsim + 1))
pvalue_FC4_down <- ((nrow(permutations[permutations$FC4_down <= obs_FC4_down, ])+1)/(nsim + 1))
pvalue_core_up <- ((nrow(permutations[permutations$core_up >= obs_core_up, ])+1)/(nsim + 1))
pvalue_core_down <- ((nrow(permutations[permutations$core_down <= obs_core_down, ])+1)/(nsim + 1))


summary <- data.frame()
for (i in colnames(permutations)){
    overlap <- as.factor(i)
    range <- as.integer(range(permutations[, c(i)]))
    pvalue <- as.numeric(get(paste("pvalue_", i, sep = "")))
    row <- c(overlap, range, pvalue)
    summary <- rbind(summary, row)
}
colnames(summary) <- c("overlap", "range_start","range_end", "pvalue")

write.table(summary, summaryfile, quote = FALSE, row.names = FALSE)



results <- permutations
#####  GRAPH PERMUTATIONS

up_2FC <- ggplot(results) +
            geom_bar(aes(x = FC2_up)) +
            geom_vline(xintercept = obs_FC2_up, color = "red", size=1)

down_2FC <- ggplot(results) +
              geom_bar(aes(x = FC2_down)) +
              geom_vline(xintercept = obs_FC2_down, color = "red", size=1)


up_4FC <- ggplot(results) +
            geom_bar(aes(x = FC4_up)) +
            geom_vline(xintercept = obs_FC4_up, color = "red", size=1)

down_4FC <- ggplot(results) +
              geom_bar(aes(x = FC4_down)) +
              geom_vline(xintercept = obs_FC4_down, color = "red", size=1)

up_core <- ggplot(results) +
              geom_bar(aes(x = core_up)) +
              geom_vline(xintercept = obs_core_up, color = "red", size=1)

down_core <- ggplot(results) +
              geom_bar(aes(x = core_down)) +
              geom_vline(xintercept = obs_core_down, color = "red", size=1)

# see combined fig example http://www.sthda.com/english/articles/32-r-graphics-essentials/126-combine-multiple-ggplots-in-one-graph/
figure <- ggarrange(up_2FC, down_2FC, up_4FC, down_4FC, up_core, down_core,
                    labels = c("A", "B", "C", "D", "E", "F"),
                    ncol = 2, nrow = 3)
figure


ggsave(graphout)





#### EXPERIMENT to see if you can iterate through names of columns using a basename

#rows <- c("FC2", "FC4", "core")

#summary <- data.frame()
#for (i in rows) {
#  OL <- eval(paste(i, "_overlaps", sep = ""))
#  obs_up <- eval(paste(i, "_obs_up", sep = ""))
#  obs_down <- eval(paste(i, "_obs_down", sep = ""))
#  up <- eval(paste(i, "_up", sep = ""))
#  down <- eval(paste(i, "_down", sep = ""))#
#  obs_up <- nrow(filter(data, (OL == 1) & (dir_FC == "+")))
#  obs_down <- nrow(filter(data, (OL == 1) & (dir_FC == "-")))
#  pvalue_up <- ((nrow(permutations[permutations$up >= obs_up, ])+1)/(nsim +1))
#  pvalue_down <- ((nrow(permutations[permutations$down <= obs_down, ])+1)/(nsim +1))
#  tmp <- nrow(permutations[permutations$up >= obs_up, ])
#  print(permutations$up)
#  print(tmp)
#  range_up <- range(permutations$up, na.rm = TRUE)
#  #range_down <- range(permutations$down, na.rm = FALSE)
#  row <- c(obs_up, range_up, pvalue_up, obs_down, range_down, pvalue_down)
#  summary <- rbind(summary, row)
#}

#range(permutations$FC2_up, na.rm = TRUE)

#colnames(summary) <- c("obs_up", "range_up", "pvalue_up", "obs_down", "range_down", "pvalue_down")
#summary
