# randomise DE and see if you can still identify clusters

#install.packages("doParallel")
library(doParallel)
library(tidyverse)
library(data.table)


###############################################################
#           GET FILENAMES AND NUMBER OF SIMULATIONS NEEDED
###############################################################

#resdir <- "/media/kate/Massey_linux_onl/projects/results/stromata/Epichloe_clusters/id_DE_clusters/STR_PS/"
#scriptdir <- "/media/kate/Massey_linux_onl/projects/stromata_analysis/Epichloe_clusters/id_DE_clusters/"
#setwd(resdir)
#args = c("STR_PS", 10, "cluster_fun.R")

args = commandArgs(trailingOnly=TRUE)
dir.create("cluster_perm_results")

if (length(args)<1) {
  stop("At least one argument must be supplied: <base_filename>, <number of simulations [default = 10,000]>.n", call.=FALSE)
} else {
  basefilename <- args[1]
  nsim <- as.numeric(args[2])
  functionfile <- args[3]
  datadir <- args[4]
  resdir <- args[5]
  infile = paste(datadir, "core_genes_", args[1], "_rfmt.txt", sep = "")
  cluster_res_file = paste("cluster_perm_results/", basefilename, "_cluster_results.txt", sep = "")
  perm_res_file = paste("cluster_perm_results/", basefilename, "_permuted_DE_cluster_results.txt", sep = "")
  sum_perm_file = paste("cluster_perm_results/", basefilename, "_permuted_cluster_summary.txt", sep = "")
  graphout_file = paste("cluster_perm_results/", basefilename, "_permuted_DE_clusters.pdf", sep = "")
}

setwd(resdir)
source(functionfile)

# read in reformatted data
df <- read.delim(infile, header = TRUE, sep = "\t")

# setup data to find clusters (order by contig and position, rm duplicates,
# bin DE into 0, 1 (>=2fold), 2 (>=4fold), identify groups of consecutive
# genes on the same contig with the same direction of fold change)
#total_groups <- 0
df <- setup_cluster_assessment(df)

# identify clusters, and list the genes/orthologs in them
clusters <- find_assessment_blocks_with_clusters(df)

# put cluster information back into main data frame
df <- genes_in_clusters(df, clusters)
head(df[df$clust_len != 0,])

# use cluster-annotated df to summarise cluster information
cluster_sum <- sum_clusters(df)
head(cluster_sum)
write.table(cluster_sum, cluster_res_file, row.names = FALSE, quote = FALSE, sep = "\t")

# number of observed clusters
obs_clusters <- nrow(clusters)

# write out python scripts to graph clusters
make_pygraphs(df)

# run permutation function with foreach parellisation
n.cores <- detectCores()
registerDoParallel(n.cores)
num_clusters <- foreach(k = 1:nsim, .combine = c) %dopar% permute_clusters(df, perm_res_file)

# summarise permutations
summarise_perm(num_clusters, obs_clusters, sum_perm_file)

# plot the permutations as a histogram
permutation_plot <- ggplot() +
  geom_bar(aes(x = num_clusters)) +
  geom_vline(xintercept = obs_clusters, color = "red", size=1)
ggsave(graphout_file, plot = permutation_plot)
