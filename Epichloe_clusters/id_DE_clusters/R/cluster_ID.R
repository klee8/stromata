# randomise DE and see if you can still identify clusters

#install.packages("doParallel")
library(tidyverse)
library(data.table)


###############################################################
#           GET FILENAMES AND NUMBER OF SIMULATIONS NEEDED
###############################################################

# TESTRUNS:
basefilename <- "STR_PS"
functionfile <- "/media/kate/Massey_linux_onl/projects/stromata_analysis/Epichloe_clusters/id_DE_clusters/R/cluster_fun.R"
datadir <- "/media/kate/Massey_linux_onl/projects/results/stromata/Epichloe_clusters/rfmt_core_gene_sets/"
resdir <- "/media/kate/Massey_linux_onl/projects/results/stromata/Epichloe_clusters/id_DE_clusters/STR_PS"
args <- c(basefilename, functionfile, datadir, resdir)
setwd(resdir)


#args = commandArgs(trailingOnly=TRUE)
#dir.create("cluster_perm_results")

if (length(args)<1) {
  stop("At least one argument must be supplied: <base_filename>, <number of simulations [default = 10,000]>.n", call.=FALSE)
} else {
  basefilename <- args[1]
  functionfile <- args[2]
  datadir <- args[3]
  resdir <- args[4]
  infile = paste(datadir, "core_genes_", args[1], "_rfmt.txt", sep = "")
  cluster_res_file = paste("cluster_perm_results/", basefilename, "_cluster_results.txt", sep = "")
  gene_res_file = paste("cluster_perm_results/", basefilename, "_per_gene_cluster_results.txt", sep = "")
  perm_res_file = paste("cluster_perm_results/", basefilename, "_permuted_DE_cluster_results.txt", sep = "")
  sum_perm_file = paste("cluster_perm_results/", basefilename, "_permuted_cluster_summary.txt", sep = "")
  graphout_file = paste("cluster_perm_results/", basefilename, "_permuted_DE_clusters.pdf", sep = "")
}
checks <- c(basefilename, functionfile, datadir, resdir)
print(checks)

#setwd(paste(resdir, basefilename, sep = ""))
source(functionfile)

# read in reformatted data
df <- read.delim(infile, header = TRUE, sep = "\t")

# setup data to find clusters (order by contig and position, rm duplicates,
# bin DE into 0, 1 (>=2fold), 2 (>=4fold), identify groups of consecutive
# genes on the same contig with the same direction of fold change)
# total_groups <- 0
df <- setup_cluster_assessment(df)

# identify clusters, and list the genes/orthologs in them
clusters <- find_assessment_blocks_with_clusters(df)

# put cluster information back into main data frame
df <- genes_in_clusters(df, clusters)
#head(df[df$clust_len != 0,])

write.table(df, "all_genes_cluster_annotated.txt", quote = FALSE, row.names = FALSE)

# use cluster-annotated df to summarise cluster information
cluster_sum <- sum_clusters(df)
#head(cluster_sum)
write.table(cluster_sum, cluster_res_file, row.names = FALSE, quote = FALSE, sep = "\t")

# number of observed clusters
obs_clusters <- nrow(clusters)

# check for shared genes among clusters
df <- shared_cluster_genes(df)

# add SMURF backbone data
smurf_bkbone <- read.delim("/media/kate/Massey_linux_onl/projects/results/stromata/Epichloe_clusters/SMURF/SMURF_flat_results_by_orthologs.txt")
smurf <- smurf_bkbone[, c("ortho_group", "el_SMURF_backbone_gene_prediction", "fes_SMURF_backbone_gene_prediction", "typ_SMURF_backbone_gene_prediction")]
smurf <- smurf[!is.na(smurf$el_SMURF_backbone_gene_prediction) | !is.na(smurf$fes_SMURF_backbone_gene_prediction) | !is.na(smurf$typ_SMURF_backbone_gene_prediction),]
smurf$bkbone <- paste(smurf$el_SMURF_backbone_gene_prediction, smurf$fes_SMURF_backbone_gene_prediction, smurf$typ_SMURF_backbone_gene_prediction, sep = ",")
smurf$bkbone <- gsub( "NA,", "", smurf$bkbone)
smurf$bkbone <- gsub( ",NA", "", smurf$bkbone)
smurf$bkbone <- as.factor(smurf$bkbone)
smurf <- smurf[!is.na(smurf$bkbone),]
smurf$bkbone <- as.factor(smurf$bkbone)
smurf$ortho_group <- factor(smurf$ortho_group, levels=levels(df$orthogroup)) 

df$smurf_bk <- sapply(df$orthogroup, function(x) ifelse(x %in% smurf$ortho_group, smurf[smurf$ortho_group == x, c("bkbone")], NA))
temp <- smurf$bkbone[df$smurf_bk] 
df$smurf_bk <- temp
df
# write out per-gene results data frame
write.table(data.frame(df), gene_res_file, row.names = FALSE, quote = FALSE, sep = "\t")

  # write out python scripts to graph clusters
dir.create("cluster_perm_results/Rgraphs")
make_pygraphs(df)
