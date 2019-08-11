# randomise DE and see if you can still identify clusters


library(tidyverse)

  
args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
  stop("At least one argument must be supplied: <base_filename>, <number of simulations [default = 10,000]>.n", call.=FALSE)
} else {
  infile = paste("../../core_gene_sets/core_genes_", args[1], "_rfmt.txt", sep = "") #| "../core_gene_sets/core_genes_INF_PS_rfmt.txt"
  outfile1 = paste(args[1], "_cluster_results.txt", sep = "") #| "observed_clusters.txt"
  outfile2 = paste(args[1], "_permuted_DE_cluster_results.txt", sep = "")   #| "permuted_data.txt"
  graphout = paste(args[1], "_permuted_DE_clusters.pdf", sep = "")
  nsim = args[2] #|| 10000
}

print(infile)
print(outfile1)
print(outfile2)
print(graphout)
print(nsim)

#infile <- "../core_gene_sets/core_genes_INF_PS_rfmt.txt"
#outfile1 <- "INF_PS_observed_clusters.txt"
#outfile2 <- "INF_PS_permuted_DE_cluster_numbers.txt"

# read in reformatted data
df <- read.delim(infile, header = TRUE, sep = "\t")

# remove any duplicate rows
df <- df %>% distinct

# set DE at 1 for log2FC >= 1 and svalue <= 0.005, set DE at 2 for log2FC >= 2 and svalue <= 0.005
df <- df %>% mutate(DE = ifelse(((log2fc >= 2) | (log2fc <= -2)) & (svalue_1 <= 0.005), 2, 
                                        ifelse(((log2fc >= 1) | (log2fc <= -1)) & (svalue_1 <= 0.01), 1 , 0)))

# flag direction of fold change
df <- df %>% mutate(dir_FC = ifelse((log2fc < 0), "-", "+"))

# order dataframe
df <- df[order(df$species, df$contig, df$start),]

# identify groups of genes, changing group number when DE changes direction or contig changes
changes <- function(a, b, x, y) {
  if (is.na(y) | is.na(x) | is.na(a) | is.na(b)) {counter <<- counter + 1}
  else if (x != y) {counter <<- counter + 1}     
  else if (a != b) {counter <<- counter + 1}
  return(counter)
}
counter <- 0
df$temp_contig <- lag(df$contig)
df$temp_dirFC <- lag(df$dir_FC)
df$change <- mapply(changes, df$contig, df$temp_contig, df$dir_FC, df$temp_dirFC)
df$temp_contig <- NULL
df$temp_dirFC <- NULL

# identify clusters in each spp by collapsing each group of gene DE results into a string 
# assess with regex
clusters <- data.frame()
#my_orths <- as.character()
#my_genes <- as.character()
for (group in 1:counter) {
  DEstring <- paste(df[df$change == group, c("DE")], collapse = "")
  #DEorths <- df[df$change == group, c("orthogroup")]  
  #DEgenes <- df[df$change == group, c("gene_id")] 
  if (grepl("[1,2]+[0]{0,1}[1,2]+", DEstring, perl = TRUE)) {
    temp <- (gregexpr("[1,2]+[0]{0,1}[1,2]+", DEstring, perl = TRUE))
    for (i in temp) {
      pos <- temp[1][[1]][1]
      len <- (attr(temp[1][[1]], "match.length"))
      clust <- substr(DEstring, pos, (pos + len -1))
      if ((len > 3) || ( (len == 3) & (grepl('0', clust, perl = TRUE) == FALSE) ) ){
        #orths <- paste(DEorths[pos:(pos +len -1)], collapse = ",")
        #genes <- paste(DEgenes[pos:(pos +len -1)], collapse = ",")
        #my_orths <- c(my_orths, orths)
        #my_genes <- c(my_genes, genes)
        row <- c(as.numeric(group), as.numeric(clust), as.numeric(pos), as.numeric(len))
        clusters <- rbind(clusters, row)
      } 
    }
  }
}
# add orth and gene lists
#clusters <- cbind(clusters, my_orths)
#clusters <- cbind(clusters, my_genes)
#colnames(clusters) <- c("group", "cluster_DE", "pos", "len", "orths", "genes")


colnames(clusters) <- c("group", "cluster_DE", "pos", "len")

# number of observed clusters
obs_clusters <- nrow(clusters)

#####      Flag the genes that are in a cluster (and not just in a group of genes which has a cluster)
prev <- 0
flag <- 1
inclust <- as.character()
for (x in df$change) {
  if (x != prev) { flag <-  1 }
  if (x %in% clusters$group) {
    start <- clusters[clusters$group == x, c("pos")]
    end <- clusters[clusters$group == x, c("pos")] + clusters[clusters$group == x, c("len")] - 1
    if ( (start <= flag) & (flag <= end ) ) {
      inclust <- c(inclust, x)
  }
    else { inclust <- c(inclust, 0) }
  } 
  else { inclust <- c(inclust, 0) }
  flag <- flag + 1
  prev <- x
}

# add cluster information to main dataframe  
df$in_cluster <- inclust
df$clust_pos <- sapply(df$in_cluster, function(x) ifelse( x %in% clusters$group, clusters[clusters$group == x, c("pos")], 0 ))
df$clust_len <- sapply(df$in_cluster, function(x) ifelse( x %in% clusters$group, clusters[clusters$group == x, c("len")], 0 ))
df$clust_cig <- sapply(df$in_cluster, function(x) ifelse( x %in% clusters$group, clusters[clusters$group == x, c("cluster_DE")], 0 ))

# write to file
write.table(df, outfile1, row.names = FALSE, quote = FALSE, sep = "\t")


# permute DE and check number of clusters found
num_clusters <- as.numeric()
for (i in 1:nsim) {
  ## standard approach: scramble response value
  bdat <- data.frame()
  for (spp in unique(df$species)) {
    perm <- sample(nrow(df[df$species == spp,]))
    tmp <- transform(df[df$species == spp,], DE = DE[perm])
    tmp <- transform(tmp, DE2 = ifelse((DE == 1) | (DE == 2), 1, 0), DE4 = ifelse((DE == 2), 1, 0))
    bdat <- rbind(bdat, tmp)
  }

  # check number of clusters found in this permutation
  tmp_clusters <- data.frame()
  for (group in 1:counter) {
    DEstring <- paste(bdat[bdat$change == group, c("DE")], collapse = "")
    if (grepl("[1,2]+[0]{0,1}[1,2]+", DEstring, perl = TRUE)) {
      temp <- (gregexpr("[1,2]+[0]{0,1}[1,2]+", DEstring, perl = TRUE))
      for (i in temp) {
        pos <- temp[1][[1]][1]
        len <- (attr(temp[1][[1]], "match.length"))
        clust <- substr(DEstring, pos, (pos + len -1))
        if ((len > 3) || ( (len == 3) & (grepl('0', clust, perl = TRUE) == FALSE) ) ){
          row <- c(as.numeric(group), as.numeric(clust), as.numeric(pos), as.numeric(len))
          tmp_clusters <- rbind(tmp_clusters, row)
        } 
      }
    }
  }
  num_clusters <- c(num_clusters, nrow(tmp_clusters))
}

write.table(num_clusters, outfile2, quote = FALSE, row.names = FALSE, sep = "\t")

ggplot() +
  geom_bar(aes(x = num_clusters)) +
  geom_vline(xintercept = obs_clusters, color = "red", size=1) 

ggsave(graphout)
