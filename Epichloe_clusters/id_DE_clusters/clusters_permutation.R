# randomise DE and see if you can still identify clusters


library(tidyverse)


# read in reformatted data
df <- read.delim("../core_gene_sets/core_genes_INF_PS_rfmt.txt", header = TRUE, sep = "\t")

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

# identify clusters in each spp by collapsing each group of gene DE results into a string 
# assess with regex
clusters <- data.frame()
for (group in 1:counter) {
  DEstring <- paste(df[df$change == group, c("DE")], collapse = "")
  if (grepl("[1,2]+[0]{0,1}[1,2]+", DEstring, perl = TRUE)) {
    temp <- (gregexpr("[1,2]+[0]{0,1}[1,2]+", DEstring, perl = TRUE))
    for (i in temp) {
      pos <- temp[1][[1]][1]
      len <- (attr(temp[1][[1]], "match.length"))
      clust <- substr(DEstring, pos, (pos + len -1))
      if ((len > 3) || ( (len == 3) & (grepl('0', clust, perl = TRUE) == FALSE) ) ){
        #print(paste(group, DEstring, pos, len, sep = " "))
        row <- c(as.numeric(group), as.numeric(clust), as.numeric(pos), as.numeric(len))
        clusters <- rbind(clusters, row)
      } 
    }
  }
}

# add cluster information to main dataframe
colnames(clusters) <- c("group", "cluster_DE", "pos", "len")
df$clust_pos <- sapply(df$change, function(x) ifelse( x %in% clusters$group, clusters[clusters$group == x, c("pos")], 0 ))
df$clust_len <- sapply(df$change, function(x) ifelse( x %in% clusters$group, clusters[clusters$group == x, c("len")], 0 ))
df$clust_cig <- sapply(df$change, function(x) ifelse( x %in% clusters$group, clusters[clusters$group == x, c("cluster_DE")], 0 ))


write.table(df, "cluster_id_in_R.txt", row.names = FALSE, quote = FALSE, sep = "\t")
