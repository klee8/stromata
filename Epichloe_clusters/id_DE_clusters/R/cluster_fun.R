

##############################################################################3
#                Setup to identify clusters
###############################################################################

#' remove any duplicate rows
rm_dup_rows <- function(df) {
  df <- df %>% distinct
  return(df)
} 

#' set DE at 1 for log2FC >= 1 and svalue <= 0.005, set DE at 2 for log2FC >= 2 and svalue <= 0.005
#' 
bin_DE <- function(df) { 
  df <- df %>% mutate(DE = ifelse(((log2fc >= 2) | (log2fc <= -2)) & (svalue_1 <= 0.005), 2,
                                  ifelse(((log2fc >= 1) | (log2fc <= -1)) & (svalue_1 <= 0.01), 1 , 0)))
  return(df)
}

#' flag direction of fold change
dir_foldchange <- function(df) {
  df <- df %>% mutate(dir_FC = ifelse((log2fc < 0), "-", "+"))
  return(df)
}

#' order dataframe
order_contig_start <- function(df) {
  df <- df[order(df$species, df$contig, df$start),]
  return(df)
}

#' identify groups of genes, changing group number when DE changes direction or contig changes
#' used by identify_assessment_blocks function only
DE_change_groups <- function(contig, lag_contig, direction_FC, lag_direction_FC) {
  if (is.na(direction_FC) | is.na(lag_direction_FC) | is.na(contig) | is.na(lag_contig)) { total_groups <<- total_groups + 1 }
  else if (contig != lag_contig) {total_groups <<- total_groups + 1}
  else if (direction_FC != lag_direction_FC) {total_groups <<- total_groups + 1}
  return(total_groups)
}


#' identifies blocks of genes on the same contig with the same DE fold change direction
#' takes in df from setup_cluster_assessment
identify_assessment_blocks <- function(df) {
  total_groups <<- 0
  df$temp_contig <- lag(df$contig)
  df$temp_dirFC <- lag(df$dir_FC)
  df$change <- mapply(DE_change_groups, df$contig, df$temp_contig, df$dir_FC, df$temp_dirFC)
  df$temp_contig <- NULL
  df$temp_dirFC <- NULL
  return(df)
}


#' setup dataframe for cluster assessment
#' remove duplicate rows, bin DE into those above 2 and 4 fold
#' flag the direction of fold change (positive or negative)
#' order the data rame and flag groups of genes for assessment
#' (new group started when DE changes or contig changes)
#' 
setup_cluster_assessment <- function(df) {
  df <- order_contig_start(df)
  df <- rm_dup_rows(df)
  df <- bin_DE(df)
  df <- dir_foldchange(df)
  df <- identify_assessment_blocks(df)
  return(df)
}



##############################################################################################
#       identify clusters in observed data (includes gene, ortho output)
##############################################################################################


#' identify groups of genes with a cluster
#' takes groups of genes identified in DE_change_groups
#' collapses each set of gene DE results into a string
#' assess with regex for potential clusters
find_assessment_blocks_with_clusters <- function(df) {
  clusters <- data.frame()
  my_orths <- as.character()
  my_genes <- as.character()
  my_clust_num <- 0
  for (group in 1:max(df$change)) {
    DEstring <- paste(df[df$change == group, c("DE")], collapse = "")
    spp <- as.character(unique(df[df$change == group, c("species")]))
    DEorths <- df[df$change == group, c("orthogroup")]
    DEgenes <- df[df$change == group, c("gene_id")]
    
    # for clusters with 2 fold change DE, can have a gap, must have one 4 fold change DE
    #if ((grepl("[1,2]+[0]{0,1}[1,2]+", DEstring, perl = TRUE)) && (grepl("2", DEstring, perl = TRUE))) {  
    
    # for clusters with 2 fold change DE, can have a gap
    if (grepl("[1,2]+[0]{0,1}[1,2]+", DEstring, perl = TRUE)) {
      temp <- (gregexpr("[1,2]+[0]{0,1}[1,2]+", DEstring, perl = TRUE))
      for (i in temp) {
        pos <- temp[1][[1]][1]
        len <- (attr(temp[1][[1]], "match.length"))
        clust <- substr(DEstring, pos, (pos + len -1))
        if ((len > 3) || ( (len == 3) & (grepl('0', clust, perl = TRUE) == FALSE) ) ){
          my_clust_num <- my_clust_num + 1
          orths <- paste(DEorths[pos:(pos + len -1)], collapse = ",")
          genes <- paste(DEgenes[pos:(pos + len -1)], collapse = ",")
          my_orths <- c(my_orths, orths, stringsAsFactors = FALSE)
          my_genes <- c(my_genes, genes, stringsAsFactors = FALSE)
          row <- c(spp, my_clust_num, group, clust, pos, len, orths, genes)
          clusters <- rbind(clusters, row, stringsAsFactors = FALSE)
        }
      }
    }
  }
  colnames(clusters) <- c("species", "cluster_number", "group", "cluster_DE", "pos", "len", "orths", "genes")
  return(clusters)
}


#' identify exact locations of clusters in data
#' takes in data from find_assessment_blocks_with_clusters
#' flags genes that are in a cluster (and not just in a group with a cluster)
genes_in_clusters <- function(df, clusters) {
  prev <- 0
  flag <- 1
  inclust <- as.character()
  for (x in df$change) {
    if (x != prev) { flag <-  1 }
    if (x %in% clusters$group) {
      start <- as.numeric(clusters[clusters$group == x, c("pos")])
      end <- as.numeric(clusters[clusters$group == x, c("pos")]) + as.numeric(clusters[clusters$group == x, c("len")]) - 1
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
  df$clust_number <- sapply(df$in_cluster, function(x) ifelse( x %in% clusters$group, clusters[clusters$group == x, c("cluster_number")], 0 ))
  df$clust_pos <- sapply(df$in_cluster, function(x) ifelse( x %in% clusters$group, clusters[clusters$group == x, c("pos")], 0 ))
  df$clust_len <- sapply(df$in_cluster, function(x) ifelse( x %in% clusters$group, clusters[clusters$group == x, c("len")], 0 ))
  df$clust_cig <- sapply(df$in_cluster, function(x) ifelse( x %in% clusters$group, clusters[clusters$group == x, c("cluster_DE")], 0 ))
  
  ### remove superflous columns
  return(df)
}


#' run 
annotate_genes_with_cluster_info <- function(df) {
  clusters <- find_assessment_blocks_with_clusters(df)
  df <- genes_in_clusters(df, clusters)
}


#####################################################################################
#               Summarising observed clusters 
######################################################################################

# Gives nucleotide distance and number of genes to previous cluster 
# in case it may be necessary to relax the number of gaps allowed

# Flags when the cluster is the first/last gene on a contig
# cluster may be broken by sequencing (you may find the other part on 
# another contig or it may be missing)


#' calculates the number of genes since the begining of a contig
#' or the last cluster
#' 
#' note only called from within 'get_contig_gene_position'
contig_gene_pos <- function(contig, lag_contig, gene_pos){
  if (is.na(contig) | is.na(lag_contig)) { gene_pos <- 0}
  else if (contig != lag_contig) { gene_pos <- 0 }
  else { gene_pos <- gene_pos + 1 }
  return(gene_pos)
}


#' creates a column with the number of genes since the begining 
#' of a contig or since the last cluster
get_contig_gene_pos <- function(df) {
  # contig column present
  # check gene order
  gene_pos <- 0
  df$lag_contig <- lag(df$contig)
  df$gene_pos <- mapply(contig_gene_pos, df$contig, df$lag_contig, gene_pos)
  df$lag_contig <- NULL
  return(df)
}


#' calculates the physical distance in nuclotides since the begining 
#' of the contig or since the last cluster
#' calculates the number of genes between a cluster and the begining 
#' of the contig or since the last cluster
distance_between_clusters <- function(df) {
  genecount <- 0
  flag <- 0
  clustnum <- 0
  laststop <- 0
  df$prev_contig <- lag(df$contig)
  df$prev_contig[1] <- "Eel_0001"
  for (row in 1:nrow(df)) {
    # flag changes in contig
    if (df$contig[row] != df$prev_contig[row]) {
      flag <- 0
      genecount <- 0 
      laststop <- 0
      #rint(paste("start of contig here:", flag, genecount, laststop, sep = " "))
    }
    # flag end of cluster
    if((flag >= 1) &  (df$clust_cig[row] == 0) ) { 
      flag <- 0
      genecount <- 0
    }
    # within cluster
    if (df$clust_cig[row] != 0 ) {
      flag <- flag + 1
      # for the first position only
      if (flag == 1) {
        df$distance[row] <- df$start[row] - laststop
        df$gene_dist[row] <- genecount 
        df$clustnum[row] <- clustnum + 1
        clustnum <- clustnum + 1
        #print(paste("##########   counting cluster", clustnum, "flag: ", flag, "startpos: ", df$start[row], 
        #            "lastpos: ", laststop, "distance: ", df$distance[row], "contig: ", df$contig[row], 
        #            "gene distance: ", df$gene_dist[row], "cluster number:",  df$clustnum[row], sep = " "))
      }
      else {
        df$distance[row] <- NA
        df$gene_dist[row] <- NA
        df$clustnum[row] <- clustnum
      }
      laststop <- df$stop[row]
    }
    if(flag == 0) {
      genecount <- genecount + 1
      df$distance[row] <- NA
      df$gene_dist[row] <- NA
      df$clustnum[row] <- NA
    }
    
  }
  df$prev_contig <- NULL
  return(df)
}


#' flag the end of the last gene in a contig
#' use this prior to summarising data
end_genes <- function(df) {
  df <- df %>% group_by(contig) %>%
    mutate(contig_end = max(stop)) 
  return(df)
}

#' summarises cluster information in a data frame with one row per cluster
#' 
final_res <- function(df) {
  summary <- df %>% filter(clust_cig != 0) %>% group_by(change) %>%
  summarise(species = unique(species),
            cluster = unique(clustnum), 
            contig = unique(contig), 
            cluster_size = unique(clust_len), 
            start = min(start), 
            end = max(stop), 
            dist_to_last_cluster = unique(distance[!is.na(distance)]), 
            genes_between_clusters = unique(gene_dist[!is.na(gene_dist)]), 
            genes = paste(gene_id, collapse = ","), 
            ortholog_numbers = paste(orthogroup, collapse = ","),
            start_of_contig = ifelse(min(gene_pos) == 0 , TRUE, FALSE), 
            end_of_contig = ifelse(max(gene_pos) == unique(contig_end), TRUE, FALSE))
  return(summary)
}

#' Takes in raw cluster results, summarises data (including gene distances, nucleotide distance, 
#' cluster number) and outputs a wide file with one row per cluster
sum_clusters <- function(df) {
  df <- get_contig_gene_pos(df)
  df <- distance_between_clusters(df)
  df <- end_genes(df)
  clusters <- final_res(df)
  return(clusters)
}



########################################################################################
#     Permutation of DE (within species) to see how many clusters are found by chance
########################################################################################

#' Permutation of DE by species 
#' count how many clusters found
#' 
permute_clusters <- function(df, permutation_counts_outfile) {
  ## standard approach: scramble response value
  write("number_of_clusters", permutation_counts_outfile)
  bdat <- data.frame()
  for (spp in unique(df$species)) {
    perm <- sample(nrow(df[df$species == spp,]))
    tmp <- transform(df[df$species == spp,], DE = DE[perm])
    tmp <- transform(tmp[tmp$species == spp,], log2fc = log2fc[perm])
    tmp <- transform(tmp[tmp$species == spp,], dir_FC = dir_FC[perm])
    tmp <- transform(tmp, DE2 = ifelse((DE == 1) | (DE == 2), 1, 0), DE4 = ifelse((DE == 2), 1, 0))
    bdat <- rbind(bdat, tmp)
  }
  
  # order dataframe
  bdat <- bdat[order(bdat$species, bdat$contig, bdat$start),]
  
#  # flag direction of fold change
  bdat <- bdat %>% mutate(dir_FC = ifelse((log2fc < 0), "-", "+"))
  # flag when direction of DE or contig (including across species) changes
  bdat$temp_contig <- lag(bdat$contig)
  bdat$temp_dirFC <- lag(bdat$dir_FC)
  bdat$change <- mapply(DE_change_groups, bdat$contig, bdat$temp_contig, bdat$dir_FC, bdat$temp_dirFC)
  
  # check number of clusters found in this permutation
  tmp_clusters <- data.frame()
  for (group in 1:max(bdat$change)) {
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
  write.table(nrow(tmp_clusters), permutation_counts_outfile , quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t", append = TRUE)
  return(nrow(tmp_clusters)) 
}


summarise_perm <- function(permutations, observed, outfile) {
  pvalue <- 2*mean(permutations>=observed)
  write(paste("pvalue:", pvalue), outfile)
  range <- as.integer(range(permutations))
  write(paste("range:", range[1], range[2]), outfile, append = TRUE)
}

###############################################################################
#           Check if genes are shared among clusters/species
###############################################################################

#' reduce data frame to genes in clusters
#' identify which genes are shared with other clusters (usually across species)
#' 
shared_cluster_genes <- function(df) {
    df <- df[df$in_cluster != 0, ]
    df$shared <- sapply(df$orthogroup, function(x) sum(df$orthogroup == x))
    df$spp_groups <- sapply(df$orthogroup, function(x) paste(df[df$orthogroup == x, c("species")], collapse = ","))
    df$cluster_groups <- sapply(df$orthogroup, function(x) paste(df[df$orthogroup == x, c("clust_number")], collapse = ","))
    return(df)
}





################################################################################
#              write python function to draw clusters
################################################################################




#' generate a python script to make svg of cluster
#' 
#' 
#' 
#'
cluster_graph_pyscript <- function(df, counter) {
  
  header <- "#!usr/bin/python
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from reportlab.lib.units import cm
  
gdd = GenomeDiagram.Diagram('graph')
gdt_features = gdd.new_track(1, greytrack=False)
gds_features = gdt_features.new_set()
  
# add features
"

  # open file to write to 
  spp = unique(df$species)
  write(header, paste("cluster_perm_results/Rgraphs/graph", counter, spp, "py", sep = "."))
  
  # write out parameters for each gene
  for (orth in df$orthogroup) {
    temp <- df[df$orthogroup == orth,]
    if (temp$strand == "+") {
      strand = "+1"
      beg = temp$start
      end = temp$stop
      lab_pos = "start"
      angle = 0
    }
    else {
      strand = "-1"
      beg = temp$start
      end = temp$stop
      lab_pos = "end"
      angle = 180
    }
    if ((temp$shared ==1) && (temp$DE == 0)) { colour = "gainsboro"}
    if ((temp$shared ==1) && (temp$DE != 0)) { colour = "grey"}
    if ((temp$shared !=1) && (temp$DE == 0)) { colour = "indianred"}
    if ((temp$shared !=1) && (temp$DE != 0)) { colour = "darkred"}
    line1 = paste("feature = SeqFeature(FeatureLocation(", beg, ",",  end, "), strand=", strand, ")", sep = "")
    write(line1, paste("cluster_perm_results/Rgraphs/graph", counter, spp, "py", sep = "."), append = TRUE)
    line2 = paste("gds_features.add_feature(feature, name=\"", orth, "\", label=\"True\", color=\"", colour, "\", label_size=10, label_position=\"", lab_pos ,"\", label_angle=", angle, ", sigil=\"BIGARROW\")", sep = "")
    write(line2, paste("cluster_perm_results/Rgraphs/graph", counter, spp, "py", sep = "."), append = TRUE)
  }
  
  # formatting
  clust_start = min(min(df$start), min(df$stop))
  clust_end = max(max(df$start), max(df$stop))
  length = (abs(clust_end - clust_start))/1000
  line3 = paste("gdd.draw(format='linear', pagesize=(", length, "*cm,4*cm), fragments=1, start=", clust_start, ", end=", clust_end, ")", sep = "")
  write(line3, paste("cluster_perm_results/Rgraphs/graph", counter, spp, "py", sep = "."), append = TRUE)
  line4 = paste("gdd.write(\"", paste( spp, "diagram", counter,"svg", sep = "."), "\", \"SVG\")", sep = "")
  write(line4, paste("cluster_perm_results/Rgraphs/graph", counter, spp, "py", sep = "."), append = TRUE)
  
  counter <<- counter + 1
}  

#' group df by cluster and run cluster_graph_pyscript on each one
#' 
make_pygraphs <- function(df) {
  df <- df[df$in_cluster != 0, ]
  counter <<- 1
  clustdf <- shared_cluster_genes(df)
  for (clust in unique(df$in_cluster)) {
   cluster_graph_pyscript(df[df$in_cluster == clust, ], counter)
  }
}

