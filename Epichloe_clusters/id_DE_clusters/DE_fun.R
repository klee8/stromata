


#' add column which sets DE at 1 for log2FC >= 1 and svalue <= 0.005, set DE at 2 for log2FC >= 2 and svalue <= 0.005
#' Requires a dataframe with log2fc output (e.g. from DESeq) per gene
DE_categories <- function(data) {
  res <- data %>% mutate(DE = ifelse((log2fc >= 2) & (svalue_2 <= 0.005), 2,
                                    ifelse((log2fc >= 1) &(svalue_1 <= 0.005), 1, 
                                           ifelse((log2fc <= -2) &(svalue_2 <= 0.005), -2, 
                                                  ifelse((log2fc <= -1) &(svalue_1 <= 0.005), -1, 0)))))

}

#' Returns total genes and total number of DE genes per species
#' Requries df with 'species' and 'DE' columns per gene
species_raw_DE_counts <- function(data) {raw_counts <- data %>% group_by(species) %>% 
  summarise(total = n(), 
            total_DE2_up = length(which(DE == 1)) + length(which(DE == 2)), 
            total_DE4_up = length(which(DE == 2)),
            total_DE2_down = length(which(DE == -1)) + length(which(DE == -2)),
            total_DE4_down = length(which(DE == -2)))

return(raw_counts) 
}

#' find DE overlaps across all species (core_gene_set)
#' Requries df with 'species' and 'DE' columns per gene
DE_overlaps <- function(data) {
  max_spp <- (length(unique(data$species)))
  res <- data %>% 
    group_by(orthogroup) %>%
    summarise(FC2_up = ifelse((length(unique(species)) == max_spp)   & (!0 %in% DE) & (length(unique(dir_FC)) == 1) & (sum(DE) >= max_spp) , 1, 0),
              FC4_up = ifelse((length(unique(species)) == max_spp)   & (!0 %in% DE) & (length(unique(dir_FC)) == 1) & (sum(DE) == 2*max_spp), 1, 0),
              core_up = ifelse((length(unique(species)) == max_spp)  & (!0 %in% DE) & (length(unique(dir_FC)) == 1) & (sum(DE) >= max_spp + 1), 1, 0),
              FC2_down = ifelse((length(unique(species)) == max_spp) & (!0 %in% DE) & (length(unique(dir_FC)) == 1) & (sum(DE) <= -1*max_spp), 1, 0),                     
              FC4_down = ifelse((length(unique(species)) == max_spp) & (!0 %in% DE) & (length(unique(dir_FC)) == 1) & (sum(DE) == -2*max_spp), 1, 0),
              core_down = ifelse((length(unique(species)) == max_spp)& (!0 %in% DE) & (length(unique(dir_FC)) == 1) & (sum(DE) <= -1*max_spp - 1), 1, 0))   
  return(res)
}


# count up total number of each kind of overlap
total_overlaps <- function(res, max_spp){
  TOL_FC4_up <- nrow(filter(res, (FC4_up == 1)))
  TOL_FC4_down <- nrow(filter(res, (FC4_down == 1))) 
  TOL_FC2_up <- nrow(filter(res, (FC2_up == 1))) 
  TOL_FC2_down <- nrow(filter(res, (FC2_down == 1))) 
  TOL_core_up <- nrow(filter(res, (core_up == 1))) 
  TOL_core_down <- nrow(filter(res, (core_down == 1)))  
  tmp <- c(TOL_core_up, TOL_FC2_up, TOL_FC4_up, TOL_core_down, TOL_FC2_down, TOL_FC4_down )
  return(tmp)
}

# permutation test
DE_permutation <- function(data, outfile) {
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