# Kate Lee 2019
# functions for interacting with topGO


#' Keep results from fisher test each time it is run by topGO
#' 
#' Sets up file if it doesn't exist
#' NB it will append the file if it does exist
#' Takes results from resultFisher data in topGO
#' outputs main results to a file called "fisher_results.txt"
keep_fish <- function(resultFisher) {
  if (!file.exists("fisher_results.txt")){
    names <- paste(c("Description", "GO_terms_scored", "terms_with_pvalue>0.01", 
                     "annotated_genes", "significant_annotated_genes", 
                     "min_genes_associated_with_a_GO", "nontrivial_nodes"), collapse = "\t") 
    write.table(names, "fisher_results.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  descript <- gsub(" \n.*", "", resultFisher@description)
  scored <- length(resultFisher@score)
  sig <- length(resultFisher@score[resultFisher@score < 0.01])
  ann <- unname(resultFisher@geneData[1])
  sig_ann <- unname(resultFisher@geneData[2])
  nodesize <- unname(resultFisher@geneData[3])
  sig_terms <- unname(resultFisher@geneData[4])
  row <- paste(c(descript, scored, sig, ann, sig_ann, nodesize, sig_terms), collapse = "\t")
  write.table(row, "fisher_results.txt", quote = FALSE, row.names = FALSE, append = TRUE, col.names = FALSE)
}


#' Keep results from Kolmogorov-Smirnov test in topGO
#' 
#' Sets up file if it doesn't exist
#' NB it iwll append the file if it does exist
#' Takes in results form resultKS from topGO
#' outputs main results to file
keep_ks <- function(resultKS) {
  if (!file.exists("KS_results.txt")) {
    names <- paste(c("Description", "GO_terms_scored", "terms_with_pvalue>0.01", 
                     "annotated_genes", "significant_annotated_genes", 
                     "min_genes_associated_with_a_GO", "nontrivial_nodes"), collapse = "\t") 
    write.table(names, "KS_results.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  descript <- gsub(" \n.*", "", resultKS@description)
  scored <- length(resultKS@score)
  sig <- length(resultKS@score[resultKS@score < 0.01])
  ann <- unname(resultKS@geneData[1])
  sig_ann <- unname(resultKS@geneData[2])
  nodesize <- unname(resultKS@geneData[3])
  sig_terms <- unname(resultKS@geneData[4])
  row <- paste(c(descript, scored, sig, ann, sig_ann, nodesize, sig_terms), collapse = "\t")
  write.table(row, "KS_results.txt", quote = FALSE, row.names = FALSE, append = TRUE, col.names = FALSE)
}


#' Checks the number of genes input to topGO
#' 
#' outputs the size of the geneUniverse, the genelist
#' and the number of genes that are annotated (usually the entire geneUniverse)
check_topgo_input <- function(geneUniverse, geneList, geneID2GO){ 
  if (!file.exists("inputfile_stats.txt")) {
    names <- paste(c("total_genes", "selected_goi", "genes_with_GO_ann"), collapse = "\t")
    write.table(names, "inputfile_stats.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  row <- paste(c(length(geneUniverse), length(geneList[geneList ==1]), length(geneID2GO)), collapse = "\t")
  write.table(row, "inputfile_stats.txt", quote = FALSE, row.names = FALSE, append = TRUE, col.names = FALSE)
}

#' Takes in Pannzer GO term list and converts for topGO input
#' 
#' NB you still have to save this to a file to be read in by topGO readMappings()
#' NB make sure that the gene names match the other topGO input files 
pannzer_to_topgo <- function(pannzer_GO_file) {
  pannzer_GO_file %>% 
    group_by(qpid) %>% 
    summarise(go_terms = paste(sapply(as.character(goid), function(x) paste("GO:", x, sep = "")), collapse = ", "))
}
  
