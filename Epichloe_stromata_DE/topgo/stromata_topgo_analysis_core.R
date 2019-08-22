# Kate Lee 2019
# Analysis of stromata DE genes with topgo (GO terms from pannzer)


#if (!requireNamespace("BiocManager", quietly=TRUE))+     install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install("topGO")
# BiocManager::install("Rgraphviz")
#BiocManager::install("RCurl")
library(topGO)
library(tidyverse)
library(Rgraphviz)
library(RCurl)


### Sanity check 
# have 88 up regulaged and 24 down regulated core genes.
# should have this many genes of interest for each species, but generally getting 75 up and 22 down regulated genes
# genesOfInterest list is the correct size
# geneList is not.... CONCLUSION - check that these genes are not missing in the GO annotation file => They are ;)

scriptdir <- "/media/kate/Massey_linux_onl/projects/stromata_analysis/Epichloe_stromata_DE/topgo/"
resdir <- "/media/kate/Massey_linux_onl/projects/results/stromata/Epichloe_stromata_DE/topGO/"

setwd(resdir)
dir.create("core")
setwd("core/")
source(paste(scriptdir, "topGO_fun.R", sep = "" ))

# read in differential expression data 
DE_all <- read.delim("../../3_DEseq2/STR_PS/core_gene_set_STR_PS_ann.txt", header = TRUE)


#####     FESTUCAE

# read in pannzer annotations and put into topGO datastructure geneID2GO
pannzer <- read.delim("../../Pannzer/E.festucae_E2368/E.festucae_E2368_pannzer_GO.txt", header = TRUE, sep = "\t")
name <- "festucae"
go_terms <- pannzer_to_topgo(pannzer)
go_terms$qpid <- gsub("-T1", "", go_terms$qpid )
write.table(go_terms, paste(name, "gene2GO.txt", sep = "_"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
geneID2GO <- readMappings(file = paste(name, "gene2GO.txt", sep = "_"))

# subset DE information
DE <- DE_all[, c( "festucae_gene_id", "core_up", "core_down", "top_4FC_up", "top_4FC_down", "top_2FC_up", "top_2FC_down")]



# set the gene universe (all genes considered) and the genes of interest
geneUniverse <- go_terms$qpid

# UPREGULATED

    genesOfInterest <- as.character(DE[(DE$core_up == 1), c("festucae_gene_id")])
    #fes_up_goi <- genesOfInterest

    # locations of genes of interest in the geneUniverse set in a named vector
    geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
    missing <- factor(as.integer(genesOfInterest %in% geneUniverse))
    names(missing) <- genesOfInterest
    names(geneList) <- geneUniverse
    check_topgo_input(geneUniverse, geneList, geneID2GO)
    
    #fes_up_list <- geneList[geneList == 1]
    #fes_up_missing <- names(missing[missing == 0])

    # create GOdata object
    myGOdata <- new("topGOdata", description="festucae_core_up_DE_BP", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher <- runTest(myGOdata, algorithm='weight01', statistic="fisher")
    allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
    allRes$test <- "festucae_core_up_DE_BP"
    write.table(allRes, "all_Fisher_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE,  sep = "	", append = TRUE)
    showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo ='all')
    printGraph(myGOdata, resultFisher, firstSigNodes = 5, fn.prefix = "festucae_core_up_BP", useInfo = "all", pdfSW = TRUE)
    keep_fish(resultFisher)
    resultKS <- runTest(myGOdata, algorithm = "weight01", statistic = "ks") 
    KSRes <- GenTable(myGOdata, classicKS = resultKS, orderBy = "resultKS", ranksOf = "classicFisher", topNodes = 10)
    KSRes$test <- "festucae_core_up_DE_BP"
    write.table(KSRes, "all_KS_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE,  sep = "	", append = TRUE)
    keep_ks(resultKS)
    
    myGOdata <- new("topGOdata", description="festucae_core_up_DE_CC", ontology="CC", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher <- runTest(myGOdata, algorithm='weight01', statistic="fisher")
    allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
    allRes$test <- "festucae_core_up_DE_CC"
    write.table(allRes, "all_Fisher_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "	", append = TRUE)
    showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo ='all')
    printGraph(myGOdata, resultFisher, firstSigNodes = 5, fn.prefix = "festucae_core_up_CC", useInfo = "all", pdfSW = TRUE)
    keep_fish(resultFisher)
    resultKS <- runTest(myGOdata, algorithm = "weight01", statistic = "ks") 
    KSRes <- GenTable(myGOdata, classicKS = resultKS, orderBy = "resultKS", ranksOf = "classicFisher", topNodes = 10)
    KSRes$test <- "festucae_core_up_DE_CC"
    write.table(KSRes, "all_KS_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "	", append = TRUE)
    keep_ks(resultKS)

    myGOdata <- new("topGOdata", description="festucae_core_up_DE_MF", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher <- runTest(myGOdata, algorithm='weight01', statistic="fisher")
    allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
    allRes$test <- "festucae_core_up_DE_MF"
    write.table(allRes, "all_Fisher_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "	", append = TRUE)
    showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo ='all')
    printGraph(myGOdata, resultFisher, firstSigNodes = 5, fn.prefix = "festucae_core_up_MF", useInfo = "all", pdfSW = TRUE)
    keep_fish(resultFisher)
    resultKS <- runTest(myGOdata, algorithm = "weight01", statistic = "ks") 
    KSRes <- GenTable(myGOdata, classicKS = resultKS, orderBy = "resultKS", ranksOf = "classicFisher", topNodes = 10)
    KSRes$test <- "festucae_core_up_DE_MF"
    write.table(KSRes, "all_KS_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE,  sep = "	", append = TRUE)
    keep_ks(resultKS)

# DOWNREGULATED
    
    genesOfInterest <- as.character(DE[(DE$core_down == 1), c("festucae_gene_id")])
    #fes_down_goi <- genesOfInterest    
    # locations of genes of interest in the geneUniverse set in a named vector
    geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
    names(geneList) <- geneUniverse
    check_topgo_input(geneUniverse, geneList, geneID2GO)
    
    # create GOdata object
    myGOdata <- new("topGOdata", description="festucae_core_down_DE_BP", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher <- runTest(myGOdata, algorithm='weight01', statistic="fisher")
    allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
    allRes$test <- "festucae_core_down_DE_BP"
    write.table(allRes, "all_Fisher_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "	", append = TRUE)
    showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo ='all')
    printGraph(myGOdata, resultFisher, firstSigNodes = 5, fn.prefix = "festucae_core_down_BP", useInfo = "all", pdfSW = TRUE)
    keep_fish(resultFisher)
    resultKS <- runTest(myGOdata, algorithm = "weight01", statistic = "ks") 
    KSRes <- GenTable(myGOdata, classicKS = resultKS, orderBy = "resultKS", ranksOf = "classicFisher", topNodes = 10)
    KSRes$test <- "festucae_core_down_DE_BP"
    write.table(KSRes, "all_KS_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE,  sep = "	", append = TRUE)
    keep_ks(resultKS)

    myGOdata <- new("topGOdata", description="festucae_core_down_DE_CC", ontology="CC", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher <- runTest(myGOdata, algorithm='weight01', statistic="fisher")
    allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
    allRes$test <- "festucae_core_down_DE_CC"
    write.table(allRes, "all_Fisher_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "	", append = TRUE)
    showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo ='all')
    printGraph(myGOdata, resultFisher, firstSigNodes = 5, fn.prefix = "festucae_core_down_CC", useInfo = "all", pdfSW = TRUE)
    keep_fish(resultFisher)
    resultKS <- runTest(myGOdata, algorithm = "weight01", statistic = "ks") 
    KSRes <- GenTable(myGOdata, classicKS = resultKS, orderBy = "resultKS", ranksOf = "classicFisher", topNodes = 10)
    KSRes$test <- "festucae_core_down_DE_CC"
    write.table(KSRes, "all_KS_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE,  sep = "	", append = TRUE)
    keep_ks(resultKS)

    myGOdata <- new("topGOdata", description="festucae_core_down_DE_MF", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher <- runTest(myGOdata, algorithm='weight01', statistic="fisher")
    allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
    allRes$test <- "festucae_core_down_DE_MF"
    write.table(allRes, "all_Fisher_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "	", append = TRUE)
    showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo ='all')
    printGraph(myGOdata, resultFisher, firstSigNodes = 5, fn.prefix = "festucae_core_down_MF", useInfo = "all", pdfSW = TRUE)
    keep_fish(resultFisher)
    resultKS <- runTest(myGOdata, algorithm = "weight01", statistic = "ks") 
    KSRes <- GenTable(myGOdata, classicKS = resultKS, orderBy = "resultKS", ranksOf = "classicFisher", topNodes = 10)
    KSRes$test <- "festucae_core_down_DE_MF"
    write.table(KSRes, "all_KS_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE,  sep = "	", append = TRUE)
    keep_ks(resultKS)


    
    
#####     TYPHINA
    
    # read in pannzer annotations and put into topGO datastructure geneID2GO
    
    pannzer<- read.delim("../../Pannzer/E.typhina_E8/E.typhina_E8_pannzer_GO.txt", header = TRUE, sep = "\t")
    name <- "typhina"
    go_terms <- pannzer_to_topgo(pannzer)
    go_terms$qpid <- gsub("-T1", "", go_terms$qpid )
    write.table(go_terms, paste(name, "gene2GO.txt", sep = "_"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    geneID2GO <- readMappings(file = paste(name, "gene2GO.txt", sep = "_"))
    
    # subset DE information
    DE <- DE_all[, c( "typhina_gene_id", "core_up", "core_down", "top_4FC_up", "top_4FC_down", "top_2FC_up", "top_2FC_down")]
    
    # set the gene universe (all genes considered) and the genes of interest
    geneUniverse <- go_terms$qpid
    
# UPREGULATED
    genesOfInterest <- as.character(DE[(DE$core_up == 1), c("typhina_gene_id")])
    #typ_up_goi <- genesOfInterest

    # locations of genes of interest in the geneUniverse set in a named vector
    geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
    names(geneList) <- geneUniverse
    check_topgo_input(geneUniverse, geneList, geneID2GO)
    
    # create GOdata object
    myGOdata <- new("topGOdata", description="typhina_core_up_DE_BP", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher <- runTest(myGOdata, algorithm='weight01', statistic="fisher")
    allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
    allRes$test <- "typhina_core_up_DE_BP"
    write.table(allRes, "all_Fisher_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "	", append = TRUE)
    showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo ='all')
    printGraph(myGOdata, resultFisher, firstSigNodes = 5, fn.prefix = "typhina_core_up_BP", useInfo = "all", pdfSW = TRUE)
    keep_fish(resultFisher)
    resultKS <- runTest(myGOdata, algorithm = "weight01", statistic = "ks") 
    KSRes <- GenTable(myGOdata, classicKS = resultKS, orderBy = "resultKS", ranksOf = "classicFisher", topNodes = 10)
    KSRes$test <-"typhina_core_up_DE_BP"
    write.table(KSRes, "all_KS_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE,   sep = "	", append = TRUE)
    keep_ks(resultKS)

    myGOdata <- new("topGOdata", description="typhina_core_up_DE_CC", ontology="CC", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher <- runTest(myGOdata, algorithm='weight01', statistic="fisher")
    allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
    allRes$test <- "typhina_core_up_DE_CC"
    write.table(allRes, "all_Fisher_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "	", append = TRUE)
    showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo ='all')
    printGraph(myGOdata, resultFisher, firstSigNodes = 5, fn.prefix = "typhina_core_up_CC", useInfo = "all", pdfSW = TRUE)
    keep_fish(resultFisher)
    resultKS <- runTest(myGOdata, algorithm = "weight01", statistic = "ks") 
    KSRes <- GenTable(myGOdata, classicKS = resultKS, orderBy = "resultKS", ranksOf = "classicFisher", topNodes = 10)
    KSRes$test <- "typhina_core_up_DE_CC"
    write.table(KSRes, "all_KS_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "	", append = TRUE)
    keep_ks(resultKS)

    myGOdata <- new("topGOdata", description="typhina_core_up_DE_MF", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher <- runTest(myGOdata, algorithm='weight01', statistic="fisher")
    allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
    allRes$test <- "typhina_core_up_DE_MF"
    write.table(allRes, "all_Fisher_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "	", append = TRUE)
    showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo ='all')
    printGraph(myGOdata, resultFisher, firstSigNodes = 5, fn.prefix = "typhina_core_up_MF", useInfo = "all", pdfSW = TRUE)
    keep_fish(resultFisher)
    resultKS <- runTest(myGOdata, algorithm = "weight01", statistic = "ks") 
    KSRes <- GenTable(myGOdata, classicKS = resultKS, orderBy = "resultKS", ranksOf = "classicFisher", topNodes = 10)
    KSRes$test <- "typhina_core_up_DE_MF"
    write.table(KSRes, "all_KS_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "	", append = TRUE)
    keep_ks(resultKS)

# DOWNREGULATED
    
    genesOfInterest <- as.character(DE[(DE$core_down == 1), c("typhina_gene_id")])
    #typ_down_goi <- genesOfInterest    
    # locations of genes of interest in the geneUniverse set in a named vector
    geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
    names(geneList) <- geneUniverse
    check_topgo_input(geneUniverse, geneList, geneID2GO)
    
    # create GOdata object
    myGOdata <- new("topGOdata", description="typhina_core_down_DE_BP", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher <- runTest(myGOdata, algorithm='weight01', statistic="fisher")
    allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
    allRes$test <- "typhina_core_down_DE_BP"
    write.table(allRes, "all_Fisher_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "	", append = TRUE)
    showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo ='all')
    printGraph(myGOdata, resultFisher, firstSigNodes = 5, fn.prefix = "typhina_core_down_BP", useInfo = "all", pdfSW = TRUE)
    keep_fish(resultFisher)
    resultKS <- runTest(myGOdata, algorithm = "weight01", statistic = "ks") 
    KSRes <- GenTable(myGOdata, classicKS = resultKS, orderBy = "resultKS", ranksOf = "classicFisher", topNodes = 10)
    KSRes$test <- "typhina_core_down_DE_BP"
    write.table(KSRes, "all_KS_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "	", append = TRUE)
    keep_ks(resultKS)

    myGOdata <- new("topGOdata", description="typhina_core_down_DE_CC", ontology="CC", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher <- runTest(myGOdata, algorithm='weight01', statistic="fisher")
    allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
    allRes$test <- "typhina_core_down_DE_CC"
    write.table(allRes, "all_Fisher_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "	", append = TRUE)
    showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo ='all')
    printGraph(myGOdata, resultFisher, firstSigNodes = 5, fn.prefix = "typhina_core_down_CC", useInfo = "all", pdfSW = TRUE)
    keep_fish(resultFisher)
    resultKS <- runTest(myGOdata, algorithm = "weight01", statistic = "ks") 
    KSRes <- GenTable(myGOdata, classicKS = resultKS, orderBy = "resultKS", ranksOf = "classicFisher", topNodes = 10)
    KSRes$test <- "typhina_core_down_DE_CC"
    write.table(KSRes, "all_KS_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "	", append = TRUE)
    keep_ks(resultKS)

    myGOdata <- new("topGOdata", description="typhina_core_down_DE_MF", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher <- runTest(myGOdata, algorithm='weight01', statistic="fisher")
    allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
    allRes$test <- "typhina_core_down_DE_MF"
    write.table(allRes, "all_Fisher_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "	", append = TRUE)
    showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo ='all')
    printGraph(myGOdata, resultFisher, firstSigNodes = 5, fn.prefix = "typhina_core_down_MF", useInfo = "all", pdfSW = TRUE)
    keep_fish(resultFisher)
    resultKS <- runTest(myGOdata, algorithm = "weight01", statistic = "ks") 
    KSRes <- GenTable(myGOdata, classicKS = resultKS, orderBy = "resultKS", ranksOf = "classicFisher", topNodes = 10)
    KSRes$test <- "typhina_core_down_DE_MF"
    write.table(KSRes, "all_KS_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "	", append = TRUE)
    keep_ks(resultKS)




#####     ELYMI
    
    
    # read in pannzer annotations and put into topGO datastructure geneID2GO
    pannzer <- read.delim("../../Pannzer/E.elymi_NfE728/E.elymi_NFE728_pannzer_GO.txt", header = TRUE, sep = "\t", colClasses="character")
    name <- "elymi"
    go_terms <- pannzer_to_topgo(pannzer)
    go_terms$qpid <- gsub("-T1", "", go_terms$qpid )
    write.table(go_terms, paste(name, "gene2GO.txt", sep = "_"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    geneID2GO <- readMappings(file = paste(name, "gene2GO.txt", sep = "_"))
    
    # subset DE information
    DE <- DE_all[, c( "elymi_gene_id", "core_up", "core_down", "top_4FC_up", "top_4FC_down", "top_2FC_up", "top_2FC_down")]
    
    # set the gene universe (all genes considered) and the genes of interest
    geneUniverse <- go_terms$qpid

# UPREGULATED
    
    genesOfInterest <- as.character(DE[(DE$core_up == 1), c("elymi_gene_id")])
    #ely_up_goi <- genesOfInterest
    # locations of genes of interest in the geneUniverse set in a named vector
    geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
    names(geneList) <- geneUniverse
    check_topgo_input(geneUniverse, geneList, geneID2GO)

    # create GOdata object
    myGOdata <- new("topGOdata", description="elymi_core_up_DE_BP", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher <- runTest(myGOdata, algorithm='weight01', statistic="fisher")
    allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
    allRes$test <- "elymi_core_up_DE_BP"
    write.table(allRes, "all_Fisher_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "	", append = TRUE)
    showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo ='all')
    printGraph(myGOdata, resultFisher, firstSigNodes = 5, fn.prefix = "elymi_core_up_BP", useInfo = "all", pdfSW = TRUE)
    keep_fish(resultFisher)
    resultKS <- runTest(myGOdata, algorithm = "weight01", statistic = "ks") 
    KSRes <- GenTable(myGOdata, classicKS = resultKS, orderBy = "resultKS", ranksOf = "classicFisher", topNodes = 10)
    KSRes$test <-  "elymi_core_up_DE_BP"
    write.table(KSRes, "all_KS_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "	", append = TRUE)
    keep_ks(resultKS)

    myGOdata <- new("topGOdata", description="elymi_core_up_DE_CC", ontology="CC", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher <- runTest(myGOdata, algorithm='weight01', statistic="fisher")
    allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
    allRes$test <- "elymi_core_up_DE_CC"
    write.table(allRes, "all_Fisher_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "	", append = TRUE)
    showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo ='all')
    printGraph(myGOdata, resultFisher, firstSigNodes = 5, fn.prefix = "elymi_core_up_CC", useInfo = "all", pdfSW = TRUE)
    keep_fish(resultFisher)
    resultKS <- runTest(myGOdata, algorithm = "weight01", statistic = "ks") 
    KSRes <- GenTable(myGOdata, classicKS = resultKS, orderBy = "resultKS", ranksOf = "classicFisher", topNodes = 10)
    KSRes$test <- "elymi_core_up_DE_CC"
    write.table(KSRes, "all_KS_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "	", append = TRUE)
    keep_ks(resultKS)

    myGOdata <- new("topGOdata", description="elymi_core_up_DE_MF", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher <- runTest(myGOdata, algorithm='weight01', statistic="fisher")
    allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
    allRes$test <- "elymi_core_up_DE_MF"
    write.table(allRes, "all_Fisher_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "	", append = TRUE)
    showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo ='all')
    printGraph(myGOdata, resultFisher, firstSigNodes = 5, fn.prefix = "elymi_core_up_MF", useInfo = "all", pdfSW = TRUE)
    keep_fish(resultFisher)
    resultKS <- runTest(myGOdata, algorithm = "weight01", statistic = "ks") 
    KSRes <- GenTable(myGOdata, classicKS = resultKS, orderBy = "resultKS", ranksOf = "classicFisher", topNodes = 10)
    KSRes$test <- "elymi_core_up_DE_MF"
    write.table(KSRes, "all_KS_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "	", append = TRUE)
    keep_ks(resultKS)

# DOWNREGULATED
    
    genesOfInterest <- as.character(DE[(DE$core_down == 1), c("elymi_gene_id")])
    #ely_down_goi <- genesOfInterest    
    # locations of genes of interest in the geneUniverse set in a named vector
    geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
    names(geneList) <- geneUniverse
    check_topgo_input(geneUniverse, geneList, geneID2GO)
    
    # create GOdata object
    myGOdata <- new("topGOdata", description="elymi_core_down_DE_BP", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher <- runTest(myGOdata, algorithm='weight01', statistic="fisher")
    allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
    allRes$test <- "elymi_core_down_DE_BP"
    write.table(allRes, "all_Fisher_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "	", append = TRUE)
    showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo ='all')
    printGraph(myGOdata, resultFisher, firstSigNodes = 5, fn.prefix = "elymi_core_down_BP", useInfo = "all", pdfSW = TRUE)
    keep_fish(resultFisher)
    resultKS <- runTest(myGOdata, algorithm = "weight01", statistic = "ks") 
    KSRes <- GenTable(myGOdata, classicKS = resultKS, orderBy = "resultKS", ranksOf = "classicFisher", topNodes = 10)
    KSRes$test <- "elymi_core_down_DE_BP"
    write.table(KSRes, "all_KS_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "	", append = TRUE)
    keep_ks(resultKS)

    myGOdata <- new("topGOdata", description="elymi_core_down_DE_CC", ontology="CC", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher <- runTest(myGOdata, algorithm='weight01', statistic="fisher")
    allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
    allRes$test <- "elymi_core_down_DE_CC"
    write.table(allRes, "all_Fisher_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "	", append = TRUE)
    showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo ='all')
    printGraph(myGOdata, resultFisher, firstSigNodes = 5, fn.prefix = "elymi_core_down_CC", useInfo = "all", pdfSW = TRUE)
    keep_fish(resultFisher)
    resultKS <- runTest(myGOdata, algorithm = "weight01", statistic = "ks") 
    KSRes <- GenTable(myGOdata, classicKS = resultKS, orderBy = "resultKS", ranksOf = "classicFisher", topNodes = 10)
    KSRes$test <- "elymi_core_down_DE_CC"
    write.table(KSRes, "all_KS_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "	", append = TRUE)
    keep_ks(resultKS)

    myGOdata <- new("topGOdata", description="elymi_core_down_DE_MF", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher <- runTest(myGOdata, algorithm='weight01', statistic="fisher")
    allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
    allRes$test <- "elymi_core_down_DE_MF"
    write.table(allRes, "all_Fisher_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "	", append = TRUE)
    showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo ='all')
    printGraph(myGOdata, resultFisher, firstSigNodes = 5, fn.prefix = "elymi_core_down_MF", useInfo = "all", pdfSW = TRUE)
    keep_fish(resultFisher)
    resultKS <- runTest(myGOdata, algorithm = "weight01", statistic = "ks") 
    KSRes <- GenTable(myGOdata, classicKS = resultKS, orderBy = "resultKS", ranksOf = "classicFisher", topNodes = 10)
    KSRes$test <- "elymi_core_down_DE_MF"
    write.table(KSRes, "all_KS_results_tests_top_ten.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "	", append = TRUE)
    keep_ks(resultKS)






