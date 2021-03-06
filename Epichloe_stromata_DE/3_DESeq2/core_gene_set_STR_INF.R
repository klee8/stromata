# "core_gene_set"

#install.packages("tidyverse")
library("tidyverse")

resdir <- ("/media/kate/Massey_linux_onl/projects/stromata_analysis/Epichloe_stromata_DE/3_DESeq2/")


### get_data from elymi and festucae (typhina does not have infloresence data)
elymi <- read.delim(paste(resdir, "E.elymi_NfE728/result_tables/E.elymi_NfE728_STR_INF.txt", sep = ""), header = TRUE)
festucae <- read.delim(paste(resdir, "E.festucae_E2368/result_tables/E.festucae_E2368_STR_INF.txt", sep = ""), header = TRUE)

# subset apeglm results
keep <- c("gene_id", "apeglm_log2FC", "apeglm_lfcSE", "apeglm_1_svalue", "apeglm_2_svalue", "mean_STR_TPM", "mean_INF_TPM", "up_2FC", "down_2FC")
elymi <- elymi[, keep]
festucae <- festucae[, keep]

# set colnames for each set of results
colnames(elymi) <- paste("elymi", colnames(elymi), sep = "_")
colnames(festucae) <- paste("festucae", colnames(festucae), sep = "_")

# rename_genes_by_ortholog_names
orthology <- read.delim("flat_ortho_file.txt", header = TRUE)
orthology$gene_id <- gsub("-T1", "", orthology$gene_id)

elymi <- merge(elymi, orthology[, c("ortho_group",  "gene_id")], by.x = "elymi_gene_id", by.y = "gene_id")
festucae <- merge(festucae, orthology[, c("ortho_group",  "gene_id")], by.x = "festucae_gene_id", by.y = "gene_id")

# merge result tables
core_set <- merge(elymi, festucae, by = "ortho_group", all = TRUE)

# add number of species the ortho_group is present in (i.e. 1/2/3 of elymi, festucae and typhina)
num_orth <- read.table("lookup_ortholog_groups.txt", header = TRUE)
num_orth <- num_orth %>% select(orthogroup, elymi, fescuae, typhina) %>% mutate(num_spp = rowSums(. != "*") - 1)
core_set <- merge(num_orth[, c("orthogroup", "num_spp")], core_set, by.x = "orthogroup", by.y = "ortho_group", all = TRUE)

# add m3 number so Chris' annotations can be looked up easily
m3 <- read.delim("lookup_ortholog_groups.txt", header = TRUE, sep = "\t")
m3 <- m3[, c("orthogroup", "m3")]
core_set <- merge(core_set, m3, by.x = "orthogroup", by.y = "orthogroup", all = FALSE)

rownames(core_set) <- core_set$group
colnames(core_set)


##### CORE: flag if up by 4fold in one data set and by up by at least 2 fold in the other
#targets_up <-  rownames(subset(core_set, ( ((elymi_apeglm_2_svalue <.005 & elymi_apeglm_log2FC >= 2 ) 
#                                          | (festucae_apeglm_2_svalue <.005 & festucae_apeglm_log2FC >= 2 ) )
#                                          & (elymi_apeglm_1_svalue <.005 & elymi_apeglm_log2FC >= 1 ) 
#                                          & (festucae_apeglm_1_svalue <.005 & festucae_apeglm_log2FC >= 1 ) )))

#core_set$core_up <- ifelse(rownames(core_set) %in% targets_up, 1, 0)

#targets_down <-  rownames(subset(core_set, ( ((elymi_apeglm_2_svalue <.005 & elymi_apeglm_log2FC <= -2 ) 
#                                          | (festucae_apeglm_2_svalue <.005 & festucae_apeglm_log2FC <= -2 ) )
#                                          & (elymi_apeglm_1_svalue <.005 & elymi_apeglm_log2FC <= -1 ) 
#                                          & (festucae_apeglm_1_svalue <.005 & festucae_apeglm_log2FC <= -1 ) )))

#core_set$core_down <- ifelse(rownames(core_set) %in% targets_down, 1, 0)


##### 4 fold change
#targets_topup <-  rownames(subset(core_set, ( (elymi_apeglm_2_svalue <.005 & elymi_apeglm_log2FC >= 2 ) 
#                                          & (festucae_apeglm_2_svalue <.005 & festucae_apeglm_log2FC >= 2 ) )))

#core_set$top_4FC_up <- ifelse(rownames(core_set) %in% targets_topup, 1, 0)

#targets_topdown <-  rownames(subset(core_set, ( (elymi_apeglm_2_svalue <.005 & elymi_apeglm_log2FC <= -2 ) 
#                                          & (festucae_apeglm_2_svalue <.005 & festucae_apeglm_log2FC <= -2 ) )))

#core_set$top_4FC_down <- ifelse(rownames(core_set) %in% targets_topdown, 1, 0)



##### 2 fold change
targets_topup <-  rownames(subset(core_set, ( (elymi_apeglm_1_svalue <.005 & elymi_apeglm_log2FC >= 1 ) 
                                          & (festucae_apeglm_1_svalue <.005 & festucae_apeglm_log2FC >= 1 ) )))

core_set$core_2FC_up <- ifelse(rownames(core_set) %in% targets_topup, 1, 0)

targets_topdown <-  rownames(subset(core_set, ( (elymi_apeglm_1_svalue <.005 & elymi_apeglm_log2FC <= -1 ) 
                                          & (festucae_apeglm_1_svalue <.005 & festucae_apeglm_log2FC <= -1 ) )))

core_set$core_2FC_down <- ifelse(rownames(core_set) %in% targets_topdown, 1, 0)


write.table(core_set, paste(resdir, "STR_INF/core_set_STR_INF.txt", sep = ""), quote = FALSE, col.names = TRUE, row.names = FALSE)




# Graph overlaps with UpSetR package
# UpSet_diagram_apeglm_normalised_results}
#install.packages("UpSetR")
library(UpSetR)

rownames(core_set) <- core_set$gene_id

listInput <- list(elymi_2FC_up = rownames(core_set[core_set$elymi_apeglm_log2FC >= 1 & core_set$elymi_apeglm_1_svalue < 0.005 & !is.na(core_set$elymi_apeglm_log2FC), ]),
                 # elymi_4Fold_up = rownames(core_set[core_set$elymi_apeglm_log2FC >= 2 & core_set$elymi_apeglm_2_svalue < 0.005 & !is.na(core_set$elymi_apeglm_log2FC), ]),
                  festucae_2FC_up = rownames(core_set[core_set$festucae_apeglm_log2FC >= 1 & core_set$festucae_apeglm_1_svalue < 0.005 & !is.na(core_set$festucae_apeglm_log2FC), ]),
                 # festucae_4Fold_up = rownames(core_set[core_set$festucae_apeglm_log2FC >= 2 & core_set$festucae_apeglm_2_svalue < 0.005 & !is.na(core_set$festucae_apeglm_log2FC), ]),

                  
                  elymi_2FC_down = rownames(core_set[core_set$elymi_apeglm_log2FC <= -1 & core_set$elymi_apeglm_1_svalue < 0.005 & !is.na(core_set$elymi_apeglm_log2FC), ]),
                 # elymi_4Fold_down = rownames(core_set[core_set$elymi_apeglm_log2FC <= -2 & core_set$elymi_apeglm_2_svalue < 0.005 & !is.na(core_set$elymi_apeglm_log2FC), ]),
                  festucae_2FC_down = rownames(core_set[core_set$festucae_apeglm_log2FC <= -1 & core_set$festucae_apeglm_1_svalue < 0.005 & !is.na(core_set$festucae_apeglm_log2FC), ]),
                 # festucae_4Fold_down = rownames(core_set[core_set$festucae_apeglm_log2FC <= -2 & core_set$festucae_apeglm_2_svalue < 0.005 & !is.na(core_set$festucae_apeglm_log2FC), ]),

                  
                 # core_up = rownames(core_set[core_set$core_up == 1 & !is.na(core_set$elymi_apeglm_log2FC) & !is.na(core_set$festucae_apeglm_log2FC), ]),
                 # core_down = rownames(core_set[core_set$core_down == 1 & !is.na(core_set$elymi_apeglm_log2FC) & !is.na(core_set$festucae_apeglm_log2FC), ]),
                  core_2FC_up = rownames(core_set[core_set$core_2FC_up == 1 & !is.na(core_set$elymi_apeglm_log2FC) & !is.na(core_set$festucae_apeglm_log2FC), ]),
                  core_2FC_down = rownames(core_set[core_set$core_2FC_down == 1 & !is.na(core_set$elymi_apeglm_log2FC) & !is.na(core_set$festucae_apeglm_log2FC), ])
                  
                  )

sets <- names(listInput)
colors <- c("blue", "blue", "blue", "purple", "purple", "purple")
metadata <- as.data.frame(cbind(sets, colors))

set_colors <- c(festucae_2FC_down = "darkcyan",  elymi_2FC_down = "darkcyan",  core_2FC_down = "darkblue", festucae_2FC_up = "purple", elymi_2FC_up = "purple",  core_2FC_up = "red") 

pdf(paste(resdir, "STR_INF/UpSet_core_set_STR_INF.pdf", sep = ""), onefile = FALSE)
upset(fromList(listInput), sets = c("core_2FC_down", "core_2FC_up"), keep.order = TRUE, order.by = "freq")
dev.off()

pdf(paste(resdir, "STR_INF/UpSet_all_FC_data_STR_INF.pdf", sep = ""), onefile = FALSE)
upset(fromList(listInput), sets = c("festucae_2FC_down",  "elymi_2FC_down", "core_2FC_down", "festucae_2FC_up", "elymi_2FC_up", "core_2FC_up"), 
			  keep.order = TRUE, order.by = "freq", nintersects = 50, set.metadata = list(data = metadata, 
                          plots = list(list(type = "matrix_rows", column = "sets", colors = set_colors, alpha = 0.5))))
dev.off()

pdf(paste(resdir, "STR_INF/UpSet_all_FC_data_up_STR_INF.pdf", sep = ""), onefile = FALSE)
upset(fromList(listInput), sets = c("festucae_2FC_up", "elymi_2FC_up", "core_2FC_up"), keep.order = TRUE, order.by = "freq")
dev.off()

pdf(paste(resdir, "STR_INF/UpSet_all_FC_data_down_STR_INF.pdf", sep = ""), onefile = FALSE)
upset(fromList(listInput), sets = c("festucae_2FC_down",  "elymi_2FC_down", "core_2FC_down"), keep.order = TRUE, order.by = "freq")
dev.off()



# add_annoatations}

#### GET ANNOTATIONS FROM annotations.txt file from David's funannotate pipeline file
elymi_ann <- read.delim(paste(resdir, "../0_raw_data/genome_assembly/Eel_728/Epichloe_elymi.annotations.txt", sep = ""), header = TRUE, comment.char = "#", blank.lines.skip = TRUE, fill = TRUE, sep = "\t")
festucae_ann <- read.delim(paste(resdir, "../0_raw_data/genome_assembly/Efe_2368/Epichloe_festucae.annotations.txt", sep = ""), header = TRUE, comment.char = "#", blank.lines.skip = TRUE, fill = TRUE, sep = "\t")

elymi_ann$protein_length <- nchar(as.character(elymi_ann$Translation))
festucae_ann$protein_length <- nchar(as.character(festucae_ann$Translation))

# tag each spp annotation
colnames(elymi_ann) <- paste("elymi", colnames(elymi_ann), sep = "_")
colnames(festucae_ann) <- paste("festucae", colnames(festucae_ann), sep = "_")

# merge DESeq output 
core_set <- merge(core_set, festucae_ann, by.x = "festucae_gene_id", by.y = "festucae_GeneID", all.x = TRUE)
core_set <- merge(core_set, elymi_ann, by.x = "elymi_gene_id", by.y = "elymi_GeneID", all.x = TRUE)

# (note, removed the apeglm_2_svalues from this set)
core_set <- core_set[, c("orthogroup", "num_spp", "m3", "elymi_gene_id", "festucae_gene_id", 
                         "elymi_apeglm_log2FC", "elymi_apeglm_lfcSE", "elymi_apeglm_1_svalue", "elymi_mean_STR_TPM", "elymi_mean_INF_TPM", 
                         "festucae_apeglm_log2FC", "festucae_apeglm_lfcSE", "festucae_apeglm_1_svalue", "festucae_mean_STR_TPM", "festucae_mean_INF_TPM",
                         "core_2FC_up", "core_2FC_down", "elymi_up_2FC", "elymi_down_2FC", "festucae_up_2FC", "festucae_down_2FC",
                         "elymi_Feature", "elymi_Contig", "elymi_Start", "elymi_Stop", "elymi_Strand", "elymi_Name", "elymi_Product", 
                         "elymi_BUSCO", "elymi_PFAM", "elymi_InterPro", "elymi_EggNog", "elymi_COG", "elymi_GO.Terms", "elymi_Secreted",  
                         "elymi_Membrane", "elymi_Protease", "elymi_CAZyme", "elymi_Notes", "elymi_Translation", "elymi_protein_length",  
                         "festucae_Feature", "festucae_Contig", "festucae_Start", "festucae_Stop", "festucae_Strand", "festucae_Name", "festucae_Product",
                         "festucae_BUSCO", "festucae_PFAM", "festucae_InterPro", "festucae_EggNog", "festucae_COG", "festucae_GO.Terms", "festucae_Secreted",
                         "festucae_Membrane", "festucae_Protease", "festucae_CAZyme", "festucae_Notes", "festucae_Translation", "festucae_protein_length" ) ]


write.table(core_set, paste(resdir, "STR_INF/core_gene_set_STR_INF_ann.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")
  

# pannzer_ann
# (note, removed the apeglm_2_svalues from this set)
core_set_pann <- core_set[, c("orthogroup", "num_spp", "m3", "elymi_gene_id", "festucae_gene_id", 
                              "elymi_apeglm_log2FC", "elymi_apeglm_lfcSE", "elymi_apeglm_1_svalue", "elymi_mean_STR_TPM", "elymi_mean_INF_TPM",
                              "festucae_apeglm_log2FC", "festucae_apeglm_lfcSE",  "festucae_apeglm_1_svalue", "festucae_mean_STR_TPM", "festucae_mean_INF_TPM",
                              "core_2FC_up", "core_2FC_down", "elymi_up_2FC", "elymi_down_2FC", "festucae_up_2FC", "festucae_down_2FC")]      

# add SignalP annotations
signalP <- read.delim(paste(resdir, "../data/SignalP/E.elymi_NfE728_signalP.txt", sep = ""), header = FALSE, blank.lines.skip = TRUE, fill = TRUE, sep = "\t", comment.char = "#")
colnames(signalP) <- c("gene_id", "Pr_signalP", "sigP_pos", "sigP_clev", "Pr_sigP_clev")
signalP$gene_id <- gsub("-T1", "", signalP$gene_id)
colnames(signalP) <- paste("el", colnames(signalP), sep = "_")
core_set_pann <- merge(core_set_pann, signalP, by.x = "elymi_gene_id", by.y = "el_gene_id", all = TRUE)

signalP <- read.delim(paste(resdir, "../data/SignalP/E.festucae_E2368_signalP.txt", sep = ""), header = FALSE, blank.lines.skip = TRUE, fill = TRUE, sep = "\t", comment.char = "#")
colnames(signalP) <- c("gene_id", "Pr_signalP", "sigP_pos", "sigP_clev", "Pr_sigP_clev")
signalP$gene_id <- gsub("-T1", "", signalP$gene_id)
colnames(signalP) <- paste("fe", colnames(signalP), sep = "_")
core_set_pann <- merge(core_set_pann, signalP, by.x = "festucae_gene_id", by.y = "fe_gene_id", all = TRUE)

# add EffectorP annotations
effectorP <- read.delim(paste(resdir, "../data/EffectorP/E.elymi_NfE728_effectorP.txt", sep = ""), header = TRUE, blank.lines.skip = TRUE, fill = TRUE, sep = "\t", comment.char = "#")
colnames(effectorP) <- c("gene_id", "effectorP", "efP_Pr")
effectorP$gene_id <- gsub("-T1.*", "", effectorP$gene_id)
colnames(effectorP) <- paste("el", colnames(effectorP), sep = "_")
core_set_pann <- merge(core_set_pann, effectorP,  by.x = "elymi_gene_id", by.y = "el_gene_id", all = TRUE)

effectorP <- read.delim(paste(resdir, "../data/EffectorP/E.festucae_E2368_effectorP.txt", sep = ""), header = TRUE, blank.lines.skip = TRUE, fill = TRUE, sep = "\t", comment.char = "#")
colnames(effectorP) <- c("gene_id", "effectorP", "efP_Pr")
effectorP$gene_id <- gsub("-T1.*", "", effectorP$gene_id)
colnames(effectorP) <- paste("fe", colnames(effectorP), sep = "_")
core_set_pann <- merge(core_set_pann, effectorP,  by.x = "festucae_gene_id", by.y = "fe_gene_id", all = TRUE)

# add PANNZER DE annoations
pannzer <- read.delim(paste(resdir, "../data/Pannzer/E.elymi_NfE728/E.elymi_NFE728_pannzer_DE.txt", sep = ""), header = TRUE, comment.char = "#", blank.lines.skip = TRUE, fill = TRUE, sep = "\t")
pannzer$gene_id <- gsub("-T1", "", pannzer$qpid)
pan_keep <- pannzer[, c("gene_id", "desc", "genename")]
colnames(pan_keep) <- c("gene_id", "pannzer_desc", "pannzer_genename")
colnames(pan_keep) <- paste("el", colnames(pan_keep), sep = "_")
core_set_pann <- merge(core_set_pann, pan_keep, by.x = "elymi_gene_id", by.y = "el_gene_id", all = TRUE)

pannzer <- read.delim(paste(resdir, "../data/Pannzer/E.festucae_E2368/E.festucae_E2368_pannzer_DE.txt", sep = ""), header = TRUE, comment.char = "#", blank.lines.skip = TRUE, fill = TRUE, sep = "\t")
pannzer$gene_id <- gsub("-T1", "", pannzer$qpid)
pan_keep <- pannzer[, c("gene_id", "desc", "genename")]
colnames(pan_keep) <- c("gene_id", "pannzer_desc", "pannzer_genename")
colnames(pan_keep) <- paste("fe", colnames(pan_keep), sep = "_")
core_set_pann <- merge(core_set_pann, pan_keep, by.x = "festucae_gene_id", by.y = "fe_gene_id", all = TRUE)

core_set_pann <- core_set_pann[, c( "orthogroup", "num_spp","m3", "elymi_gene_id", "festucae_gene_id", 
                                    "elymi_apeglm_log2FC", "elymi_apeglm_lfcSE", "elymi_apeglm_1_svalue", "elymi_mean_STR_TPM", "elymi_mean_INF_TPM", 
                                    "festucae_apeglm_log2FC", "festucae_apeglm_lfcSE", "festucae_apeglm_1_svalue", "festucae_mean_STR_TPM", "festucae_mean_INF_TPM", 
                                    "core_2FC_up", "core_2FC_down", "elymi_up_2FC", "elymi_down_2FC", "festucae_up_2FC", "festucae_down_2FC",
                                    "el_Pr_signalP", "el_sigP_pos", "el_sigP_clev", "el_Pr_sigP_clev", 
                                    "fe_Pr_signalP", "fe_sigP_pos", "fe_sigP_clev", "fe_Pr_sigP_clev", 
                                    "el_effectorP", "el_efP_Pr", 
                                    "fe_effectorP", "fe_efP_Pr", 
                                    "el_pannzer_desc","el_pannzer_genename", 
                                    "fe_pannzer_desc", "fe_pannzer_genename")]


write.table(core_set_pann, paste(resdir, "STR_INF/core_gene_set_STR_INF_pannzer.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")


