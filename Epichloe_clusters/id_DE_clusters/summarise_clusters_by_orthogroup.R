# Kate Lee 2019
# summarise the groups of clusters and arrange by orthogroup

library(tidyverse)

setwd("/media/kate/Massey_linux_onl/projects/STROMATA/results/Epichloe_clusters/id_DE_clusters/")

# read in ortholog flat file
orth <- read.delim("../../../../data/gene_model_sets/proteinortho/flat_ortho_file.txt", sep = " ", header = TRUE)
orth$spp <- NULL
orth$gene_id <- NULL
orth_list <- data.frame(unique(orth$ortho_group))
colnames(orth_list) <- c("ortho_group")

# read in cluster groups for each interaction
str_ps <- read.delim("STR_PS/cluster_groups.txt", header = TRUE, sep = "\t")
str_inf <- read.delim("STR_INF/cluster_groups.txt", header = TRUE, sep = "\t")
inf_ps <- read.delim("INF_PS/cluster_groups.txt", header = TRUE, sep = "\t")

# divide files by species
str_ps_el <- str_ps[ str_ps$species == "elymi", ]
str_ps_fs <- str_ps[ str_ps$species == "festucae", ]
str_ps_ty <- str_ps[ str_ps$species == "typhina", ]

str_inf_el <- str_inf[ str_inf$species == "elymi", ]
str_inf_fs <- str_inf[ str_inf$species == "festucae", ]
str_inf_ty <- str_inf[ str_inf$species == "typhina", ]

inf_ps_el <- inf_ps[ inf_ps$species == "elymi", ]
inf_ps_fs <- inf_ps[ inf_ps$species == "festucae", ]
inf_ps_ty <- inf_ps[ inf_ps$species == "typhina", ]


# annotate column headers with species identifier
colnames(str_ps_el) <- paste("el", colnames(str_ps_el), sep = "_")
colnames(str_ps_fs) <- paste("fs", colnames(str_ps_fs), sep = "_")
colnames(str_ps_ty) <- paste("ty", colnames(str_ps_ty), sep = "_")

colnames(str_inf_el) <- paste("el", colnames(str_inf_el), sep = "_")
colnames(str_inf_fs) <- paste("fs", colnames(str_inf_fs), sep = "_")
colnames(str_inf_ty) <- paste("ty", colnames(str_inf_ty), sep = "_")

colnames(inf_ps_el) <- paste("el", colnames(inf_ps_el), sep = "_")
colnames(inf_ps_fs) <- paste("fs", colnames(inf_ps_fs), sep = "_")
colnames(inf_ps_ty) <- paste("ty", colnames(inf_ps_ty), sep = "_")


# merge with ortholog list
str_ps_el <- merge(orth_list, str_ps_el, by.x = "ortho_group", by.y = "el_ortholog", all = TRUE)
str_ps_fs <- merge(orth_list, str_ps_fs, by.x = "ortho_group", by.y = "fs_ortholog", all = TRUE)
str_ps_ty <- merge(orth_list, str_ps_ty, by.x = "ortho_group", by.y = "ty_ortholog", all = TRUE)
str_ps_groups <- merge(str_ps_el, str_ps_fs, by = "ortho_group")
str_ps_groups <- merge(str_ps_groups, str_ps_ty, by = "ortho_group")

str_inf_el <- merge(orth_list, str_inf_el, by.x = "ortho_group", by.y = "el_ortholog", all = TRUE)
str_inf_fs <- merge(orth_list, str_inf_fs, by.x = "ortho_group", by.y = "fs_ortholog", all = TRUE)
str_inf_ty <- merge(orth_list, str_inf_ty, by.x = "ortho_group", by.y = "ty_ortholog", all = TRUE)
str_inf_groups <- merge(str_inf_el, str_inf_fs, by = "ortho_group")
str_inf_groups <- merge(str_inf_groups, str_inf_ty, by = "ortho_group")

inf_ps_el <- merge(orth_list, inf_ps_el, by.x = "ortho_group", by.y = "el_ortholog", all = TRUE)
inf_ps_fs <- merge(orth_list, inf_ps_fs, by.x = "ortho_group", by.y = "fs_ortholog", all = TRUE)
inf_ps_ty <- merge(orth_list, inf_ps_ty, by.x = "ortho_group", by.y = "ty_ortholog", all = TRUE)
inf_ps_groups <- merge(inf_ps_el, inf_ps_fs, by = "ortho_group")
inf_ps_groups <- merge(inf_ps_groups, inf_ps_ty, by = "ortho_group")


# add group and number of species columns
str_ps_groups$cluster_group <- ifelse(!is.na(str_ps_groups$el_group), str_ps_groups$el_group,
                              ifelse(!is.na(str_ps_groups$fs_group), str_ps_groups$fs_group,
                              ifelse(!is.na(str_ps_groups$ty_group), str_ps_groups$ty_group, "NA")))
str_ps_groups$cluster_num_spp <- ifelse(!is.na(str_ps_groups$el_num_spp), str_ps_groups$el_num_spp,
                              ifelse(!is.na(str_ps_groups$fs_num_spp), str_ps_groups$fs_num_spp,
                                     ifelse(!is.na(str_ps_groups$ty_num_spp), str_ps_groups$ty_num_spp, "NA")))

str_inf_groups$cluster_group <- ifelse(!is.na(str_inf_groups$el_group), str_inf_groups$el_group,
                                      ifelse(!is.na(str_inf_groups$fs_group), str_inf_groups$fs_group,
                                             ifelse(!is.na(str_inf_groups$ty_group), str_inf_groups$ty_group, "NA")))
str_inf_groups$cluster_num_spp <- ifelse(!is.na(str_inf_groups$el_num_spp), str_inf_groups$el_num_spp,
                                        ifelse(!is.na(str_inf_groups$fs_num_spp), str_inf_groups$fs_num_spp,
                                               ifelse(!is.na(str_inf_groups$ty_num_spp), str_inf_groups$ty_num_spp, "NA")))

inf_ps_groups$cluster_group <- ifelse(!is.na(inf_ps_groups$el_group), inf_ps_groups$el_group,
                                      ifelse(!is.na(inf_ps_groups$fs_group), inf_ps_groups$fs_group,
                                             ifelse(!is.na(inf_ps_groups$ty_group), inf_ps_groups$ty_group, "NA")))
inf_ps_groups$cluster_num_spp <- ifelse(!is.na(inf_ps_groups$el_num_spp), inf_ps_groups$el_num_spp,
                                        ifelse(!is.na(inf_ps_groups$fs_num_spp), inf_ps_groups$fs_num_spp,
                                               ifelse(!is.na(inf_ps_groups$ty_num_spp), inf_ps_groups$ty_num_spp, "NA")))

keep <- c("ortho_group", "el_cluster", "el_cluster_size", "el_gene_id", "el_start", "el_end", "el_strand", "fs_cluster",
"fs_cluster_size", "fs_gene_id", "fs_start", "fs_end", "fs_strand", "ty_cluster", "ty_cluster_size",
"ty_gene_id", "ty_start", "ty_end", "ty_strand", "cluster_group", "cluster_num_spp")

str_ps_groups <- str_ps_groups[, keep]

keep <- c("ortho_group", "el_cluster", "el_cluster_size", "el_gene_id", "el_start", "el_end", "el_strand", "fs_cluster",
          "fs_cluster_size", "fs_gene_id", "fs_start", "fs_end", "fs_strand", "cluster_group", "cluster_num_spp")

str_inf_groups <- str_inf_groups[, keep]
inf_ps_groups <- str_inf_groups[, keep]

write.table(str_ps_groups, "STR_PS/str_ps_groups_wide_file_for_ann.txt", quote = FALSE, row.names = FALSE)
write.table(str_inf_groups, "STR_INF/str_inf_groups_wide_file_for_ann.txt", quote = FALSE, row.names = FALSE)
write.table(inf_ps_groups, "INF_PS/inf_ps_groups_wide_file_for_ann.txt", quote = FALSE, row.names = FALSE)
