# Kate Lee 2019

# read in ortholog flat file
orth <- read.delim("flat_ortho_file.txt", sep = " ", header = TRUE) 
orth$spp <- NULL

# get SMURF secondary metabolite clusters
ely <- read.delim("E.elymi_NfE728_results/SMURF_ely_rfmt.txt", sep = "\t", header = TRUE)
fes <- read.delim("E.festucae_E2368_results/SMURF_fes_rfmt.txt", sep = "\t", header = TRUE)
typ <- read.delim("E.typhina_E8_results/SMURF_typ_rfmt.txt", sep = "\t", header = TRUE)

# get backbone proteins (that SMURF clusters are based on)
ely_bk <- read.delim("E.elymi_NfE728_results/Backbone-genes.txt", sep = "\t", header = TRUE)
fes_bk <- read.delim("E.festucae_E2368_results/Backbone-genes.txt", sep = "\t", header = TRUE)
typ_bk <- read.delim("E.typhina_E8_results/Backbone-genes.txt", sep = "\t", header = TRUE)

# reduce backbone protein table to just gene id and SMURF predicted function
ely_bk <- ely_bk[, c("Backbone_gene_id", "SMURF_backbone_gene_prediction")]
fes_bk <- fes_bk[, c("Backbone_gene_id", "SMURF_backbone_gene_prediction")]
typ_bk <- typ_bk[, c("Backbone_gene_id", "SMURF_backbone_gene_prediction")]

# replace Annotated_gene_function in full output to just the SMURF backbone gene prediction
ely$Annotated_gene_function <- NULL
fes$Annotated_gene_function <- NULL
typ$Annotated_gene_function <- NULL

ely <- merge(ely, ely_bk, by = "Backbone_gene_id")
fes <- merge(fes, fes_bk, by = "Backbone_gene_id")
typ <- merge(typ, typ_bk, by = "Backbone_gene_id")

# add in spp name to colum headers
colnames(ely) <- paste("el", colnames(ely), sep = "_")
colnames(fes) <- paste("fes", colnames(fes), sep = "_")
colnames(typ) <- paste("typ", colnames(typ), sep = "_")

# add in ortholog numbers
ely <- merge(ely, orth, by.x = "el_Gene_id", by.y = "gene_id")
fes <- merge(fes, orth, by.x = "fes_Gene_id", by.y = "gene_id")
typ <- merge(typ, orth, by.x = "typ_Gene_id", by.y = "gene_id")

# remove gene id (will be elsewhere in spreadsheet)
ely$el_Gene_id <- NULL
fes$fes_Gene_id <- NULL
typ$typ_Gene_id <- NULL

# merge everything
orth$gene_id <- NULL
orth_list <- data.frame(unique(orth$ortho_group))
colnames(orth_list) <- c("ortho_group")

SMURF <- merge(orth_list, ely, by = "ortho_group", all = TRUE)
SMURF <- merge(SMURF, fes, by = "ortho_group", all = TRUE)
SMURF <- merge(SMURF, typ, by = "ortho_group", all = TRUE)

head(SMURF)
write.table(SMURF, "SMURF_flat_results_by_orthologs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

