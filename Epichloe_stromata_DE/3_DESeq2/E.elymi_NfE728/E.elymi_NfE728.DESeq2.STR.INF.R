# "E.elymi_NfE728_STR_vs_INF"
# Following example from https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html


# package installation via BiocManager if required
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DESeq2", version = "3.9")
#BiocManager::install("tximport", version = "3.9")
#BiocManager::install("readr", version = "3.9")
#BiocManager::install("tidyverse", version = "3.9")
#BiocManager::install("hexbin", version = "3.9")
#BiocManager::install("pheatmap", version = "3.9")
#BiocManager::install("RColorBrewer", version = "3.9")
#BiocManager::install("apeglm", version = "3.9")

library("DESeq2")
library("tximport")
library("readr")
library("tidyverse")
library("hexbin")
library("pheatmap")
library("RColorBrewer")
library("apeglm")
##BiocManager::install("dplyr", version = "3.9")
##BiocManager::install("ggplot2", version = "3.9")
#library("dplyr")
#library("ggplot2")



### Use Tximport to import data from Salmon
# POP = species used
# center = sequencing center
# assay = sequencing
# run_number = sequencing run (sent some samples back for re-sequencing)
# lane = sequencing lane number
# sample = your name for the sample
# experiment = name for the experiment
# run = Illumina base-name for the fastq file, or something similar which indicates sample, lane etc (must be unique for each row)
# condition = treatment (factor levels)


### get_data
# import sample list and assign column and row names
samples <- read.table("../../config.txt", header = FALSE)
colnames(samples) <- c("pop", "center", "assay", "run_number", "lane", "sample", "experiment", "run", "condition")
rownames(samples) <- samples$run

# subset only one species
samples <- samples[samples$pop == "E.elymi_NfE728",]

# subset conditions of interest
samples <- samples[samples$condition == "INF" | samples$condition == "STR",]

# get_salmon_quant_counts_filenames}
dir <- "../../2_salmon"
files <- file.path(dir, "salmon_quant", samples$pop, samples$run, "quant.sf")
names(files) <- samples$run

# get_transcript_and_gene_names}
tx2gene <- read.table(file.path(dir, "transcriptome/E.elymi_NfE728/transcriptome.headers.txt"))
colnames(tx2gene) <- c("TXNAME", "GENEID")

# import salmon count data to R
txi <- tximport(files, type="salmon", tx2gene=tx2gene, abundanceCol = "TPM")



### calculate Salmon TPM abundance
# general solution, get list of samples for each condition
STR_list <- as.vector(samples[samples$condition == "STR", "run"])
INF_list <- as.vector(samples[samples$condition == "INF", "run"])

# divide by the number of unique samples (some samples are duplicated over multiple lanes)
STR_samplenum <- length(as.vector(unique(samples[samples$condition == "STR", "sample"])))
INF_samplenum <- length(as.vector(unique(samples[samples$condition == "INF", "sample"])))

salmon_TPM <- as.data.frame (txi$abundance) %>% 
  mutate( mean_STR_TPM = rowSums(select(., STR_list)/STR_samplenum),
          mean_INF_TPM = rowSums(select(., INF_list)/INF_samplenum),
          gene_id = rownames(txi$abundance)) %>%
  select(gene_id, mean_STR_TPM, mean_INF_TPM) 



### create initial DESeq data structure (dds)
dds <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)

# collapse_technical_replicates}
dds <- collapseReplicates(dds, groupby = samples$sample, run = samples$run, renameCols = TRUE)



### Initial data exploratory analysis
# Use regularized-logarithm transformation or rlog (Love, Huber, and Anders 2014) to 
# normalise raw data prior to createing PCA and checking sample distance. 
# An alternative transformation is vst, which is faster)
# specified blind = FALSE, which means that differences between cell lines and treatment
# (the variables in the design) will not contribute to the expected variance-mean trend 
# of the experiment.
# Don't use this transform in the DESeq analysis this because for differential testing 
# the DESeq function applied to raw counts is recommended

# rlog_trans
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)

# plot_rlog_trans_effects
dds <- estimateSizeFactors(dds)

df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
         mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

pdf("plots/STR_vs_INF_rlog2_normalised_log2FC.pdf", onefile = FALSE)
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  

dev.off()


# calculate sample distances
sampleDists <- dist(t(assay(rld)))

# heatmap of sample distances
sampleDistMatrix <- as.matrix( sampleDists )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf("plots/STR_vs_INF_heat_map_sample_dist_rld-normalised_data.pdf", onefile = FALSE)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

dev.off()


# Plot PCA of rlog transformed data with plotPCA
# removed coord-fixed from plotPCA function ggplot
plotPCA <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE)
{
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
        length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
        stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup,
        drop = FALSE])
    group <- if (length(intgroup) > 1) {
        factor(apply(intgroup.df, 1, paste, collapse = ":"))
    }
    else {
        colData(object)[[intgroup]]
    }
    d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group,
        intgroup.df, name = colnames(object))
    if (returnData) {
        attr(d, "percentVar") <- percentVar[1:2]
        return(d)
    }
    ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) +
      geom_point(size = 3) +
      xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
      ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) #+
      #coord_fixed()
}

pdf("plots/STR_v_INF_PCA_plot_rld-normalised_data.pdf", onefile = FALSE)
plotPCA(rld, intgroup="condition")
while (!is.null(dev.list()))  dev.off()



###  Run DESeq algorithm
# If you never tell the DESeq2 functions which level you want to compare
# against (e.g. which level represents the control group), the comparisons
# will be based on the alphabetical order of the levels.
# Here just set the reference condition to wild type
dds$condition <- relevel(dds$condition, ref = "INF")


# DESeq consists of three steps
#   estimating size factors (gene level differences, these will be overridden by avgTxLength from Tximport)
#   estimating dispersions (differences between biological replicates)
#   fit data to negative binomial and conduct Wald test

# run DESeq
dds <- DESeq(dds)
res <- results(dds)

# diff_gene_exp
summary(res) #summary of results

# summary sorted by p vlaue adjusted
res <- res[order(res$padj),]

# plot genes with top 6 size diffs
pdf("plots/STR_vs_INF_top_DE_genes.pdf" )
top_counts <- rownames(res)
par(mfrow=c(2,3))
plotCounts(dds, gene=top_counts[1], intgroup="condition")
plotCounts(dds, gene=top_counts[2], intgroup="condition")
plotCounts(dds, gene=top_counts[3], intgroup="condition")
plotCounts(dds, gene=top_counts[4], intgroup="condition")
plotCounts(dds, gene=top_counts[5], intgroup="condition")
plotCounts(dds, gene=top_counts[6], intgroup="condition")
dev.off()

# plot_dispersions
pdf("plots/STR_vs_INF_plot_dispersions.pdf", onefile = FALSE)
plotDispEsts(dds)
dev.off()


### Prior to plotting DE's shrink the data with apeglm
# shrink results using thresholds of log2FC of 1 or 2
resApeT_thresh1  <- lfcShrink(dds, coef="condition_STR_vs_INF", type="apeglm", lfcThreshold = 1)
resApeT_thresh2  <- lfcShrink(dds, coef="condition_STR_vs_INF", type="apeglm", lfcThreshold = 2)

pdf("plots/STR_vs_INF_apeglm_shrinkage_MA_plot.pdf", onefile = FALSE)
drawLines <- function() abline(h=c(-2,2),col="dodgerblue",lwd=2)
plotMA(resApeT_thresh2, ylim=c(-8,8), cex=.8); drawLines()
dev.off()

colnames(resApeT_thresh1)  <- c("baseMean", "apeglm_1_log2FC", "apeglm_1_lfcSE", "apeglm_1_svalue" )
colnames(resApeT_thresh2)  <- c("baseMean", "apeglm_2_log2FC", "apeglm_2_lfcSE", "apeglm_2_svalue" )


# volcano plot
pdf("plots/STR_vs_INF_volcano_plot.pdf", onefile = FALSE)
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()



### summarise results
# set results as a data frame
res$gene_id <- rownames(res)
all_results <- as.data.frame(res)

# merge data after apeglm shrinkage 
resApeT_thresh1$gene_id <- rownames(resApeT_thresh1)
resApeT_thresh2$gene_id <- rownames(resApeT_thresh2)
all_apeglmT <- merge(as.data.frame(resApeT_thresh1), as.data.frame(resApeT_thresh2), by = "gene_id")
all_results <- merge(all_results, all_apeglmT, by = "gene_id")

# merge TPM from salmon
all_results <- merge(all_results, salmon_TPM, by = "gene_id")
all_results$apeglm_log2FC <- all_results$apeglm_1_log2FC
all_results$apeglm_lfcSE <- all_results$apeglm_1_lfcSE
goi <- all_results %>% select(gene_id, apeglm_log2FC, apeglm_lfcSE, apeglm_1_svalue, apeglm_2_svalue, mean_INF_TPM, mean_STR_TPM)

write.table(goi, "result_tables/E.elymi_NfE728_STR_INF.txt", quote = FALSE, row.names = FALSE, sep = "\t")

2FConly <- all_results %>% select(gene_id, apeglm_log2FC, apeglm_lfcSE, apeglm_1_svalue, mean_INF_TPM, mean_STR_TPM)

# identify significant DE











