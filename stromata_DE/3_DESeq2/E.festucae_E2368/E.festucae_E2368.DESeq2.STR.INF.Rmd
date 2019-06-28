---
title: "E.festucae_E2368_STR_vs_INF"
output: html_document
---
Following example from https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)

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
```


#Use Tximport to import data from Salmon
POP = species used
center = sequencing center
assay = sequencing
run_number = sequencing run (sent some samples back for re-sequencing)
lane = sequencing lane number
sample = your name for the sample
experiment = name for the experiment
run = Illumina base-name for the fastq file, or something similar which indicates sample, lane etc (must be unique for each row)
condition = treatment (factor levels)

```{r get_data}

# import sample list and assign column and row names
samples <- read.table("../../config.txt", header = FALSE)
colnames(samples) <- c("pop", "center", "assay", "run_number", "lane", "sample", "experiment", "run", "condition")
rownames(samples) <- samples$run

# subset only one species
samples <- samples[samples$pop == "E.festucae_E2368",]

# subset conditions of interest
samples <- samples[samples$condition == "INF" | samples$condition == "STR",]
```


```{r get_salmon_quant_counts_filenames}

dir <- "../../2_salmon"
files <- file.path(dir, "salmon_quant", samples$pop, samples$run, "quant.sf")
names(files) <- samples$run
```


```{r get_transcript_and_gene_names}

tx2gene <- read.table(file.path(dir, "transcriptome/E.festucae_E2368/transcriptome.headers.txt"))
colnames(tx2gene) <- c("TXNAME", "GENEID")
```


```{r import_data}

txi <- tximport(files, type="salmon", tx2gene=tx2gene, abundanceCol = "TPM")
```

```{r get_Salmon_TPM_abundance}

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
```


```{r create_initial_dds}

dds <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)
```

```{r collapse_technical_replicates}

dds <- collapseReplicates(dds, groupby = samples$sample, run = samples$run, renameCols = TRUE)
```


#Initial data exploratory analysis

Use regularized-logarithm transformation or rlog (Love, Huber, and Anders 2014) to normalise raw data prior to createing PCA and checking sample distance. An alternative transformation is vst, which is faster)
specified blind = FALSE, which means that differences between cell lines and treatment (the variables in the design) will not contribute to the expected variance-mean trend of the experiment.
Don't use this transform in the DESeq analysis this because For differential testing the DESeq function applied to raw counts is recommended

```{r rlog_trans}

rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)
```

```{r plot_rlog_trans_effects}

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
```

```{r sample_distances}

sampleDists <- dist(t(assay(rld)))
```

```{r heatmap_of_sample_dist}
sampleDistMatrix <- as.matrix( sampleDists )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf("plots/STR_vs_INF_heat_map_sample_dist_rld-normalised_data.pdf", onefile = FALSE)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

dev.off()
```

Plot PCA of rlog transformed data with plotPCA
```{r PCA}

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
```


#Run DESeq algorithm

If you never tell the DESeq2 functions which level you want to compare against (e.g. which level represents the control group), the comparisons will be based on the alphabetical order of the levels. Here just set the reference condition to wild type
```{r set_ref_condition}

dds$condition <- relevel(dds$condition, ref = "INF")
```

# DESeq consists of three steps
   estimating size factors (gene level differences, these will be overridden by avgTxLength from Tximport)
   estimating dispersions (differences between biological replicates)
   fit data to negative binomial and conduct Wald test

```{r run_DESeq}

dds <- DESeq(dds)
res <- results(dds)
```

```{r diff_gene_exp}

summary(res) #summary of results
```

```{r summary_sorted_by_p}

res <- res[order(res$padj),]
head(res)
```

```{r plot_top_siz_diffs}

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
```

```{r plot_dispersions}
pdf("plots/STR_vs_INF_plot_dispersions.pdf", onefile = FALSE)
plotDispEsts(dds)
dev.off()
```

Prior to plotting DE's shrink the data with apeglm
```{r shrinkage}

# shrink results using thresholds of log2FC of 1 or 2
resApeT_thresh1  <- lfcShrink(dds, coef="condition_STR_vs_INF", type="apeglm", lfcThreshold = 1)
resApeT_thresh2  <- lfcShrink(dds, coef="condition_STR_vs_INF", type="apeglm", lfcThreshold = 2)


pdf("plots/STR_vs_INF_apeglm_shrinkage_MA_plot.pdf", onefile = FALSE)
drawLines <- function() abline(h=c(-2,2),col="dodgerblue",lwd=2)
plotMA(resApeT_thresh2, ylim=c(-8,8), cex=.8); drawLines()
dev.off()

colnames(resApeT_thresh1)  <- c("baseMean", "apeglm_1_log2FC", "apeglm_1_lfcSE", "apeglm_1_svalue" )
colnames(resApeT_thresh2)  <- c("baseMean", "apeglm_2_log2FC", "apeglm_2_lfcSE", "apeglm_2_svalue" )
```


```{r volcano_plot}

pdf("plots/STR_vs_INF_volcano_plot.pdf", onefile = FALSE)

#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

dev.off()

```


```{r summarise_results}

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



write.table(goi, "result_tables/E.festucae_E2368_STR_INF.txt", quote = FALSE, row.names = FALSE, sep = "\t")
```

















################  ARCHIVE

```{r look_at_dds, eval=FALSE, include=FALSE}

# re-set dds
dds <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)
# collapse technical replicates
dds <- collapseReplicates(dds, groupby = samples$sample, run = samples$run, renameCols = TRUE)

assays(dds)
head(assays(dds)[["counts"]])   # from Tximport
head(assays(dds)[["avgTxLength"]])  # from Tximport

dds <- estimateSizeFactors(dds)
assays(dds)
head(assays(dds)[["normalizationFactors"]])  # from estimateSizeFactors

# extract normalised counts
normalized_counts <- as.data.frame(counts(dds, normalized=TRUE))
normalized_counts$gene_id <- rownames(normalized_counts)
head(normalized_counts)

dds <- estimateDispersions(dds)
assays(dds)
head(assays(dds)[["mu"]])   # from estimateDispersions

dds <- nbinomWaldTest(dds)
assays(dds)
head(assays(dds)[["H"]])    # from nbinomWaldTest (hat matrix diagonals)
head(assays(dds)[["cooks"]])   # from nbinomWaldTest (cooks distance - identifies samples that are outliers if there are three or more samples. DESeq will automatically act on this information if there are 7 or more samples/biological replicates)
# note TPM is not taken in by DESeq and is only available prior to technical replicate collapse
```


filter out data with less than 10 read counts - skip this, Barry and Daniel want to see reads that are present in one sample and not in another (most of the secondary metabolites show this pattern)
```{r filter_out_low_counts, eval=FALSE, include=FALSE}

#keep <- rowSums(counts(dds)) >= 10
#dds <- dds[keep,]
```

```{r check_normalisation, eval=FALSE, include=FALSE}

head(normalizationFactors(dds))
head(assays(dds)[["avgTxLength"]])   # gene level normalisation for length, hex primer and (beta) gc content
head(assays(dds)[["mu"]])            # dispersion level between samples of the same condition
head(counts(dds))
head(counts(dds, normalized = TRUE))

#7411.1042*0.9524896      # normalised count X normalisation factor = raw count
```

look at effects of shrinkage
```{r effects_of_apeglm_shrinkage, eval=FALSE, include=FALSE}
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("apeglm", version = "3.8")
library(apeglm)

resultsNames(dds)

resGA <- results(dds, lfcThreshold=2, altHypothesis="greaterAbs")
resApeT <- lfcShrink(dds, coef=2, type="apeglm", lfcThreshold = 2)

resGA
resApeT


pdf("plots/clrD_vs_wt_apeglm_shrinkage_MA_plot.pdf", onefile = FALSE)
par(mfrow=c(2,1),mar=c(2,2,1,1))
drawLines <- function() abline(h=c(-2,2),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=c(-8,8), cex=.8); drawLines()
plotMA(resApeT, ylim=c(-8,8), cex=.8); drawLines()
dev.off()
```

```{r alt_hypothesis, eval=FALSE, include=FALSE}

par(mfrow=c(2,2),mar=c(2,2,1,1))
ylim <- c(-4,4)
resGA <- results(dds, lfcThreshold=2, altHypothesis="greaterAbs")
resLA <- results(dds, lfcThreshold=2, altHypothesis="lessAbs")
resG <- results(dds, lfcThreshold=2, altHypothesis="greater")
resL <- results(dds, lfcThreshold=2, altHypothesis="less")
drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()
plotMA(resLA, ylim=ylim); drawLines()
plotMA(resG, ylim=ylim); drawLines()
plotMA(resL, ylim=ylim); drawLines()
```

Visualise counts in a heat map
```{r heat_map_of_count_matrix, eval=FALSE, include=FALSE}
rld <- rlog(dds, blind=FALSE)
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("sample", "condition")])
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
```

```{r add_raw_count_data, eval=FALSE, include=FALSE}
# merge in count data from counts(dds) <<<<<<<<  found some high log2 fold results that came from very low counts in one or two samples
count_data <- as.data.frame(counts(dds))
count_data$gene_id <- rownames(count_data)
group colnames by config file 'condition' and get average.
conditions <- as.vector(unique(samples$condition))
mean_counts <- sapply(conditions, function(x) as.integer(rowMeans(count_data[,colnames(count_data) %in%
                                                    as.vector(unique(samples[samples$condition == x, "sample"]))], na.rm=TRUE)))
colnames(mean_counts) <- paste(colnames(mean_counts), "counts", sep = "_")
mean_counts <- as.data.frame(base::cbind(gene_id = count_data$gene_id, mean_counts, deparse.level = 0))
all_results <- merge(all_results, mean_counts, by="gene_id")

```

```{r add_DESeq2_normalised_count_data}

# extract normalised counts
normalized_counts <- as.data.frame(counts(dds, normalized=TRUE))
normalized_counts$gene_id <- rownames(normalized_counts)

# merge in mean normalised count data
# i.e. group colnames by config file 'condition' and get average count data.
conditions <- as.vector(unique(samples$condition))
mean_counts <- sapply(conditions, function(x) as.integer(rowMeans(normalized_counts[,colnames(normalized_counts) %in%
                                                    as.vector(unique(samples[samples$condition == x, "sample"]))], na.rm=TRUE)))
colnames(mean_counts) <- paste(colnames(mean_counts), "norm_counts", sep = "_")
mean_counts<- as.data.frame(base::cbind(gene_id = normalized_counts$gene_id, mean_counts, deparse.level = 0))
all_results <- merge(all_results, mean_counts, by="gene_id")
```

```{r merge_fpkm, eval=FALSE, include=FALSE}

# merge in mean fpkm count data
fpkm_counts <- as.data.frame(fpkm(dds))
fpkm_counts$gene_id <- rownames(fpkm_counts)
fpkm_counts
mean_fpkm <- sapply(conditions, function(x) as.integer(rowMeans(fpkm_counts[,colnames(fpkm_counts) %in%
                                                    as.vector(unique(samples[samples$condition == x, "sample"]))], na.rm=TRUE)))
colnames(mean_fpkm) <- paste(colnames(mean_fpkm), "fpkm", sep = "_")
mean_fpkm <- as.data.frame(base::cbind(gene_id = fpkm_counts$gene_id, mean_fpkm, deparse.level = 0))
all_results <- merge(all_results, mean_fpkm, by="gene_id")
```

```{r merge_DESeq_normalised_counts, eval=FALSE, include=FALSE}

# extract normalised counts
normalized_counts <- as.data.frame(counts(dds, normalized=TRUE))
normalized_counts$gene_id <- rownames(normalized_counts)

# merge in mean normalised count data
# i.e. group colnames by config file 'condition' and get average count data.
conditions <- as.vector(unique(samples$condition))
mean_counts <- sapply(conditions, function(x) as.integer(rowMeans(normalized_counts[,colnames(normalized_counts) %in%
                                                    as.vector(unique(samples[samples$condition == x, "sample"]))], na.rm=TRUE)))
colnames(mean_counts) <- paste(colnames(mean_counts), "norm_counts", sep = "_")
mean_counts<- as.data.frame(base::cbind(gene_id = normalized_counts$gene_id, mean_counts, deparse.level = 0))
all_results <- merge(all_results, mean_counts, by="gene_id")
```

```{r flag_significant_results, eval=FALSE, include=FALSE}

# flag results that have an adjusted pvalue less than .01 and a log2FoldChange >1
targets <- rownames(subset(GOI, padj<.01 & abs(log2FoldChange)>=1 & !is.na(log2FoldChange) ))
GOI$log2FC_1 <- ifelse(rownames(GOI) %in% targets, 1, 0)

# flag results that have an adjusted pvalue less than .01 and a log2FoldChange >2
targets <- rownames(subset(GOI, padj<.01 & abs(log2FoldChange)>=2 & !is.na(log2FoldChange) ))
GOI$log2FC_2 <- ifelse(rownames(GOI) %in% targets, 1, 0)
```

```{r add_annotations, eval=FALSE, include=FALSE}

# do this after defining core gene set

##### GET ANNOTATIONS FROM GFF3 FILE
positions <- read.delim("../../0_raw_data/genome/EfeFl.gff3", header = FALSE, comment.char = "#", blank.lines.skip = TRUE, fill = TRUE, sep = "\t")
colnames(positions) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

# select out mRNA (these are the ones with product described)
mRNA <- positions[positions$type == "mRNA",]

# split attributes to get c("trans_id", "gene_id", "product")
temp <- str_split(mRNA$attributes, ";", n=3)
df <- as.data.frame(t(as.data.frame(temp)))
rownames(df)<-NULL
mRNA[, c("trans_id", "gene_id", "product")] <- df
mRNA$trans_id <- str_replace(mRNA$trans_id, "ID=", "")
mRNA$gene_id <- str_replace(mRNA$gene_id, "Parent=", "")
mRNA$product <- str_replace(mRNA$product, "product=", "")
mRNA$attributes <- NULL
mRNA$score <- NULL
mRNA$phase <- NULL
head(mRNA)
all_results <- merge(tempdf, mRNA, by = "gene_id")

#### GET ANNOTATIONS FROM annotations.txt file from David's funannotate pipeline file
annotations <- read.delim("../../0_raw_data/genome/EfeFl_annotations.txt", header = TRUE, comment.char = "#", blank.lines.skip = TRUE, fill = TRUE, sep = "\t")
head(positions)

# merge DESeq output (res) and GFF3 data for mRNA (where the annotations are)
all_results <- merge(tempdf, annotations, by.x = "gene_id", by.y = "GeneID")

```
