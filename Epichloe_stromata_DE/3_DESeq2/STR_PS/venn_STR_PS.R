# Kate October 2020
# create a venn diagram of stromata DE genes
# adapted from https://flowingdata.com/2018/01/10/how-to-make-venn-diagrams-in-r/

#install.packages("VennDiagram")
library(VennDiagram)

# could use polyclip afterwards to manually set colours of intersections
#install.packages("polyclip")
#library(polyclip)


# grab data
data <- read.delim("core_set_STR_PS.txt", header = TRUE, stringsAsFactors = FALSE, sep = " ", na.strings = "NA")
head(data)
str(data)
all_int <- length(which(data$num_spp == 3))
ely_fes <- nrow(data[(data$num_spp == 2) & (data$typhina_gene_id == "NA"),]) + all_int
fes_typ <- nrow(data[(data$num_spp == 2) & (data$elymi_gene_id == "NA"),]) + all_int
typ_ely <- nrow(data[(data$num_spp == 2) & (data$festucae_gene_id == "NA"),]) + all_int
nrow(data[(data$num_spp == 1), ])
nrow(data[(data$num_spp == 1) & !is.na(data$elymi_gene_id),])
ely <- nrow(data[!is.na(data$elymi_gene_id),])
fes <- nrow(data[!is.na(data$festucae_gene_id),])
typ <- nrow(data[!is.na(data$typhina_gene_id),])


# Construct Venn
# Triple Venn - labeled from topleft clockwise
# first three numbers set total in each circle from top left clockwise
# second three numbers set total in each pair from top center, clockwise
# last number is the number at the intersection between all three
svg(filename = "orth_venn.svg" ,
    width = 7, height = 7
)
plot.new()
categories <- c("Epichloe elymi", "Epichloe festucae", "Epichloe typhina")
purplehues <- c("#8856a7", "#9ebcda", "#e0ecf4")
greenhues <- c("#edf8b1", "#7fcdbb", "#2c7fb8")
venn3 <- draw.triple.venn(ely, fes, typ, ely_fes, fes_typ, typ_ely, all_int, category = categories, 
                          fill = greenhues, col = NA, euler.d = TRUE, scaled = TRUE, 
                          fontface = rep("plain", 7), fontfamily = rep("serif", 7),
                          cat.fontface = rep("plain", 3), cat.fontfamily = rep("serif", 3))
dev.off()

