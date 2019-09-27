## Identify shared DE genes and co-regulated gene clusters in differentially expressed data

#### Kate Lee 2019


### Aim:
Identify core set of genes that are differentially expressed (DE) in stromata versus stem
Identify clusters of DE genes that are shared between species.

### Run:
    bash run_reformat_core_gene_set.sh 
    bash run_DE_permutations.sh       # permutes foldchange results within specices and checks for core gene sets
    bash run_cluster_ID.sh            # outputs clusters and python graph scripts
    bash run_cluster_permutations.sh  # same as run_cluster_ID.sh with permutations also

### Reformatting:
Input has columns for each gene in analysis (reformatted using reformat_core_gene_set.Rmd):  
"contig", "start", "stop", "gene_id","orthogroup", "log2fc", "lfcSE", "svalue_1", "svalue_2", "species"

### shared DE genes
Number of DE genes shared between all species 
- 'top_2FC' : >= 2-fold-change in all species
- 'top_4FC' : >= 4-fold-change in all speices### Identify clusters:
- 'core_genes' : at least 2-fold-change in all species at least 4-fold-change in one species* Number of clusters found in each species

### Co-regulated gene clusters 
- consecutive genes with >= 2-fold-change DE
- either 3 genes in a row, or more than three genes with one gap.

### Permutation tests
R scripts which randomise DE values among genes for each species and then test for 
- shared DE genes
- co-regulated gene clusters




    
    
    
