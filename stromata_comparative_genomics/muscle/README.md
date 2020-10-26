
### elsewhere:
+ in exonerate folder: 
+ in blastn folder: Took E. elymi NfE728 nuclotide sequences for core genes (positions from gff file) and used blastn to find homologs in all genomes.
+ added new gene models to E. elymi gff3 file

### setup files for alignment
+ note change list of genes in for loop as appropriate (e.g. core_gene_list, control_list, last_10)
+ note scripts assume 'gene' annotation will be the first annotation for each gene in the gff file (note- genious outputs CDS first)

```
bash setup_muscle.sh
```
+ Preferentially took exonerate results, but where these were missing took blastn results (previously used tblastn of protein sequence - see old alignments)

or

```
bash 
```
+ used only blastn results for alignment setup



### align sequence with muscle 
```
bash run_muscle.sh
```
+ aligned with muscle (just runs last group of genes that were setup using one of the above)


### annotation
+ note change list of genes in 3 for loop as appropriate (e.g. core_gene_list, control_list, last_10)

```
bash exon_info_blastn.sh
```
+ annotated exons with E. elymi NfE728 gff file
+ added 99 nt upstream and downstream of aligned sequence (to be able to to check for start and stop codons that are outside of the E. elymi NfE728 sequences)
+ results in blastn_ann_alns 



#### TO DO
Niceties: automate file in for loop from bash command line
