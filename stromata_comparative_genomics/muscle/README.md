
### elsewhere:
+ in exonerate folder: 
+ in blastn folder: Took E. elymi NfE728 nuclotide sequences for core genes (positions from gff file) and used blastn to find homologs in all genomes.


### setup files for alignment
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
```

```
+ annotated exons with E. elymi NfE728 gff file
+ added 99 nt upstream and downstream of aligned sequence (to be able to to check for start and stop codons that are outside of the E. elymi NfE728 sequences)

