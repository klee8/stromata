# map reads to genome with STAR

#!/usr/bin/env bash

# make a genome index
#STAR --runMode genomeGenerate --runThreadN 4 --genomeDir genomes/Eel_728 --genomeLoad LoadAndKeep --genomeFastaFiles genomes/Eel_728/Eel_masked.fasta --outFileNamePrefix Eel_728_
#STAR --runMode genomeGenerate --runThreadN 4 --genomeDir genomes/Efe_2368 --genomeLoad LoadAndKeep --genomeFastaFiles genomes/Eel_728/Eel_masked.fasta --outFileNamePrefix Efe_2368_
#STAR --runMode genomeGenerate --runThreadN 4 --genomeDir genomes/Ety_8 --genomeLoad LoadAndKeep --genomeFastaFiles genomes/Eel_728/Eel_masked.fasta --outFileNamePrefix Ety_8_


# map reads to genome
