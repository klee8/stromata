# setup RNA-seq analysis for Epichloe spp

# require config file with info on fastq files grouped in biological replicates by tissue and spp
# config file with 7 columns (as in DESeq), no headers in file:
# species centre  assay run_num(when samples sent for re-sequencing) lane sample  experiment  run(illumina run fq base filename) condition(tissue)
# E.elymi_NfE728	NA	ILLUM	first_run	Lane6	E.elymi_NfE728_INF1	STROMA	Ee-Inf1_S14_L006_R1_001	INF

#!usr/bin/perl
use strict;
use warnings;


my $config = $ARGV[0];
my $deseq_compare = $ARGV[1];
my $raw_dir = $ARGV[2] || "/opt/stroma_RNAseq";


open(IN, "<$config") || die "ERROR, please provide config file: $!";

my %config;
while(<IN>){
  #print $_;
  chomp;
  my ($species, $centre, $assay, $run_num, $lane, $sample, $experiment, $run, $condition) = split("\t", $_);
  $config{$species}{$run}{$condition} = $_;
  $config{$species}{$run}{'condition'} = $condition;
  $config{$species}{$run}{'run_number'} = $run_num;
  $config{$species}{$run}{'lane'} = $lane;
}

################################# make Quality control folders and set up QC scripts
`mkdir 1_QC`;

open(FQCR, ">1_QC/fastqc_raw.sh") || die "ERROR, couldn't open 1_QC/fastqc_raw.sh: $!";
print FQCR "# fastQC of reads prior to trimming

#!/usr/bin/env bash
source /home/kate/.bash_profile
mkdir fastQC_raw\n";

open(FQCT, ">1_QC/fastqc_trimmed.sh") || die "ERROR, couldn't open 1_QC/fastqc_trimmed.sh: $!";
print FQCT "# fastQC of reads after trimming

#!/usr/bin/env bash
source /home/kate/.bash_profile
mkdir fastQC_trimmed\n";

open(TRIM, ">1_QC/trimmomatic.sh") || die "ERROR, couldn't open 1_QC/fastqc_trimmed.sh: $!";
print TRIM "# trim reads including adapters

#!/usr/bin/env bash
source /home/kate/.bash_profile

# Trim data

# - trimlog records readname and lengths before and after trimming
#ILLUMINACLIP:TruSeq2-PE.fa:2:30:10   <fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>
#SLIDINGWINDOW:5:20   <windowSize>:<requiredQuality>
#LEADING:5         removes low quality bases from the start of a read
#TRAILING:5        removes low quality bases from the end of a read
#MINLEN:40         removes reads shorter than this
# note: stopped outputting the trimlog, I don't use it and it takes up a LOT of space.
# note: stdin was overwriting itself in the fastqc_raw/\$species output directory - added \$run as an extra folder in the path
mkdir trimmed
\n";

foreach my $species (sort keys %config){
  print FQCR "mkdir fastQC_raw/$species\n";
  print FQCT "mkdir fastQC_trimmed/$species\n";
  print TRIM "mkdir trimmed/$species trimmed/$species/logs\n";
  foreach my $run (sort keys %{$config{$species}}){
    my $cond = $config{$species}{$run}{'condition'} ;
    my $run_number = $config{$species}{$run}{'run_number'} ;
    my $lane_num = $config{$species}{$run}{'lane'} ;
    print FQCR "zcat $raw_dir/$run_number/$lane_num/$run.fastq.gz  | fastqc stdin --outdir=fastQC_raw/$species/$run\n";
    print FQCT "fastqc trimmed/$species/$run.trim.fastq.gz --outdir=fastQC_trimmed/$species\n";
    print TRIM "java -jar /home/kate/bin/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 4 $raw_dir/$run_number/$lane_num/$run.fastq.gz trimmed/$species/$run.trim.fastq.gz ILLUMINACLIP:/home/kate/bin/Trimmomatic-0.38/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:40\n";
  }
}

close FQCR;
close FQCT;
close TRIM;


################################# setup Salmon mapping and quant
`mkdir 2_salmon`;

open(QUANT, ">2_salmon/salmon_quant.sh") || die "ERROR, couldn't open salmon_quant file: $!";


print QUANT "# Salmon mapping and read count quantification
# make sure Salmon software is in the PATH

#!/usr/bin/env bash
source /home/kate/.bash_profile

#
mkdir salmon_quant

# Salmon usage
# build an indices of the transcriptomes
# -k 31 suggested for reads over 75bp
salmon index -t transcriptome/E.elymi_NfE728/Epichloe_elymi_NfE728.transcripts.fa -i transcriptome/E.elymi/E.elymi_NfE728.trans_index -k 31
salmon index -t transcriptome/E.festucae_E2368/Epichloe_festucae_E2368.transcripts.fa -i transcriptome/E.festucae/E.festucae_E2368.trans_index -k 31
salmon index -t transcriptome/E.typhina_E8/Epichloe_typhina_E8.transcripts.fa -i transcriptome/E.typhina_E8/E.typhina_E8.trans_index -k 31


# grab headers from transriptome files for use in mapping transcripts to genes later
for i in transcriptome/*; do grep '>' \$i/*.fa > \$i/transcriptome.headers.txt; sed -i s'/>//g' \$i/transcriptome.headers.txt; done

##########  QUANT OPTION (quasi-mapping and counting)
# -l U for single-end (unstranded) librarytype
# -r for SE fq file
# --validateMappings makes mapping more sensitive
# --seqBias corrects for random hexamer bias in pcr priming
# --gcBias
# --posBias enables modeling of a position-specific fragment start distribution (bias from position in read)
# --writeUnmappedNames   # this takes up loads of space - skip it
\n";

foreach my $species (sort keys %config){
  print QUANT "mkdir salmon_quant/$species\n";
  foreach my $run (sort keys %{$config{$species}}){
    my $cond = $config{$species}{$run}{'condition'} ;
    print QUANT "salmon quant -p 8 -i transcriptome/$species/$species.trans_index -l U -r <(zcat ../1_QC/trimmed/$species/$run.trim.fastq) --validateMappings --seqBias --gcBias --posBias -o salmon_quant/$species/$run\n";
  }
}

print "CAUTION: Ensure the transcriptomes in 2_salmon/salmon_quant.sh script match the species you are looking at\n";

################################# setup DEseq2
`mkdir 3_DESeq2`;

# read in list of DESeq2 comparisons to make
# file with one line for each comparison with three tab delimited columns
# (species, condition to compare (SOI), reference condition (REF))
open(DESEQ, "<$deseq_compare") || die "ERROR, couldn't open $deseq_compare: $!";
my %comparisons;
while(<DESEQ>) {
  chomp;
  my ($spp, $cond, $ref) = split("\t", $_);
  my $job = join("_", $cond, $ref);
  $comparisons{$spp}{$job} = 1;
}

# set up file structure and print Rmd files for each comparison
open(TEMP, ">3_DESeq2/temp.sh") || die "ERROR, couldn't open temporary file for DESeq2: $!";
foreach my $species (sort keys %comparisons){
  print TEMP "mkdir $species\n";
  foreach my $job (sort keys %{$comparisons{$species}}) {
    my ($cond, $ref) = split("_", $job);
    print "COND = $cond; REF = $ref\n";
    print TEMP "cp ../Generic_DESeq2.Rmd $species/$species.DESeq2.$cond.$ref.Rmd\n";
    print TEMP "sed -i \"s/SPECIES/$species/g\" $species/$species.DESeq2.$cond.$ref.Rmd\n";
    print TEMP "sed -i \"s/SOI/$cond/g\" $species/$species.DESeq2.$cond.$ref.Rmd\n";
    print TEMP "sed -i \"s/REF/$ref/g\" $species/$species.DESeq2.$cond.$ref.Rmd\n";
  }
}
