# id_DE_clusters.pl Kate D Lee July 2019
# usef file with columns for:
#  species, chr, start, end, orthogroup, log2fc, log2fcSE, svalue
# identify clusters in each species with DE info and positional info
# use ortholog info to compare clusters across species
#############################################################

#!usr/bin/perl -w
use strict;

my $file = $ARGV[0] || "sub_core_gene_set_STR_PS_ann.txt";
open(IN, "<$file") || die "ERROR: couldn't open file: $!";

my $minlog2fc = 1;
my $maxsvalue = 0.001;
my $minlog2fcEXTN = 0.5;
my $maxsvalueEXTN = 0.001;
my $gaps_allowed = 1;
my %raw;


# read in data to hash
while(<IN>) {
  if ($_=~ m/orthogroup/) {
      next;
  }
  else {
    my ($contig, $start, $end, $gene_id, $orthogroup, $l2fc, $l2fcSE, $svalue_1, $svalue_2, $spp) = split("\t", $_);
    $raw{$spp}{$contig}{$start}{'end'} = $end;
    $raw{$spp}{$contig}{$start}{'gene_id'} = $gene_id;
    $raw{$spp}{$contig}{$start}{'orthogroup'} = $orthogroup;
    $raw{$spp}{$contig}{$start}{'l2fc'} = $l2fc;
    $raw{$spp}{$contig}{$start}{'l2fcSE'} = $l2fcSE;
    $raw{$spp}{$contig}{$start}{'svalue_1'} = $svalue_1;
    $raw{$spp}{$contig}{$start}{'svalue_2'} = $svalue_2;
  }
}


my $spp = "";
my $chr = "";
my $pos = 0;
my $gene_id = "";
my $orthogroup = "";
my $l2fc;
my $l2fcSE;
my $svalue_1;
my $svalue_2;
my $prev;
my $gapflag = 0;

# identify initial genes in clusters
for my $species (sort keys %raw){
  for my $chromosome (sort keys %{$raw{$spp}}){
    my $prev = 0; my $svalue_1; my $svalue_2;
    for my $start_pos (sort keys %{$raw{$spp}{$chromosome}}){
      $raw{$spp}{$chromosome}{$start_pos}{'prev'} = $prev;
      $prev = $start_pos;
      $l2fc = $raw{$spp}{$chromosome}{$start_pos}{'l2fc'};
      $svalue_1 = $raw{$spp}{$chromosome}{$start_pos}{'svalue_1'};
      $svalue_2 = $raw{$spp}{$chromosome}{$start_pos}{'svalue_2'};
      if ( ( ($l2fc  > 0) && ($l2fc > $minlog2fc) ) || \\
           ( ($l2fc  < 0) && ($l2fc  < -$minlog2fc) ) &&
           ($svalue_1 < $maxsvalue) ) {
        $raw{$spp}{$chromosome}{$start_pos}{'clust'} = 1;
        $l2fc = $raw{$spp}{$chromosome}{$prev}{'l2fc'};
        $svalue_1 = $raw{$spp}{$chromosome}{$prev}{'svalue_1'};
        $svalue_2 = $raw{$spp}{$chromosome}{$prev}{'svalue_2'};
        $prev = $raw{$spp}{$chromosome}{$prev}{'prev'};
        $gapflag = 0;
        # if you find an initial cluster, roll back to see if previous ones can be included
        while ($gapflag <= $gaps_allowed) {
          if ( ( ($l2fc  > 0) && ($l2fc > $minlog2fcEXTN) ) || \\
          ( ($l2fc  < 0) && ($l2fc  < -$minlog2fcEXTN) ) && \\
          ($svalue_1 < $maxsvalueEXTN) ) {
            $raw{$spp}{$chromosome}{$prev}{'ext'} = 1;
            $l2fc = $raw{$spp}{$chromosome}{$prev}{'l2fc'};
            $svalue_1 = $raw{$spp}{$chromosome}{$prev}{'svalue_1'};
            $svalue_2 = $raw{$spp}{$chromosome}{$prev}{'svalue_2'};
            $prev = $raw{$spp}{$chromosome}{$prev}{'prev'};
          }
          else {
            $gapflag ++;
            $raw{$spp}{$chromosome}{$prev}{'ext'} = 1;
            $l2fc = $raw{$spp}{$chromosome}{$prev}{'l2fc'};
            $svalue_1 = $raw{$spp}{$chromosome}{$prev}{'svalue_1'};
            $svalue_2 = $raw{$spp}{$chromosome}{$prev}{'svalue_2'};
            $prev = $raw{$spp}{$chromosome}{$prev}{'prev'};
          }
        }
      }
    }
  }
}
