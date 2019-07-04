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
my $maxsvalue = 0.005;
my $minlog2fcEXTN = 0.5;
my $maxsvalueEXTN = 0.001;
my $gaps_allowed = 1;
my %raw;


# read in data to hash
while(<IN>) {
  chomp;
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
    #print "spp: $spp\tsvalue_2: $svalue_2\tgene_id:$raw{$spp}{$contig}{$start}{'gene_id'}\n";                #:)
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
my $current_start;
my $gapflag = 0;

# identify initial genes in clusters
foreach my $spp (sort keys %raw){
  foreach my $chromosome (sort keys %{$raw{$spp}}){
    #print "chromosome: $chromosome\n";                     #### ;)
    my $prev = 0; my $svalue_1; my $svalue_2;
    foreach my $start_pos (sort keys %{$raw{$spp}{$chromosome}}){
      $raw{$spp}{$chromosome}{$start_pos}{'prev'} = $prev;
      #$prev = $current_start;
      $current_start = $start_pos;
      $l2fc = $raw{$spp}{$chromosome}{$start_pos}{'l2fc'};
      $svalue_1 = $raw{$spp}{$chromosome}{$start_pos}{'svalue_1'};
      $svalue_2 = $raw{$spp}{$chromosome}{$start_pos}{'svalue_2'};
      print "$spp\t$chromosome\t$start_pos\t$l2fc\t$svalue_1\n";
      # check whether to initiate cluster
      if ( ( ($l2fc  > 0) && ($l2fc > $minlog2fc) ) || \\
           ( ($l2fc  < 0) && ($l2fc  < -$minlog2fc) ) &&
           ($svalue_1 < $maxsvalue) ) {
        $raw{$spp}{$chromosome}{$start_pos}{'clust'} = 1;
        $gapflag = 0;
        # if you find an initial cluster, roll back to see if previous ones can be included
        while ($gapflag <= $gaps_allowed) {
          print "prev: $prev\n";
          $l2fc = $raw{$spp}{$chromosome}{$prev}{'l2fc'};
          $svalue_1 = $raw{$spp}{$chromosome}{$prev}{'svalue_1'};
          $svalue_2 = $raw{$spp}{$chromosome}{$prev}{'svalue_2'};
          $prev = $raw{$spp}{$chromosome}{$prev}{'prev'};
          print "prev values: $l2fc\t$svalue_1\tnext previous start: $prev\n";
          if ( ($svalue_1 < $maxsvalueEXTN)  && ( ( ($l2fc  > 0) && ($l2fc > $minlog2fcEXTN) ) || \\
          ( ($l2fc  < 0) && ($l2fc  < -$minlog2fcEXTN) ) ) ) {
            if ($l2fc  > 0) {print "greater than 0\n";}
            if ($l2fc > $minlog2fcEXTN)  {print "greater then min log2fc for extn $minlog2fcEXTN)\n";}
            if ($l2fc  < 0) {print "less than 0\n";}
            if ($l2fc  < -$minlog2fcEXTN) {print "less than negative min log2fc for extn -$minlog2fcEXTN\n";}
            if ($svalue_1 < $maxsvalueEXTN) {print "svalue $svalue_1 is less than the max for extn $maxsvalueEXTN\n";}
            print "extended $l2fc\t$svalue_1\n";
            print "svalue is $svalue_1, max svalue for extension is $maxsvalueEXTN\n";
            exit;
            $raw{$spp}{$chromosome}{$prev}{'ext'} = 1;
          }
          else {
            $gapflag ++;
            $raw{$spp}{$chromosome}{$prev}{'ext'} = 1;
            $l2fc = $raw{$spp}{$chromosome}{$prev}{'l2fc'};
            $svalue_1 = $raw{$spp}{$chromosome}{$prev}{'svalue_1'};
            $svalue_2 = $raw{$spp}{$chromosome}{$prev}{'svalue_2'};
            $prev = $raw{$spp}{$chromosome}{$prev}{'prev'};
            print "fell below threshold, gapflag = $gapflag\n";
          }
        }
        # forward through sequences for extending clusters
        # reset start position
      }
      else {
        $prev = $current_start;
        print "No cluster intitiation prev: $prev\n";
      }
    }
  }
}
