# format_SMURF_output.#!/usr/bin/env perl
# Kate Lee 2019
# Creates a table that can be read into R
#
# usage: perl format_SMURF_output.pl <SMURF_output_file> <reformatted_filename>
###########################################################

#!usr/bin/perl
use strict;
use warnings;


my $smurf = $ARGV[0];
my $outfile = $ARGV[1];

open(SMURF, "<$smurf") || die "ERROR: cannot open $smurf: $!";
open(OUT, ">$outfile") || die "ERROR: cannot open $outfile: $!";

print OUT "Cluster\tBackbone_gene_id\tGene_id\tGene_positions\tChromosome-Contig\t";
print OUT "Gene_order\t5_end\t3_end\tGene_distance\tDomain_score\t";
print OUT "Annotated_gene_function\n";


my $cluster;
while(<SMURF>) {
  $_ =~ s/\r/\t/g;
  chomp;
  my $backbone;	my $gene_id;	my $gene_pos;
  my $chr_contig; my $gene_order; my $five_end; my	$three_end;
  my $gene_dist; my	$domain_score; my	$ann;
  if ($_ =~ /Backbone_gene_id/) { next; }
  elsif ($_ =~ /^\s+/) { next; }
  elsif ($_ =~ /^Cluster/)
  {

      my @temp = split(":", $_);
      $cluster = $temp[1];
      $cluster =~ s/\t//g;
  }
  else
  {
    ($backbone, $gene_id, $gene_pos, $chr_contig, $gene_order, $five_end, $three_end, $gene_dist, $domain_score, $ann) = split("\t", $_);
    print OUT "$cluster\t$backbone\t$gene_id\t$gene_pos\t$chr_contig\t$gene_order\t$five_end\t$three_end\t$gene_dist\t$domain_score\t$ann\n";
  }
}
