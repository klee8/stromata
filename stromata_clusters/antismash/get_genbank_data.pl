# Kate Lee June 2019 get_genbank_data.pl
# retrive region, start and end of region, most similar known cluster, description, similarity%,
# gene_name, gene_start, gene_end, len_ex_introns, antiSMASH_gene_type, aa_seq
# usage: perl get_genbank_data.pl <antismash_file.gbk>

#!/usr/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;

my $usage = "get_genbank_data.pl <file>";
my $file = shift or die $usage;
my $format = "genbank";

# antismash puts region number into the filename
my $region = $file;
$region =~ s/\.gbk//;
$region =~ s/.*region//;

# antismash adds region location on chr into COMMENT section of genbank file
my $origin_start = "";
my $origin_end = "";
open(IN, "<$file") || die "ERROR: couldn't open $file: $!";
while(<IN>){
  chomp;
  if ($_=~ m/Orig.*start/){
    $origin_start = $_ ;
    $origin_start =~ s/Orig. start  :: //;
    $origin_start =~ s/\s*//;
  }
  if ($_=~ m/Orig.*end/){
    $origin_end = $_;
    $origin_end =~ s/Orig. end    :: //;
    $origin_end =~ s/\s*//;
  }
  if ($_=~ m/##antiSMASH-Data-END##/){ last;}
}
close IN;

# read in genbank file to Bio::seqIO
my $inseq = Bio::SeqIO->new(-file   => "$file",
                            -format => $format, );

# use Bio::SeqIO to parse file for Features
while (my $seq = $inseq->next_seq) {
  for my $feat ($seq->get_SeqFeatures) {
    my $primary = $feat->primary_tag;
    if ($primary eq 'CDS') {
      print $seq->accession_number,"\t", $region, "\t", $origin_start, "\t", $origin_end, "\t";
      print $feat->get_tag_values('Name'), "\t";
      my $beg = "";
      my $end = "";
      foreach my $location ( $feat->location->each_Location ) {
        print $location->start, "..", $location->end, ",";
        if ($beg eq "") {$beg = $location->start; }
        $end = $location->end;
      }
      print "\t$beg\t$end\t";
      if ($feat->has_tag('gene_functions')) {
        print $feat->get_tag_values('gene_functions'), "\t";
      }
      else { print 'NA', "\t"; }
      if ($feat->has_tag('translation')) {
        print $feat->get_tag_values('translation');
      }
      else { print 'NA'; }
      print "\n";
    }
  }
}
