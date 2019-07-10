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
my $maxsvalueEXTN = 0.01;
my $gaps_allowed = 1;
my $distthresh = 10000;
my $mingenenum = 3;
my %raw;

print "maxsvalue: $maxsvalue\nminlog2fc: $minlog2fc\n";
print "maxsvalueEXTN: $maxsvalueEXTN\nminlog2fcEXTN: $minlog2fcEXTN\n";
print "distance threshold: $distthresh\n";
print "reading in data...\n";
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

#my $gene_id = "";
#my $l2fc;
#my $svalue_1;
#my $svalue_2;
my $prev;

print "calculating gene characteristics...\n";
# first pass to calculate characteristics of genes (cluster/extension thresholds; distances to next gene etc.)
foreach my $spp (sort keys %raw){
#  print "species: $spp\n";
  foreach my $contig (sort keys %{$raw{$spp}}){
#    print "\tchromosome: $contig\n";                     #### ;)
    my $prev = 0; my $svalue_1; my $svalue_2;
    foreach my $start_pos (sort { $a <=> $b } keys %{$raw{$spp}{$contig}}){
      $raw{$spp}{$contig}{$start_pos}{'prev'} = $prev;
      my $gene_id = $raw{$spp}{$contig}{$start_pos}{'gene_id'};
      my $l2fc = $raw{$spp}{$contig}{$start_pos}{'l2fc'};
      my $svalue_1 = $raw{$spp}{$contig}{$start_pos}{'svalue_1'};
      my $svalue_2 = $raw{$spp}{$contig}{$start_pos}{'svalue_2'};

      # check whether to initiate cluster
      if ( ($svalue_1 < $maxsvalue) && ( ( ($l2fc  > 0) && ($l2fc > $minlog2fc) ) || ( ($l2fc  < 0) && ($l2fc  < -$minlog2fc) ) ) )
      {
        $raw{$spp}{$contig}{$start_pos}{'clust'} = 2;
      }
      # check whether to extend cluster
      elsif ( ($svalue_1 < $maxsvalueEXTN)  && ( ( ($l2fc  > 0) && ($l2fc > $minlog2fcEXTN) ) || \\
            ( ($l2fc  < 0) && ($l2fc  < -$minlog2fcEXTN) ) ) )
      {
        $raw{$spp}{$contig}{$start_pos}{'clust'} = 1;
      }
      else
      {
        $raw{$spp}{$contig}{$start_pos}{'clust'} = 0;
      }

      # pos/neg DE
      if ($l2fc < 0) {$raw{$spp}{$contig}{$start_pos}{'posDE'} = 0;}
      if ($l2fc > 0) {$raw{$spp}{$contig}{$start_pos}{'posDE'} = 1;}

      # dist to prev gene
      $raw{$spp}{$contig}{$start_pos}{'distprev'} = $start_pos - $raw{$spp}{$contig}{$prev}{'end'};

      # reset for next gene
      $prev = $start_pos;

      #print "$gene_id\t$l2fc\t$svalue_1\t$raw{$spp}{$contig}{$start_pos}{'clust'}\n";
      #print "$gene_id\t$raw{$spp}{$contig}{$start_pos}{'clust'}\n";
    }
  }
}


my $clustsign = 10;   # used as a cluster flag and to track cluster sign
my %cluster;
my $initiate = 0;
my $clust_flag = 0;
my $clusternumber = 1;
my $gapflag = 0;
my %all_clusters;
my $clust_contig;
my $clust_start;
my $last_clust_end = 0;
my $gene_dist = 0;
my $larger_gap_needed;

open(CLUST, ">clusters.txt") || die "ERROR: couldn't open clusters.txt outfile: $!";
print CLUST "cluster\tcontig\tstart\tend\tdist\tgene_dist\tgenes\n";
print "identifying clusters...\n";

# identify initial genes in clusters
foreach my $spp (sort keys %raw){
  foreach my $contig (sort keys %{$raw{$spp}}){
    print "\n########### $contig\n";
    my $prev = 0; $gapflag = 0; $last_clust_end = 0; $gene_dist = 0;
    foreach my $start_pos (sort { $a <=> $b } keys %{$raw{$spp}{$contig}}){
      print "$start_pos\t";
      # check for gaps
      if ( ($clustsign != 10 ) && ($raw{$spp}{$contig}{$start_pos}{'clust'} == 0) ) { $gapflag++;}

      # flag genes that coud initiate a cluster
      if ($raw{$spp}{$contig}{$start_pos}{'clust'} == 2)
      {
        $initiate = 1;
      }

      # check if the contig has changed
      #if ($contig ne $clust_contig) { print "changed contig mid-cluster\n"; }

      # temp test print checks
      print "$raw{$spp}{$contig}{$start_pos}{'gene_id'}\t$raw{$spp}{$contig}{$start_pos}{'clust'}\t";
      print "posDE: $raw{$spp}{$contig}{$start_pos}{'posDE'}\tclustsign: $clustsign\tgap: $gapflag\t";
      print "dist: $raw{$spp}{$contig}{$start_pos}{'distprev'}\tdthresh: $distthresh\tcontig: $contig\tclust_contig: $clust_contig\t";

      # if you are in a cluster and pass more gaps than are allowed or
      # the distance to previous gene is larger than the thresholds
      # or the direction of DE sign changes
      # or you get to the next contig
      # evaluate potential cluster
      if ( ($clustsign != 10) && ( ($gapflag > $gaps_allowed) || ($raw{$spp}{$contig}{$start_pos}{'distprev'} > $distthresh) || ($clustsign != $raw{$spp}{$contig}{$start_pos}{'posDE'}) || ($contig ne $clust_contig) ) )
      {
        print "evaluate cluster...";

        # TRIM TRAILING GAPS
          foreach my $start_pos (sort { $b <=> $a } keys %cluster )
          {
            if ($raw{$spp}{$clust_contig}{$start_pos}{'clust'} == 0)
            {
              undef $cluster{$start_pos};
              print "\n\n#trimming\t$start_pos\n\n";
            }
            else {last;}
          }

        # if there is a gene that meets the 'intiate cluster thresholds'
        # and there are enough genes for a cluster
        if ( ($initiate == 1) && ($mingenenum <= keys %cluster) )
        {
          print "\n$clusternumber:\t";
          foreach my $start_pos (sort { $a <=> $b } keys %cluster )
          {
            push(@{$all_clusters{$spp}{$clusternumber}{'all_starts'}}, $start_pos);
            push(@{$all_clusters{$spp}{$clusternumber}{'all_gene_ids'}}, $raw{$spp}{$clust_contig}{$start_pos}{'gene_id'});
            $all_clusters{$spp}{$clusternumber}{'clust_dist'} = $clust_start - $last_clust_end;
            $all_clusters{$spp}{$clusternumber}{'gene_dist'} = $gene_dist;
            $all_clusters{$spp}{$clusternumber}{'contig'} = $clust_contig;
          }

          print CLUST "$clusternumber\t";
          print CLUST "$clust_contig\t";
          print CLUST "$clust_start\t";
          print CLUST "$raw{$spp}{$clust_contig}{$start_pos}{'end'}\t";
          print CLUST "$all_clusters{$spp}{$clusternumber}{'clust_dist'}\t";
          print CLUST "$gene_dist\t";
          print CLUST "@{$all_clusters{$spp}{$clusternumber}{'all_gene_ids'}}\n";
          print "@{$all_clusters{$spp}{$clusternumber}{'all_gene_ids'}}\n";
          $clusternumber++;
        }
        # reset for next cluster
        $clustsign = 10;
        undef %cluster;
        undef $clust_start;
        $initiate = 0;
        $clust_flag = 0;
        $gapflag = 0;
        $gene_dist = 0;
        $last_clust_end = $raw{$spp}{$clust_contig}{$start_pos}{'end'};
      }
      # note (continue to evaluate from the same position in case there are more than
      # one cluster beside one another - e.g. if the direction of DE changes)


      # flag if there is a potential cluster start (clust > 0 and clustsign == 10)
      if ( ($clustsign == 10) && ($raw{$spp}{$contig}{$start_pos}{'clust'} > 0) )
      {
        print "initiating cluster...";
        $clust_flag = 1;
        $cluster{$start_pos} = $raw{$spp}{$contig}{$start_pos}{'clust'};
        $clustsign = $raw{$spp}{$contig}{$start_pos}{'posDE'};
        $clust_start = $start_pos;
        $clust_contig = $contig;
      }
      # flag if there is a potential cluster extention (clust > 0; gapflag <= threshold;  clustsign == posDE)
      elsif ( ($clustsign == 0) || ($clustsign == 1) && ($gapflag <= $gaps_allowed) && ($raw{$spp}{$contig}{$start_pos}{'distprev'} < $distthresh) && ($clustsign == $raw{$spp}{$contig}{$start_pos}{'posDE'}) )
      {
        print "extending cluster...";
        $cluster{$start_pos} = $raw{$spp}{$contig}{$start_pos}{'clust'};
      }

      # reset for next gene
      $prev = $start_pos;
      unless ( ($clust_flag == 1) || ($prev == 0) ) {$gene_dist++;}
      print "\n";
    }
  }
}
