# id_DE_clusters.pl Kate D Lee July 2019
# usef file with columns for:
#  species, chr, start, end, orthogroup, log2fc, log2fcSE, svalue
# identify clusters in each species with DE info and positional info
# use ortholog info to compare clusters across species
#############################################################

#!usr/bin/perl -w
use strict;

my $file = $ARGV[0] || "sub_core_gene_set_STR_PS_ann.txt";
open(IN, "<$file") || die "ERROR: couldn't open file $file: $!";

my $minlog2fc = 2;
my $maxsvalue = 0.005;
my $minlog2fcEXTN = 1;
my $maxsvalueEXTN = 0.01;
my $gaps_allowed = 1;
my $distthresh = 50000;
my $mingenenum = 4;
my %raw;
my %orths;

open(LOG, ">log_cluster_finding.txt") || die "ERROR, couldn't open log_cluster_finding.txt: $!";

print LOG "maxsvalue: $maxsvalue\nminlog2fc: $minlog2fc\n";
print LOG "maxsvalueEXTN: $maxsvalueEXTN\nminlog2fcEXTN: $minlog2fcEXTN\n";
print LOG "distance threshold: $distthresh\n";
print LOG "reading in data...\n";

my $prev;
# read in data to hash and calculate characteristics of gene
while(<IN>) {
  chomp;
  if ($_=~ m/orthogroup/) {
      next;
  }
  else {
    my ($contig, $start, $end, $strand, $gene_id, $orthogroup, $l2fc, $l2fcSE, $svalue_1, $svalue_2, $spp) = split("\t", $_);
    #if ($l2fc eq "NA") { next;}
    $raw{$spp}{$contig}{$start}{'end'} = $end;
    $raw{$spp}{$contig}{$start}{'strand'} = $strand;
    $raw{$spp}{$contig}{$start}{'gene_id'} = $gene_id;
    $raw{$spp}{$contig}{$start}{'orthogroup'} = $orthogroup;
    $raw{$spp}{$contig}{$start}{'l2fc'} = $l2fc;
    $raw{$spp}{$contig}{$start}{'l2fcSE'} = $l2fcSE;
    $raw{$spp}{$contig}{$start}{'svalue_1'} = $svalue_1;
    $raw{$spp}{$contig}{$start}{'svalue_2'} = $svalue_2;
    #print "spp: $spp\tsvalue_2: $svalue_2\tgene_id:$raw{$spp}{$contig}{$start}{'gene_id'}\n";                #:)
    $orths{$spp}{$orthogroup} = "$contig:$start";
  }
}

my $prev;

print LOG "calculating gene characteristics...\n";
# first pass to calculate characteristics of genes (cluster/extension thresholds; distances to next gene etc.)
foreach my $spp (sort keys %raw){
  foreach my $contig (sort keys %{$raw{$spp}}){
    my $prev = 0; my $svalue_1; my $svalue_2;
    foreach my $start_pos (sort { $a <=> $b } keys %{$raw{$spp}{$contig}}){
      $raw{$spp}{$contig}{$start_pos}{'prev'} = $prev;
      my $gene_id = $raw{$spp}{$contig}{$start_pos}{'gene_id'};
      my $l2fc = $raw{$spp}{$contig}{$start_pos}{'l2fc'};
      my $svalue_1 = $raw{$spp}{$contig}{$start_pos}{'svalue_1'};
      my $svalue_2 = $raw{$spp}{$contig}{$start_pos}{'svalue_2'};

      # check for missing data
      if (($l2fc eq 'NA') || ($l2fc eq '')) {$raw{$spp}{$contig}{$start_pos}{'clust'} = 0;}

      # check whether to initiate cluster
      elsif ( ($svalue_1 < $maxsvalue) && ( ( ($l2fc  > 0) && ($l2fc > $minlog2fc) ) || ( ($l2fc  < 0) && ($l2fc  < -$minlog2fc) ) ) )
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
    }
  }
}

my $clustsign = 10;   # used as a cluster flag and to track cluster sign
my %check_cluster;
my %count_cluster;
my $initiate = 0;
my $clusternumber = 1;
my $cluster_end;
my $gapflag = 0;
my %all_clusters;
my $clust_contig;
my $clust_start;
my $prev_cluster_end;
my $gene_dist = 0;
my $trimmed = 0;
my %orthogroups;
my $same_contig = 1;

open(CLUST, ">clusters.txt") || die "ERROR: couldn't open clusters.txt outfile: $!";
print CLUST "species\tcluster\tcontig\tcluster_size\tstart\tend\tdist_to_last_cluster\tgenes_between_clusters\tgenes\tortholog_numbers\tstart_of_contig\tend_of_contig\n";

print LOG "# identifying clusters...\n";
print LOG "contig\tgene_id\tgene_start\gene_end\tcluster_eval\tDE_direction\t";
print LOG "cluster_sign\tgapflag\tbp_dist_to_next_cluster\tdistance_threshold\tgene_count_to_next_cluster\tcluster_contig\n";

# identify initial genes in clusters
foreach my $spp (sort keys %raw){
  $clusternumber = 1;
  foreach my $contig (sort keys %{$raw{$spp}}){
    my $prev = 0; $gapflag = 0;  $prev_cluster_end = 0; unless ($clustsign != 10 ) {$gene_dist = 0;}
    foreach my $start_pos (sort { $a <=> $b } keys %{$raw{$spp}{$contig}}){

      # check for gaps
      if ( ($clustsign != 10 ) &&  ($raw{$spp}{$contig}{$start_pos}{'clust'} == 0) ) { $gapflag++;}

      # flag genes that coud initiate a cluster
      if ($raw{$spp}{$contig}{$start_pos}{'clust'} == 2)
      {
        $initiate = 1;
      }

      # temp test print checks
      unless ($start_pos == 0) {
        print LOG "$contig\t$raw{$spp}{$contig}{$start_pos}{'gene_id'}\t";
        print LOG "$start_pos\t$raw{$spp}{$contig}{$start_pos}{'end'}\t";
        print LOG "$raw{$spp}{$contig}{$start_pos}{'clust'}\t$raw{$spp}{$contig}{$start_pos}{'posDE'}\t$clustsign\t$gapflag\t";
        print LOG "$raw{$spp}{$contig}{$start_pos}{'distprev'}\t$distthresh\t$gene_dist\t$clust_contig\t";
      }

      # if you are in a cluster and pass more gaps than are allowed or
      # the distance to previous gene is larger than the thresholds
      # or the direction of DE sign changes
      # or you get to the next contig
      # evaluate potential cluster
      if ( ($clustsign != 10) && ( ($gapflag > $gaps_allowed) || ($raw{$spp}{$contig}{$start_pos}{'distprev'} > $distthresh) || ($clustsign != $raw{$spp}{$contig}{$start_pos}{'posDE'}) || ($contig ne $clust_contig) ) )
      {
        print LOG "evaluate cluster...\n";
        if ($raw{$spp}{$contig}{$start_pos}{'clust'} == 0) { $gapflag = $gapflag - 1;}

        # TRIM TRAILING GAPS
        foreach my $clust_gene_start (sort { $b <=> $a } keys %check_cluster )
        {
          # if the last gene of a cluster is a 'gap'
          if ($raw{$spp}{$clust_contig}{$clust_gene_start}{'clust'} == 0)
          {
            $trimmed++;
            my $prev_start = $raw{$spp}{$clust_contig}{$clust_gene_start}{'prev'};
            $cluster_end = $raw{$spp}{$clust_contig}{$prev_start}{'end'};
            print LOG "\nTRIMMING: previous start = $prev_start\tcluster end = $cluster_end\tprevious cluster end = $prev_cluster_end\n";
            $gapflag = $gapflag - 1;
            delete $check_cluster{$clust_gene_start};
          }
          else
          {
              $cluster_end = $raw{$spp}{$clust_contig}{$clust_gene_start}{'end'};
              my $length_cluster = keys %check_cluster;
              print LOG "\ncluster size = $length_cluster, gapflag = $gapflag\n";
              last;
          }
        }

        # if there is a gene that meets the 'intiate cluster thresholds'
        # and there are enough genes for a cluster (i.e. minimum gene number, or less than the minimum gene number but with fewer gaps than allowed)
        if ( ($initiate == 1) && ( ($mingenenum <= keys %check_cluster)  ||  ( (($mingenenum -1) <= keys %check_cluster) && ($gapflag <= ($gaps_allowed -1) ) ) ) )
        #if ( ($initiate == 1) && ($mingenenum <= keys %check_cluster) )
        {
          print LOG "cluster initiated\n";
          # remove cluster genes from gene-distance calculation
          my $cluster_size = (keys %check_cluster);
          if ($contig = $clust_contig) { $same_contig = 1;}
          $gene_dist = $gene_dist - $cluster_size - $same_contig - $trimmed;
          $same_contig = 0;

          foreach my $start_pos (sort { $a <=> $b } keys %check_cluster )
          {
            push(@{$all_clusters{$spp}{$clusternumber}{'all_gene_ids'}}, $raw{$spp}{$clust_contig}{$start_pos}{'gene_id'});
            push(@{$all_clusters{$spp}{$clusternumber}{'all_orth_ids'}}, $raw{$spp}{$clust_contig}{$start_pos}{'orthogroup'});
            $all_clusters{$spp}{$clusternumber}{'gene_start'}{$start_pos}{'gene_id'} = $raw{$spp}{$clust_contig}{$start_pos}{'gene_id'};

            # assign clusternumbers to orthologs
            my $orthogroup = $raw{$spp}{$clust_contig}{$start_pos}{'orthogroup'};
            $orthogroups{$orthogroup}{$spp}{'clustnum'} = $clusternumber;
            $orthogroups{$orthogroup}{$spp}{'contig'} = $clust_contig;
            $orthogroups{'cluster'}{$spp}{$clusternumber}{$orthogroup} = $all_clusters{$spp}{$clusternumber}{'gene_start'}{$start_pos}{'gene_id'};
            $orthogroups{'start'}{$spp}{$clusternumber}{$orthogroup} = $start_pos;
            $orthogroups{'end'}{$spp}{$clusternumber}{$orthogroup} = $raw{$spp}{$clust_contig}{$start_pos}{'end'};

            # initiate a placeholder to loop through later
            $count_cluster{$spp}{$clusternumber}{$orthogroup} = 1;
          }
          $all_clusters{$spp}{$clusternumber}{'clust_dist'} = $clust_start - $prev_cluster_end;
          $all_clusters{$spp}{$clusternumber}{'gene_dist'} = $gene_dist;
          $all_clusters{$spp}{$clusternumber}{'contig'} = $clust_contig;
          $all_clusters{$spp}{$clusternumber}{'clust_end'} = $cluster_end;
          $prev_cluster_end = $cluster_end;
          # cluster at start of contig
          if ($gene_dist == 0) { $all_clusters{$spp}{$clusternumber}{'start_of_contig'} = "TRUE";}    #
          else { $all_clusters{$spp}{$clusternumber}{'start_of_contig'} =  "FALSE";}
          #cluster at end of contig
          if ($contig ne $clust_contig) { $all_clusters{$spp}{$clusternumber}{'end_of_contig'} =  "TRUE";}
          else {  $all_clusters{$spp}{$clusternumber}{'end_of_contig'} =  "FALSE";}
          # size of cluster
          $all_clusters{$spp}{$clusternumber}{'cluster_size'} = keys %{$all_clusters{$spp}{$clusternumber}};
          $all_clusters{$spp}{$clusternumber}{'clust_start'} = $clust_start;

          # species, cluster number and position info
          print CLUST "$spp\t$clusternumber\t$clust_contig\t$cluster_size\t$clust_start\t$cluster_end\t";
          # distance to previous cluster (bn and number of genes between them)
          print CLUST "$all_clusters{$spp}{$clusternumber}{'clust_dist'}\t$gene_dist\t";
          print CLUST join(',', @{$all_clusters{$spp}{$clusternumber}{'all_gene_ids'}}), "\t";
          print CLUST join(',', @{$all_clusters{$spp}{$clusternumber}{'all_orth_ids'}}), "\t";
          # cluster at begining or end of contig
          print CLUST "$all_clusters{$spp}{$clusternumber}{'start_of_contig'}\t$all_clusters{$spp}{$clusternumber}{'end_of_contig'}";
          print CLUST "\n";

          print LOG "CLUSTER!!!!";
          print LOG "\t$spp\t$clusternumber\t$clust_contig\t$clust_start\t$cluster_end\t";
          print LOG "$all_clusters{$spp}{$clusternumber}{'clust_dist'}\t$gene_dist\t";
          print LOG join(',', @{$all_clusters{$spp}{$clusternumber}{'all_gene_ids'}}), "\t";
          print LOG join(',', @{$all_clusters{$spp}{$clusternumber}{'all_orth_ids'}}), "\n\n";
          $clusternumber++;
          $gene_dist = 0;
        }

        # reset for next cluster or potential cluster
        $clustsign = 10;
        undef %check_cluster;
        undef $clust_start;
        $initiate = 0;
        $gapflag = 0;
        $trimmed = 0;
      }
      # note (continue to evaluate from the same position in case there are more than
      # one cluster beside one another - e.g. if the direction of DE changes)


      # flag if there is a potential cluster start (clust > 0 and clustsign == 10)
      if ( ($clustsign == 10) && ($raw{$spp}{$contig}{$start_pos}{'clust'} > 0) )
      {
        print LOG "initiating cluster...";
        $check_cluster{$start_pos} = $raw{$spp}{$contig}{$start_pos}{'clust'};
        $clustsign = $raw{$spp}{$contig}{$start_pos}{'posDE'};
        $clust_start = $start_pos;
        $clust_contig = $contig;
      }
      # flag if there is a potential cluster extention (clust > 0; gapflag <= threshold;  clustsign == posDE)
      elsif ( ($clustsign == 0) || ($clustsign == 1) && ($gapflag <= $gaps_allowed) && ($raw{$spp}{$contig}{$start_pos}{'distprev'} < $distthresh) && ($clustsign == $raw{$spp}{$contig}{$start_pos}{'posDE'}) )
      {
        print LOG "extending cluster...";
        $check_cluster{$start_pos} = $raw{$spp}{$contig}{$start_pos}{'clust'};
      }

      # reset for next gene
      $prev = $start_pos;
      $gene_dist++;
      unless ($start_pos == 0) { print LOG "\n"; }
    }
  }
}

##################################################################################################
##########          CHECK FOR OVERLAPPING CLUSTERS and collapse into $out_groups          ########
##################################################################################################

my %out_groups;
my %check;
my $group = 1;
my %overlaps;


# iterate through species, clusters and their orthologs
foreach my $species (sort keys %all_clusters) {
  print LOG "\nchecking clusters in $species....\n";
  foreach my $clusternumber (sort {$a <=> $b} keys %{$count_cluster{$species}}) {
    print LOG "PASS_1\tgroup:$group\t$species\tcluster:$clusternumber\torthologs:";

    # FIRST PASS - check orthologs in original cluster
    foreach my $orth (@{$all_clusters{$species}{$clusternumber}{'all_orth_ids'}}){
      print LOG "$orth,";
      $out_groups{$group}{$species}{$clusternumber}{$orth} = 1;
      foreach my $spp (sort {$a <=> $b} keys %all_clusters) {
        unless ($spp eq $species) {
          if (exists $orthogroups{$orth}{$spp}{'clustnum'}) {
            my $clustnum = $orthogroups{$orth}{$spp}{'clustnum'};
            $out_groups{$group}{$spp}{$clustnum}{$orth} = 1;
            $overlaps{$group}{$spp}{$clustnum}{$orth} = 1;
            $overlaps{$group}{$species}{$clusternumber}{$orth} = 1;
            $check{$group}{'PASS1'}{$spp}{$clustnum} = 1;
          }
        }
      }
    }
    # print out other species and their clusters where orthologs of the initial cluster were found
    foreach my $spp (sort keys %{$check{$group}{'PASS1'}}){
      foreach my $cluster (sort keys %{$check{$group}{'PASS1'}{$spp}}){
        print LOG "\t$spp:$cluster\t";
      }
    }
    print LOG "\n";

    # SECOND PASS
    # for any clusters found in other species in the first pass, check any orthologs not linked to the first cluster
    # identify any new clusters linked
    foreach my $spp (sort {$a <=> $b} keys %{$check{$group}{'PASS1'}}) {
      foreach my $clust (sort keys %{$check{$group}{'PASS1'}{$spp}}) {
        print LOG "PASS_2\t$spp\tcluster:$clust\textra_orthologs:";
        foreach my $orth (@{$all_clusters{$spp}{$clust}{'all_orth_ids'}}){
          unless (exists($out_groups{$group}{$spp}{$clust}{$orth})) {
            print LOG "$orth,";
            $out_groups{$group}{$spp}{$clust}{$orth} = 1;
            # check extra orthologs and if they are linked to new clusters in other species
            foreach my $checkspp (sort keys %all_clusters) {
              unless ($checkspp eq $species) {
                my $clustnum = $orthogroups{$orth}{$checkspp}{'clustnum'};
                unless ( exists($check{$group}{'PASS1'}{$checkspp}{$clustnum}) ) {
                  if (exists($orthogroups{$orth}{$checkspp}{'clustnum'})) {
                    print LOG "\t***$checkspp:$clustnum***\t";
                    $out_groups{$group}{$checkspp}{$clustnum}{$orth} = 1;
                    $overlaps{$group}{$checkspp}{$clustnum}{$orth} = 1;
                    $overlaps{$group}{$species}{$clusternumber}{$orth} = 1;
                    $check{$group}{'PASS2'}{$checkspp}{$clustnum} = 1;
                  }
                }
              }
            }
          }
        }
        # print out other species and their clusters where orthologs of the initial cluster were found
        foreach my $spp (sort keys %{$check{$group}{'PASS2'}}){
          foreach my $cluster (sort keys %{$check{$group}{'PASS2'}{$spp}}){
            print LOG "\t$spp:$cluster\t";
          }
        }
        print "\n";
      }
    }


    # THRID pass
    # on the off-chance that a cluster with no orthologs in common with the original is found, check for any other links
    foreach my $spp (sort {$a <=> $b} keys %{$check{$group}{'PASS2'}}) {
      foreach my $clust (sort keys %{$check{$group}{'PASS2'}{$spp}}){
        print LOG "PASS_3\t$spp\tcluster:$clust\textra_orthologs:";
        foreach my $orth (sort keys %{$count_cluster{$spp}{$clust}}) {
          unless (exists($out_groups{$group}{$spp}{$clust}{$orth})) {
            print LOG "$orth,";
            $out_groups{$group}{$spp}{$clust}{$orth} = 1;
            # check extra orthologs and if they are linked to new clusters in other species
            foreach my $checkspp (sort {$a <=> $b} keys %{$out_groups{$group}}) {
              my $clustnum = $orthogroups{$orth}{$checkspp}{'clustnum'};
              unless ( exists($check{$group}{'PASS1'}{$checkspp}{$clustnum}) || exists($check{$group}{'PASS2'}{$checkspp}{$clustnum})) {
                if (exists($orthogroups{$orth}{$checkspp}{'clustnum'})) {
                  print LOG "\t***$checkspp:$clustnum***\t";
                  $out_groups{$group}{$checkspp}{$clustnum}{$orth} = 1;
                  $overlaps{$group}{$checkspp}{$clustnum}{$orth} = 1;
                  $check{$group}{'PASS3'}{$checkspp}{$clustnum} = 1;
                }
              }
            }
          }
        }
        print LOG "\n";
      }
    }

    # remove original $clusternumber from $all_clusters so you don't put it into another group
    delete $count_cluster{$species}{$clusternumber};

    # remove any linked clusters found in PASS1
    foreach my $spp (sort {$a <=> $b} keys %{$out_groups{$group}}) {
      foreach my $clust (sort keys %{$check{$group}{'PASS1'}{$spp}}) {
        delete $count_cluster{$spp}{$clust};
      }
    }
    # remove any linked clusters found in PASS2
    foreach my $spp (sort {$a <=> $b} keys %{$out_groups{$group}}) {
      foreach my $clust (sort keys %{$check{$group}{'PASS2'}{$spp}}){
        delete $count_cluster{$spp}{$clust};
      }
    }
    # remove any linked clusters found in PASS3
    foreach my $spp (sort {$a <=> $b} keys %{$out_groups{$group}}) {
      foreach my $clust (sort keys %{$check{$group}{'PASS3'}{$spp}}){
        delete $count_cluster{$spp}{$clust};
      }
    }
    $group++;
  }
}

################################################################################
#####        Iterate through groups and write graph scripts                #####
################################################################################

# flat file of groups with one row for each ortholog for each species
open(GROUPS, ">cluster_groups.txt") || die "ERROR: couldn't open cluster_groups.txt: $!";
print GROUPS "group\tnum_spp\tspecies\tcluster\tcluster_size\tgene_id\tortholog\tstart\tend\tstrand\n";

# print out list of clusters that are present in more than one species
foreach my $group (sort {$a <=> $b} keys %out_groups) {
  my $species; my $cluster; my $clust_len;
  my $clust_start;
  my $clust_end;
  my $num_spp = keys %{$out_groups{$group}};
#  if ($num_spp > 1) {
    foreach my $spp (sort keys %{$out_groups{$group}}) {
      $species = $spp;
      # open a graph file
      open(GRAPH, ">graphs/graph.$group.$spp.py") || die "ERROR: couldn't open graph_$group.txt: $!";

      my $note = <<'END_LINE';

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from reportlab.lib.units import cm

gdd = GenomeDiagram.Diagram('graph')
gdt_features = gdd.new_track(1, greytrack=False)
gds_features = gdt_features.new_set()

# add features
END_LINE

      print GRAPH $note;
      foreach my $clust (sort keys %{$out_groups{$group}{$spp}}) {

        my $clust_size = @{$all_clusters{$spp}{$clust}{'all_orth_ids'}};

        my $counter = 0; my $start; my $end; my $gap;
        push (@{$all_clusters{$spp}{$clust}{'all_orth_ids'}}, "END");
        foreach my $orth ( @{$all_clusters{$spp}{$clust}{'all_orth_ids'}} ) {

          if ($orth eq "END") { next; }

          # check if begining of cluster is at the begining of a contig
          if (($counter == 0) && ( $all_clusters{$spp}{$clust}{'start_of_contig'} eq "TRUE")) {
            $start = "TRUE";
          }
          else { $start = "FALSE"; }
          # check if begining of cluster is at the begining of a contig
          if (($counter eq "END") && ( $all_clusters{$spp}{$clust}{'end_of_contig'} eq "TRUE")) {
            my $temp = ${$all_clusters{$spp}{$clust}{'all_orth_ids'}}[$counter -1];
            $out_groups{$group}{$spp}{$clust}{$temp}{'end'} = "TRUE";
            $end = "TRUE";
          }
          else { $end= "FALSE"; }

          # check presence of gaps in clusters
          my $gene_id = $orthogroups{'cluster'}{$spp}{$clust}{$orth};
          my $contig = $orthogroups{$orth}{$spp}{'contig'};
          my $start = $orthogroups{'start'}{$spp}{$clust}{$orth};
          my $end = $orthogroups{'end'}{$spp}{$clust}{$orth};
          my $orientation = $raw{$spp}{$contig}{$start}{'strand'};
          my $DE_eval = $raw{$spp}{$contig}{$start}{'clust'};
          my $color = "blue";
          if ($overlaps{$group}{$spp}{$clust}{$orth} == 1)
          {
            if ($DE_eval == 0 ) { $color = "indianred"; }
            else { $color = "darkred"; }
          }
          else
          {
            if ($DE_eval == 0 ) { $color = "gainsboro"; }
            else { $color = "gray"; }
          }
          my $label_position = "start";
          my $label_angle = 0;
          if ($orientation eq "-") {
            $label_position = "end";
            $label_angle = 180;
          }

          print GROUPS "$group\t$num_spp\t$spp\t$clust\t$clust_size\t$gene_id\t$orth\t$start\t$end\t$orientation\n";
          print GRAPH "feature = SeqFeature(FeatureLocation(".$start.", ";
          print GRAPH $end."), strand=".$orientation."1)\n";
          print GRAPH "gds_features.add_feature(feature, name=\"".$orth."\", label=\"True\", color=\"".$color."\", label_size=10, label_position=\"".$label_position;
          print GRAPH "\", label_angle=".$label_angle.", sigil=\"BIGARROW\")\n";
          $counter++;
        }
        $clust_start = $all_clusters{$spp}{$clust}{'clust_start'};
        $clust_end = $all_clusters{$spp}{$clust}{'clust_end'};
        $clust_len = ($clust_end - $clust_start)/1000;
        $cluster = $clust;
      }

      if ($all_clusters{$spp}{$cluster}{'start_of_contig'} eq "TRUE")
      {
        print GRAPH "\nfeature = SeqFeature(FeatureLocation($clust_start , $clust_start + 50))\n";
        print GRAPH "gds_features.add_feature(feature, color=\"red\",\n"; # name=\"contig_break\",\n";
        print GRAPH "                         label=True, label_size = 10)\n";
        #print GRAPH "                         label_angle = 0)\n";
      }
      if ($all_clusters{$spp}{$cluster}{'end_of_contig'} eq "TRUE")
      {
        print GRAPH "\nfeature = SeqFeature(FeatureLocation($clust_end -50 , $clust_end))\n";
        print GRAPH "gds_features.add_feature(feature, color=\"red\",\n"; # name=\"contig_break\",\n";
        print GRAPH "                         label=True, label_size = 10)\n";
        #print GRAPH "                         label_angle = 0)\n";
      }

      print GRAPH "\ngdd.draw(format=\'linear\', pagesize=($clust_len*cm,4*cm), fragments=1, start=$clust_start, end=$clust_end)\n";
      print GRAPH "gdd.write(\"$group.$spp.diagram.svg\", \"SVG\")";
      close GRAPH;
    }
#  }
}
