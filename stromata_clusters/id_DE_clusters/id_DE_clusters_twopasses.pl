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
    if ($l2fc eq "NA") { next;}
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
  foreach my $contig (sort keys %{$raw{$spp}}){
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
    }
  }
}


my $clustsign = 10;   # used as a cluster flag and to track cluster sign
my %check_cluster;
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
print CLUST "species\tcluster\tcontig\tstart\tend\tdist_to_last_cluster\tgenes_between_clusters\tgenes\tend_of_contig\n";
print "# identifying clusters...\n";
print "contig\tgene_id\tgene_start\gene_end\tcluster_eval\tDE_direction\t";
print "cluster_sign\tgapflag\tbp_dist_to_next_cluster\tdistance_threshold\tgene_count_to_next_cluster\tcluster_contig\n";

# identify initial genes in clusters
foreach my $spp (sort keys %raw){
  $clusternumber = 1;
  foreach my $contig (sort keys %{$raw{$spp}}){
    my $prev = 0; $gapflag = 0;  $prev_cluster_end = 0; unless ($clustsign != 10 ) {$gene_dist = 0;}
    foreach my $start_pos (sort { $a <=> $b } keys %{$raw{$spp}{$contig}}){

      # check for gaps
      if ( ($clustsign != 10 ) && ($raw{$spp}{$contig}{$start_pos}{'clust'} == 0) ) { $gapflag++;}

      # flag genes that coud initiate a cluster
      if ($raw{$spp}{$contig}{$start_pos}{'clust'} == 2)
      {
        $initiate = 1;
      }

      # temp test print checks
      unless ($start_pos == 0) {
        print "$contig\t$raw{$spp}{$contig}{$start_pos}{'gene_id'}\t";
        print "$start_pos\t$raw{$spp}{$contig}{$start_pos}{'end'}\t";
        print "$raw{$spp}{$contig}{$start_pos}{'clust'}\t$raw{$spp}{$contig}{$start_pos}{'posDE'}\t$clustsign\t$gapflag\t";
        print "$raw{$spp}{$contig}{$start_pos}{'distprev'}\t$distthresh\t$gene_dist\t$clust_contig\t";
      }

      # if you are in a cluster and pass more gaps than are allowed or
      # the distance to previous gene is larger than the thresholds
      # or the direction of DE sign changes
      # or you get to the next contig
      # evaluate potential cluster
      if ( ($clustsign != 10) && ( ($gapflag > $gaps_allowed) || ($raw{$spp}{$contig}{$start_pos}{'distprev'} > $distthresh) || ($clustsign != $raw{$spp}{$contig}{$start_pos}{'posDE'}) || ($contig ne $clust_contig) ) )
      {
        print "evaluate cluster...";


        # TRIM TRAILING GAPS
        foreach my $clust_gene_start (sort { $b <=> $a } keys %check_cluster )
        {
          # if the last gene of a cluster is a 'gap'
          if ($raw{$spp}{$clust_contig}{$clust_gene_start}{'clust'} == 0)
          {
            $trimmed++;
            my $prev_start = $raw{$spp}{$clust_contig}{$clust_gene_start}{'prev'};
            $cluster_end = $raw{$spp}{$clust_contig}{$prev_start}{'end'};
#             print "\nTRIMMING: previous start = $prev_start\tcluster end = $cluster_end\tprevious cluster end = $prev_cluster_end\n";
            delete $check_cluster{$clust_gene_start};
          }
          else
          {
              $cluster_end = $raw{$spp}{$clust_contig}{$clust_gene_start}{'end'};
              last;
          }
        }

        # if there is a gene that meets the 'intiate cluster thresholds'
        # and there are enough genes for a cluster
        if ( ($initiate == 1) && ($mingenenum <= keys %check_cluster) )
        {
          # remove cluster genes from gene-distance calculation
          my $cluster_size = (keys %check_cluster);
          if ($contig = $clust_contig) { $same_contig = 1;}
#          print "\ngene distance : $gene_dist = $gene_dist - $cluster_size - $same_contig - $trimmed\n";
          $gene_dist = $gene_dist - $cluster_size - $same_contig - $trimmed;
          $same_contig = 0;
#          print "\n$spp\t$clusternumber:\t";
          foreach my $start_pos (sort { $a <=> $b } keys %check_cluster )
          {
            push(@{$all_clusters{$spp}{$clusternumber}{'all_gene_ids'}}, $raw{$spp}{$clust_contig}{$start_pos}{'gene_id'});
            $all_clusters{$spp}{$clusternumber}{'gene_start'}{$start_pos}{'gene_id'} = $raw{$spp}{$clust_contig}{$start_pos}{'gene_id'};

            # assign clusternumbers to orthologs
            my $orthogroup = $raw{$spp}{$clust_contig}{$start_pos}{'orthogroup'};
            $orthogroups{$orthogroup}{$spp} = $clusternumber;
            $orthogroup{'cluster'}{$spp}{$clusternumber}{$orthogroup} = $all_clusters{$spp}{$clusternumber}{'gene_start'}{$start_pos}{'gene_id'};
            print "orthogroup:$orthogroup;cluster:$clusternumber\t";
          }
          $all_clusters{$spp}{$clusternumber}{'clust_dist'} = $clust_start - $prev_cluster_end;
#          print "calculating cluster distance $clust_start - $prev_cluster_end =  $all_clusters{$spp}{$clusternumber}{'clust_dist'}\n";
          $all_clusters{$spp}{$clusternumber}{'gene_dist'} = $gene_dist;
#          print "\ngene distance is $gene_dist\n";
          $all_clusters{$spp}{$clusternumber}{'contig'} = $clust_contig;
          $all_clusters{$spp}{$clusternumber}{'end'} = $cluster_end;
          $prev_cluster_end = $cluster_end;
          # cluster at start of contig
          if ($gene_dist == 0) { $all_clusters{$spp}{$clusternumber}{'start_of_contig'} = "TRUE\t";}    #
          else { $all_clusters{$spp}{$clusternumber}{'start_of_contig'} =  "FALSE\t";}
          #cluster at end of contig
          if ($contig ne $clust_contig) { $all_clusters{$spp}{$clusternumber}{'end_of_contig'} =  "TRUE";}
          else {  $all_clusters{$spp}{$clusternumber}{'end_of_contig'} =  "FALSE";}
          # size of cluster
          $all_clusters{$spp}{$clusternumber}{'cluster_size'} = keys %{$all_clusters{$spp}{$clusternumber}};

          # species, cluster number and position info
          print CLUST "$spp\t$clusternumber\t$clust_contig\t$clust_start\t$cluster_end\t";
          # distance to previous cluster (bn and number of genes between them)
          print CLUST "$all_clusters{$spp}{$clusternumber}{'clust_dist'}\t$gene_dist\t";
          print CLUST "@{$all_clusters{$spp}{$clusternumber}{'all_gene_ids'}}\t";
          # cluster at begining or end of contig
          print CLUST "$all_clusters{$spp}{$clusternumber}{'start_of_contig'}\t$all_clusters{$spp}{$clusternumber}{'end_of_contig'}";
          print CLUST "\n";

          print "CLUSTER!!!!";
          #print "\n\n###\t$spp\t$clusternumber\t$clust_contig\t$clust_start\t$cluster_end\t";
          #print "$all_clusters{$spp}{$clusternumber}{'clust_dist'}\t$gene_dist\t";
          #print "@{$all_clusters{$spp}{$clusternumber}{'all_gene_ids'}}\n\n";
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
        print "initiating cluster...";
        $check_cluster{$start_pos} = $raw{$spp}{$contig}{$start_pos}{'clust'};
        $clustsign = $raw{$spp}{$contig}{$start_pos}{'posDE'};
        $clust_start = $start_pos;
        $clust_contig = $contig;
      }
      # flag if there is a potential cluster extention (clust > 0; gapflag <= threshold;  clustsign == posDE)
      elsif ( ($clustsign == 0) || ($clustsign == 1) && ($gapflag <= $gaps_allowed) && ($raw{$spp}{$contig}{$start_pos}{'distprev'} < $distthresh) && ($clustsign == $raw{$spp}{$contig}{$start_pos}{'posDE'}) )
      {
        print "extending cluster...";
        $check_cluster{$start_pos} = $raw{$spp}{$contig}{$start_pos}{'clust'};
      }

      # reset for next gene
      $prev = $start_pos;
      $gene_dist++;
      unless ($start_pos == 0) { print "\n"; }
    }
  }
}

##################################################################################################
########## first pass of orthologs
#   $orthogroups{$orthogroup}{$spp} = $clusternumber;          <<<<<<<<<<<<<<<<<<<<<<<<<<<<< orthogroup filled if cluster requirements fulfilled
#   $all_clusters{$spp}{$clusternumber}{'gene_start'}{$start_pos}{'gene_id'} = $raw{$spp}{$contig}{$start_pos}{'gene_id'};
#   foreach my $orthogroup (sort $a <=> $b keys %{$orthogroups{$orthogroups}}){
my %output;
my $orthocluster;
my $group = 1;



foreach my $orthogroup (sort {$a <=> $b} keys %orthogroups) {

  # check the number of species represented
  my $ortho_size = keys %{$orthogroups{$orthogroup}};
  $output{$group}{'num_species'} = $ortho_size;
  print "$orthogroup\t#species:$ortho_size\t";

  # check the maximum size of the culster for the orthogroup in all species
  my $max_cluster_size = 0;
  my $max_clust_spp = "";
  foreach my $species (sort {$a <=> $b} keys %{$orthogroups{$orthogroup}}) {
    my $clusternumber = $orthogroups{$orthogroup}{$species};
    my $cluster_size = keys %{$all_clusters{$species}{$clusternumber}{'gene_start'}};
    $output{$group}{'species'}{$species} = $cluster_size;
    print "$species:$cluster_size; ";
    if ($cluster_size > $max_cluster_size) {
      $max_cluster_size = $cluster_size;
      $max_clust_spp = $species;
    }
  }
  $output{$group}{'max_cluster_size'} = $max_cluster_size;
  $output{$group}{'max_cluster_species'} = $max_cluster_species;

  # check all other orthologs in cluster how they match other spp
  foreach my $species (sort {$a <=> $b} keys %{$orthogroups{$orthogroup}}) {

    my $current_cluster = $orthogroups{$orthogroup}{$species};

    foreach my $ortholog (sort {$a <=> $b} keys %{$orthogroup{'cluster'}{$spp}{$current_cluster}}) {


      # if the number of species represented by an ortholog is > 1
      # get orthologs of each gene in the max cluster and check they are in the same cluster in other spp


    }

  }
  print "\n";
}
          $orthogroup{'cluster'}{$spp}{$clusternumber}{$orthogroup} = $all_clusters{$spp}{$clusternumber}{'gene_start'}{$start_pos}{'gene_id'};
