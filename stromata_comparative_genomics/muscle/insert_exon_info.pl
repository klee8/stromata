# get exon info from gff3 file
# get positions of gaps from aligned elymiNfE728 sequence from aln.fa file (first entry in each file)
# create a line to add to the alignment file which annotates where the exons are

#!usr/bin/perl -w
use strict;

my $exonpos = $ARGV[0];
my $gaps = $ARGV[1];
my $length = $ARGV[2];
my $orientation = $ARGV[3];
my $gene = $ARGV[4];

print "SeqLength: $length\nOrientation: $orientation\n";
# get lag in gene position
my @gene_pos = split(":", $gene);
my $beg = shift @gene_pos;
my $end = shift @gene_pos;
my $lag;
print "gene: $beg-$end, lag: $lag\n";

# read in gene positions and exon start and end to array
my @exons;
open(EXONS, "<$exonpos") || die "ERROR, couldn't open $exonpos: $!";
my $line = "";
while(<EXONS>) {
    print "exons: $_\n";
    chomp;
    if ($orientation eq "+") { $line = "$line $_"; }
    else { 
        my ($exonBeg, $exonEnd) = split( " ", $_);
        my $newexonB = $end - ($exonEnd - $beg);
        my $newexonE = $end - ($exonBeg - $beg);
        $_ = "$newexonB $newexonE";
        $line = "$line $_"; 
        print "$exonBeg, $exonEnd => $newexonB, $newexonE\n";
    }     
}
my @exons = split(" ", $line);


# lag compensates for seq before start of gene and perl counting from zero
# this check is not needed - gff file always seems to have smaller position first
if ($beg < $end) { $lag = $beg;}
else { $lag = $end;}
print "beg: $beg\nend: $end\nlag: $lag\n";
#print "all exons: @exons\n";

# read in gap positions to array
my @gap_pos;
open (GAPS, "<$gaps") || die "ERROR, couldn't open $gaps: $!";
while(<GAPS>){
    chomp;
    push (@gap_pos, $_);
}
print "all gaps: @gap_pos\n";

# shift
my $nextgap = shift @gap_pos;
my $annotation = "";
my $gaplength = 0;
my $exonflag = 0;
my $exonstart = shift @exons;
my $exonend = shift @exons;
my $i = 1; 
print "exon$i: $exonstart $exonend \n";
for (my $pos=$lag; $pos <= $length + $lag -1; $pos++) {
    # add in alignment gaps (gaps counted from zero)
    if ($pos == $nextgap + $lag) { 
        $annotation = $annotation."-"; 
        $nextgap = shift @gap_pos;
        $gaplength++;
        next;
    }
    # if the nucleotide is in an exon, annotate as 'C' (Dan's macVector step can't handle non nucl)
    elsif ( (($exonstart + $gaplength)<= $pos) && ($pos <= ($exonend + $gaplength)) ) {
        $annotation = $annotation."E";
    }
    # else annotate as 'N'
    else { $annotation = $annotation."N"; }
    # if you are at the last position of an exon, get the next exon info
    if ($pos == ($exonend + $gaplength)) {
        $exonstart = shift @exons;
        $exonend = shift @exons;
        $i++;
        unless ( $exonstart == "" ) {print "exon$i: $exonstart $exonend \n";}
    }    
}

open(OUT, ">temp") || die "ERROR: can't open temp file:$!"; 
print OUT ">EelyNFe728 exon positions\n";
print OUT "$annotation\n";
