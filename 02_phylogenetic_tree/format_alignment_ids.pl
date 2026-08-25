#!/usr/bin/perl
use strict;
use warnings;

my ($in,$out)=@ARGV;

if(@ARGV != 2){
    die "Usage: perl phylip_relaxed_to_fasta.pl input.phy output.fa\n";
}

open IN, "<", $in or die "Cannot open $in\n";
open OUT, ">", $out or die "Cannot write $out\n";


# First line PHYLIP header
my $header=<IN>;
chomp $header;

my ($ntax,$nsite)=split(/\s+/, $header);

print STDERR "Taxa: $ntax\n";
print STDERR "Sites: $nsite\n";


my $count=0;

while(<IN>){

    chomp;

    next if /^\s*$/;


    # relaxed PHYLIP:
    # ID      sequence

    if(/^(\S+)\s+(\S+)$/){

        my $id=$1;
        my $seq=$2;

        print OUT ">$id\n";
        print OUT "$seq\n";

        $count++;

    }else{

        warn "Skip line: $_\n";

    }

}


close IN;
close OUT;


print STDERR "Converted sequences: $count\n";

if($count != $ntax){
    print STDERR "WARNING: expected $ntax sequences but got $count\n";
}
