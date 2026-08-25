#!/usr/bin/perl -w
=head1 Name
    /DATA/05.structure/clumpp.pl
Info
    Version V
    Author: st_agric (USER), user@example.com
    Created Time: 2019-01-30 16:58:23
    Created Version: clumpp.pl
Usage
    clumpp.pl
=cut
use strict;
my $cluster = shift;
my $individual_number = shift;
my $repeat_count = shift;
my $name =shift;
my $clumpp ='/TOOLS/CLUMPP_Linux64.1.1.2/CLUMPP';
open OUT,">$name\_k$cluster.paramfile";
print OUT "# --------------- Main parameters ---------------------------------------------\n\n";
print OUT "DATATYPE 0\n\n";
print OUT "INDFILE $name\_k$cluster.indfile\n\n";
print OUT "POPFILE\n\n";
print OUT "OUTFILE $name\_k$cluster.outfile\n\n";
print OUT "MISCFILE $name\_k$cluster.miscfile\n\n";
print OUT "K $cluster\n\n";
print OUT "C $individual_number\n\n";
print OUT "R $repeat_count\n\n";
print OUT "M 3\n\n";
print OUT "W 0\n\n";
print OUT "S 2\n\n";
print OUT "# - Additional options for the Greedy and LargeKGreedy algorithm (M = 2 or 3) -\n\n";
print OUT "GREEDY_OPTION 2\n\n";
print OUT "REPEATS 100\n\n";
print OUT "PERMUTATIONFILE $name\_k$cluster.permutationfile\n\n";
print OUT "# --------------- Optional outputs --------------------------------------------\n\n";
print OUT "PRINT_PERMUTED_DATA 1\n\n";
print OUT "PERMUTED_DATAFILE $name\_k$cluster.perm_datafile\n\n";
print OUT "PRINT_EVERY_PERM 0\n\n";
print OUT "EVERY_PERMFILE $name\_k$cluster.every_permfile\n\n";
print OUT "PRINT_RANDOM_INPUTORDER 0\n\n";
print OUT "RANDOM_INPUTORDERFILE $name\_k$cluster.ran_inputorderfile\n\n";
print OUT "# --------------- Advanced options --------------------------------------------\n\n";
print OUT "OVERRIDE_WARNINGS 0\n\n";
print OUT "ORDER_BY_RUN 1\n\n";
close OUT;
open OUT1,">$name\_k$cluster.sh";
print OUT1 "$clumpp $name\_k$cluster.paramfile\n";
close OUT1;
######## Sub Routines ##### Good Luck ########
###### No error ### No bug ### No warning #####


