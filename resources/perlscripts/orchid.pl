#!/usr/bin/perl -w

# Steve Parker
# stephen.parker@nih.gov


# print the predicted cleavage pattern for a DNA sequence



use strict;
use Getopt::Long qw(GetOptions);
use PredictCleavagePatterns;





# get command line options:
my ($h, $sequence);
GetOptions ('h'   => \$h,
	    	's=s' => \$sequence);
	


my $usage = <<USAGE;


    USAGE: perl print_cleavage_pattern.pl -s <sequence>


    options:
    -h  display this help message
    -s  DNA sequence
	

    
    
    NOTES:
    This script will print out an input DNA sequence, one base per line,
    with the predicted cleavage intensity associated with that base.

    No N's allowed (only A, C, G, or T)


    example usage:
    perl print_cleavage_pattern.pl -s ACGTACGATCGACTAGCATCGACT



USAGE


# check command line arguments and display help if needed
unless ($sequence) { die "$usage"; }
if ($h) { die "$usage"; }




# split the bases of the input sequence into an array
my @seq = split '', $sequence;

# predict the cleavage pattern and return the result in an array reference
my $pattern = &predictSequenceCleavagePattern(uc $sequence);


for my $p (@{$pattern}) {
    my $base = shift @seq;
    print "$base\t$p\n";
}		

exit;

