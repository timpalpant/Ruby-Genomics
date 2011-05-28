#!/usr/bin/perl -w

# Steve Parker
# stephen.parker@nih.gov


# print the predicted cleavage pattern for a DNA sequence


use strict;
use Getopt::Long qw(GetOptions);
use PredictCleavagePatterns;
use Bio::SeqIO;




# get command line options:
my ($h, $fasta);
GetOptions ('h'   => \$h,
	    	'f=s' => \$fasta);
	


my $usage = <<USAGE;


    USAGE: perl print_cleavage_pattern_from_fasta.pl -f <fasta_file>


    options:
    -h  display this help message
    -f  fasta file name



    example usage:
    perl print_cleavage_pattern_from_fasta.pl -f fasta.fa



USAGE


# check command line arguments and display help if needed
unless ($fasta) { die "$usage"; }
if ($h) { die "$usage"; }



my $seqio_obj = Bio::SeqIO->new(-file   => "$fasta",
                                    -format => "fasta" );


while (my $seq = $seqio_obj->next_seq()) {
    my $header = $seq->id;
    my $seq    = $seq->seq;
        
    if (length($seq) < 7) {
    	print ">$header,prediction_not_possible\n";
    	next;
    }

    print ">$header";

    # predict the cleavage pattern and return the result in an array reference
    my $pattern = &predictSequenceCleavagePattern(uc $seq);




    for my $p (@{$pattern}) {
		if (defined($p)) {
	    	#$p = sprintf("%0.7f", $p);
	    	print ",$p";
		}
		else {
	    	print ",N";
		}
    }
    print "\n";
}


