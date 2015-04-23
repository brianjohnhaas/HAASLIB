#!/usr/local/bin/perl

use strict;
use Fasta_reader;

my $usage = "usage: $0 fastaFile\n";


my $fastaFile = $ARGV[0] or die $usage;

my $fastaReader = new Fasta_reader($fastaFile);
my $counter = 0;
while (my $sequence = $fastaReader->next()) {
    my ($accession, $header, $sequence) = ($sequence->{accession},
					   $sequence->{header},
					   $sequence->{sequence});
    $counter++;

    print "# $counter\nacc: $accession\n"
	. "header: $header\n" 
	. "sequence length: " . length($sequence) . "\n\n";
}

exit(0);

