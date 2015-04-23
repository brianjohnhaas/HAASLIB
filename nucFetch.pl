#!/usr/local/bin/perl

use strict;
use lib ($ENV{EUK_MODULES});
use LWP::UserAgent;
use Genbank_query;


main: {
    
    my $accession = $ARGV[0] or die "\n\nusage: $0 nucleotide_accession\n\n";
    
    my $ua = new LWP::UserAgent();
    
    my %URL_params = (db => 'nucleotide',
		      id => $accession,
		      retmode => 'text',
		      rettype => 'fasta');
    
    my $sequence = EFetch($ua, \%URL_params);
    
    if ($sequence =~ /^ERROR: Cannot/) {
	print STDERR "Sorry, couldn't retrieve sequence for $accession\n";
	exit(1);
    } else {
	print $sequence;
    }
    
    exit(0);
}
