#!/usr/local/bin/perl

use strict;
use LWP::UserAgent;
use Genbank_query;


our $SEE = 1;

my $ua = new LWP::UserAgent;


# EPost example.
#PubMed Example:
#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi?db=pubmed&id=11237011

my %vals = (db=>'pubmed',
	    id=>11237011);

my $text = EPost($ua, \%vals);

print $text;

exit(0);

