#!/usr/local/bin/perl

use strict;
use LWP::UserAgent;
use Genbank_query;

our $SEE = 1;

my $ua = new LWP::UserAgent;

#Example:
#In Protein display records for GIs 28800982 and 28628843 in xml retrieval mode:
#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=protein&id=28800982,28628843&retmode=xml

my %vals = ( db=>'protein',
	     id=>'28800982,28628843',
	     retmode=>'xml');

my $text = ESummary($ua, \%vals);

print $text;

sleep(3);




#In Nucleotide display records for GIs 28864546 and 28800981 in xml retrieval mode:
#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nucleotide&id=28864546,28800981&retmode=xml

%vals = (db=>'nucleotide',
	 id=>'28864546,28800981',
	 retmode=>'xml');


$text = ESummary($ua, \%vals);

print $text;

sleep(3);


#In Structure display records for MMDB IDs 19923 and 12120 in xml retrieval mode:
#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=structure&id=19923,12120&retmode=xml

%vals = (db=>'structure',
	 id=>'19923,12120',
	 retmode=>'xml');

$text = ESummary($ua, \%vals);

print $text;

sleep(3);


#In Taxonomy display records for TAXIDs 9913 and 30521 in xml retrieval mode:
#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy&id=9913,30521&retmode=xml

%vals = (db=>'taxonomy',
	 id=>'9913,30521',
	 retmode=>'xml');


$text = ESummary($ua, \%vals);

print $text;

sleep(3);


#In UniSTS display records for IDs 254085 and 254086 in xml retrieval mode:
#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=unists&id=254085,254086&retmode=xml

%vals = (db=>'unists',
	 id=>'254085,254086',
	 retmode=>'xml');

$text = ESummary($ua, \%vals);

print $text;


exit(0);

