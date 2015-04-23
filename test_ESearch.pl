#!/usr/local/bin/perl

use strict;
use LWP::UserAgent;
use Genbank_query;


our $SEE = 1;

my $ua = new LWP::UserAgent;

######################################3
## ESearch examples:

#Search in PubMed for the term cancer for the entrez date from the last 60 days and retrieve the first 100 IDs and translations using the history parameter:
#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=cancer&reldate=60&datetype=edat&retmax=100&usehistory=y


my %vals = (db=>'pubmed',
	    term=>'cancer',
	    reldate=>60,
	    datetype=>'edat',
	    retmax=>100,
	    usehistory=>'y');

my $text = ESearch($ua, \%vals); 

print $text;

sleep(3);

#Search in PubMed for the journal PNAS Volume 97, and retrieve 6 IDs starting at ID 7 using a tool parameter:
#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=PNAS[ta]+AND+97[vi]&retstart=6&retmax=6&tool=biomed3

%vals = (db=>'pubmed',
	 term=>"PNAS[ta] AND 97[vi]",
	 retstart=>6,
	 retmax=>6,
	 tool=>'biomed3');

$text = ESearch($ua, \%vals);

print $text;

sleep(3);


#Search in Journals for the term obstetrics:
#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=journals&term=obstetrics

%vals = ( db=>'journal',
	  term=>'obstetrics');

$text = ESearch($ua, \%vals);

print $text;


exit(0);

