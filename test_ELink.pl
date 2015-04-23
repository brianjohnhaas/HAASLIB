#!/usr/local/bin/perl

use strict;
use LWP::UserAgent;
use Genbank_query;


our $SEE = 1;

my $ua = new LWP::UserAgent;


#To retrieve IDs and relevancy scores from PubMed for PMID
#9298984 to the PubMed database:
#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&id=9298984&cmd=neighbor
 
my %vals = (dbfrom=>'pubmed',
	    id=>9298984,
	    cmd=>'neighbor');


my $text = ELink($ua, \%vals);

print $text;

sleep(3);

   
#To retrieve IDs from Nucleotide for GI 18250303, 18250301,
#18250299 to Protein:
#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nucleotide&db=protein&id=18250303,18250307


%vals = (dbfrom=>'nucleotide',
	 db=>'protein',
	 id=>'18250303,18250307');

$text = ELink($ua, \%vals);

print $text;

sleep(3);

   

    
#To retrieve PubMed related articles for PMIDs 11812492
#11774222 with a publication date from 1995 to the present:
#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&id=11812492,11774222&db=pubmed&mindate=1995&datetype=pdat

%vals = (dbfrom=>'pubmed',
	 id=>'11812492,11774222',
	 db=>'pubmed',
	 mindate=>1995,
	 datetype=>'pdat');


$text = ELink($ua, \%vals);

print $text;

sleep(3);

   

    
#To retrieve MEDLINE indexed only related articles for PMID
#12242737:
#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&id=12242737&db=pubmed&term=medline[sb]


%vals = (dbfrom=>'pubmed',
	 id=>12242737,
	 db=>'pubmed',
	 term=>'medline[sb]');


$text = ELink($ua, \%vals);

print $text;

sleep(3);

   

#To create a hyperlink to the first link available for PMID
#10611131 in PubMed:
#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&id=10611131&retmode=ref&cmd=prlinks

%vals = (dbfrom=>'pubmed',
	 id=>10611131,
	 retmode=>'ref',
	 cmd=>'prlinks');


$text = ELink($ua, \%vals);

print $text;

sleep(3);

   

#To list all available links in PubMed, except for
#libraries, for PMIDs 12085856 and 12085853:
#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&id=12085856,12085853&cmd=llinks

%vals = (dbfrom=>'pubmed',
	 id=>'12085856,12085853',
	 cmd=>'llinks');



$text = ELink($ua, \%vals);

print $text;

sleep(3);

   

#To check for the existence of a Related Articles link for
#PMIDs 0611131, 111645 and 12068369:
#http://eutils.ncbi.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&id=10611131+111645&id=12068369&cmd=ncheck

%vals = (dbfrom=>'pubmed',
	 id=>'10611131,111645,12068369',
	 cmd=>'ncheck');


$text = ELink($ua, \%vals);

print $text;

sleep(3);

   
exit(0);

