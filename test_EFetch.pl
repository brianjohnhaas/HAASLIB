#!/usr/local/bin/perl

use strict;
use LWP::UserAgent;
use Genbank_query;

our $SEE = 1;

my $ua = new LWP::UserAgent;

### Section I: Examples using the Sequence databases


#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=5

my %vals = (db=>'nucleotide',
	    id=>5
	    );

my $text = EFetch($ua, \%vals);

print $text;

sleep(3);

#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=5&complexity=0&rettype=fasta

%vals = (db=>'nucleotide',
	 id=>5,
	 complexity=>0,
	 rettype=>'fasta');


$text = EFetch($ua, \%vals);

print $text;

sleep(3);

#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=5&rettype=gb&seq_start=1&seq_stop=9

%vals = (db => 'nucleotide',
	 id=>5,
	 rettype=>'gb',
	 seq_start=>1,
	 seq_stop=>9);


$text = EFetch($ua, \%vals);

print $text;

sleep(3);


#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.cgi?db=nucleotide&id=5&rettype=fasta&seq_start=1&seq_stop=9&strand=2

%vals = (db=>'nucleotide',
	 id=>5,
	 rettype=>'fasta',
	 seq_start=>1,
	 seq_stop=>9,
	 strand=>2);

$text = EFetch($ua, \%vals);

print $text;

sleep(3);

#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=5&rettype=gb

%vals = (db=>'nucleotide',
	 id=>5,
	 rettype=>'gb');


$text = EFetch($ua, \%vals);

print $text;

sleep(3);


#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=popset&id=12829836&rettype=gp

%vals = (db=>'popset',
	 id=>12829836,
	 rettype=>'gp');


$text = EFetch($ua, \%vals);

print $text;

sleep(3);

#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=8&rettype=gp

%vals = (db=>'protein',
	 id=>8,
	 rettype=>'gp');


$text = EFetch($ua, \%vals);

print $text;

sleep(3);


#Entrez display format GBSeqXML:
#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=5&rettype=gb&retmode=xml

%vals = (db=>'nucleotide',
	 id=>5,
	 rettype=>'gb',
	 retmode=>'xml');


$text = EFetch($ua, \%vals);

print $text;

sleep(3);

#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=8&rettype=gp&retmode=xml

%vals = (db=>'protein',
	 id=>8,
	 rettype=>'gp',
	 retmode=>'xml');



$text = EFetch($ua, \%vals);

print $text;

sleep(3);


#Entrez display format TinySeqXML:

#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=5&rettype=fasta&retmode=xml 

%vals = (db=>'nucleotide',
	 id=>5,
	 rettype=>'fasta',
	 retmode=>'xml');

$text = EFetch($ua, \%vals);

print $text;

sleep(3);




#######################################################################################
## Section II: EFetch for Literature Databases


#In PubMed display PMIDs 12091962 and 9997 in html retrieval mode and abstract retrieval type:
#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=12345,9997&retmode=html&rettype=abstract

%vals = (db=>'pubmed',
	    id=>'12345,9997',
	    retmode=>'html',
	    rettype=>'abstract');


$text = EFetch($ua, \%vals);

print $text;

sleep(3);



#In PubMed display PMIDs in xml retrieval mode:
#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=11748933,11700088&retmode=xml

%vals = (db=>'pubmed',
	 id=>'11748933,11700088',
	 retmode=>'xml');

$text = EFetch($ua, \%vals);

print $text;

sleep(3);


#In Journals display records for journal IDs 22682,21698,1490:
#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=journals&id=22682,21698,1490&rettype=full 

%vals = (db=>'journal',
	 id=>'22682,21698,1490',
	 rettype=>'full');


$text = EFetch($ua, \%vals);

print $text;

sleep(3);




exit(0);
