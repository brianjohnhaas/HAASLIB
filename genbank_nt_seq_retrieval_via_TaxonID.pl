#!/usr/local/bin/perl

use strict;
use lib ($ENV{EUK_MODULES});
use Eutilities_ncbi;
use Data::Dumper;

my $usage = "usage: $0 taxonID\n\nScript retrieves all corresponding nucleotide entries from genbank and writes them to a file called taxon_{taxon_id}.fasta\n\n";

my $taxon_id = $ARGV[0] or die $usage;

## Get the list of GI's
my %params = (db => 'nucleotide',
	      term => "txid${taxon_id}" . "[orgn]",
	      usehistory => 'y',
	      tool => "$0",
	      email => "bhaas\@tigr.org");

my %results = &esearch(%params);

## Fetch the sequences:

open (OUTPUT, ">taxon_${taxon_id}.fasta") or die "Cannot write output file.\n";

my %efetch_params = (db => 'nucleotide',
		     id => $results{uids},
		     query_key => $results{query_key},
		     WebEnv => $results{WebEnv},
		     retmode => "text",
		     rettype => "FASTA",
		     tool => $results{tool},
		     email => $results{email},
		     filehandle => *OUTPUT);

&efetch_batch(%efetch_params);
		      
exit(0); ## Beware, no error checking done here.  Nothing may have worked. :(, or maybe it all worked :)



