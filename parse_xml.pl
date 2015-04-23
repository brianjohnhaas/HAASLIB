#!/usr/local/bin/perl

use lib ("/usr/local/devel/Euk_modules/src");
use XML::Parser;
use Genbankxmlparser;
use strict;
use LWP::UserAgent;
use Getopt::Std;
use vars qw($opt_h $opt_a $opt_i $SEE);

$SEE = (@ARGV) ? 1:0; 

## Parse XML
my $pl = new XML::Parser (ErrorContext => 2, Style => "Tree");
my $tree = $pl->parsefile('-');

my $genbankparser = new Genbankxmlparser (tree => $tree, SEE => $SEE);

my @genes = $genbankparser->get_genes();
my $seq = $genbankparser->get_sequence();

my $x = 0;
foreach my $gene (@genes) {
    $x++;
    print "\nGENE $x\n";
    $gene->create_CDS_sequence (\$seq);
    print $gene->toString();
}

my ($locus, $accession, $version, $title, $update_date, $chromosome) = ($genbankparser->{locus},
									$genbankparser->{accession},
									$genbankparser->{version},
									$genbankparser->{title},
									$genbankparser->{update_date},
									$genbankparser->{chromosome});
					     
print "Locus: $locus\nAccession: $accession\nVersion: $version\nTitle: $title\nUpdate date: $update_date\nchromosome: $chromosome\n";



#print "SEQUENCE:\n$seq\n";

exit;


####
sub usage {
    
	print <<_EOH_;

############################# Options ###############################
# usage: $0 [-a|-i]
#
# -a genbank accession
# -i gi number
# -h print this option menu and quit
#
###################### Process Args and Options #####################

_EOH_

}

