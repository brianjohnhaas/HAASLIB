#!/usr/local/bin/perl

use XML::Parser;
use Genbankxmlparser;
use strict;
use LWP::UserAgent;
use Getopt::Std;
use vars qw($opt_h $opt_a $opt_i); 

&getopts ('ha:i:');

if ($opt_h || !($opt_a | $opt_i)) {&usage();exit;}
my ($param, $param_type);
if ($opt_a) {
    $param = $opt_a;
    $param_type = "Search";
} else {
    $param = $opt_i;
    $param_type = "Text";
}

## Get XML from Genbank
my $ua = new LWP::UserAgent;
my $req = new HTTP::Request (GET => "http://www.ncbi.nlm.nih.gov/entrez/viewer.cgi?save=0&cmd=$param_type&cfm=on&view=xml&txt=on&val=$param");
my $res = $ua->request($req);
my $xml = $res->content;

print $xml;
exit;
 
sub usage {
    die <<_EOH_;

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

