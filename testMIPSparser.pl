#!/usr/local/bin/perl

use lib ($ENV{EUK_MODULES});
use MIPS_parser;

unless (@ARGV) {die;}

my $file = $ARGV[0];

my $parser = new MIPS_parser();
$parser->parse_dat_file ($file);
my @genes = $parser->get_genes();
my $x = 1;
foreach my $gene (@genes) {
    print "$x\n" . $gene->toString();
    $x++;
}

