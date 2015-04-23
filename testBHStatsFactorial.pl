#!/usr/local/bin/perl

use strict;
use BHStats;

my $int = $ARGV[0] or die "usage: $0 integer\n" ;

my $fact = BHStats::factorial($int);
print "\nFACT: $fact\n";




exit(0);
