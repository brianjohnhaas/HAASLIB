#!/usr/local/bin/perl

use strict;
use BHStats;

unless (@ARGV) {
    die "usage: $0 number_list\n\n";
}

my $avg = BHStats::avg(@ARGV);

print "avg(@ARGV) = $avg\n\n";

exit(0);
