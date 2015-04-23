#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Overlap_info;

my $usage = "usage: $0 lendA rendA lendB rendB\n\n";

if (scalar(@ARGV) < 4) {
	die $usage;
}

my ($lendA, $rendA, $lendB, $rendB) = @ARGV;

my $overlap_len = &Overlap_info::overlap_length([$lendA, $rendA], [$lendB, $rendB]);

print "$overlap_len\n";



exit(0);

