#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;

use lib ("$FindBin::Bin/");
use ColorGradient;


my $usage = "usage: $0 numColors\n\n";

my $num_colors = $ARGV[0] or die $usage;

my @colors = &ColorGradient::get_RGB_gradient($num_colors);

my @hex_colors = &ColorGradient::convert_RGB_hex(@colors);

print join ("\n", @hex_colors) . "\n";


exit(0);

