#!/usr/local/bin/perl

use strict;
use BHStats;

my $usage =  "usage: $0 N_observations K_successes P_probability\n" ;

unless ($#ARGV == 2) { die $usage;}

my ($n,$k,$p) = (@ARGV);

my $bin_prob = &BHStats::binomial_probability($n,$k,$p);

print "Binomial Probability: $bin_prob\n";


exit(0);


