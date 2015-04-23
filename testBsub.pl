#!/usr/bin/env perl

use strict;
use warnings;
use Bsub;

my $usage = "usage: $0 num_cmds_to_test\n\n";
my $num_cmds = $ARGV[0] or die $usage;

my @cmds;
for (1..$num_cmds) {
    my $cmd = "sleep 20";
    push (@cmds, $cmd);
}

my $bsubber = new Bsub({cmds=>\@cmds,
                    max_nodes => 200 ,
                    cmds_per_node => 1});

$bsubber->bsub_jobs();

exit(0);

