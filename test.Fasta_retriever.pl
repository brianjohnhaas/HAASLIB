#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;

use lib ("$FindBin::Bin");
use Fasta_retriever;
use Data::Dumper;

my $usage = "usage: $0 db.fasta\n\n";


my $db_fasta = $ARGV[0] or die $usage;


my $fasta_retriever = new Fasta_retriever($db_fasta);

print STDERR Dumper($fasta_retriever);


exit(0);

