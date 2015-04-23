#!/usr/local/bin/perl


use lib ($ENV{EUK_MODULES});
use strict;
use Jaccard_coefficient_cluster_resolver;

our $SEE = 0;
our $CLUSTERPATH = $ENV{EGC_UTILITIES} . "/cluster";

my $link_score = $ARGV[0];

my $input_file = $ARGV[1];

unless (defined($link_score) && ($link_score >= 0 && $link_score <= 1) && ($input_file) && -s $input_file) {
    die "usage: $0 link_score input_file_of_pairs\n";
}


open (FILE, $input_file) or die $!;
my @pairs;
while (<FILE>) {
    chomp;
    unless (/\w/) { next;} #blank line.
    my $inputLine = $_;
    $_ =~ s/^\s+//; #remove leading whitespace
    my ($a, $b, @whatever) = split (/\s+/);
    if (defined ($a) && defined($b)) {
	push (@pairs, [$a, $b]) if ($a ne $b);
    } else {
	print STDERR "ERROR processing line: $inputLine\n";
    }
}
close FILE;


my $resolver = new Jaccard_coefficient_cluster_resolver($link_score);
my @clusters = $resolver->resolve_clusters(@pairs);

if (@clusters) {
    print STDERR "Here are your clusters (resolved with link score: $link_score):\n";
    foreach my $cluster (@clusters) {
	print "@$cluster\n";
    }
} else {
    print STDERR "Sorry, no clusters were created.\n";
}

exit(0);

