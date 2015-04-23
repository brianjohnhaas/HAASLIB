#!/usr/local/bin/perl

package ConfigFileReader;
use strict;

use vars qw (@ISA @EXPORT); ## set in script using this module for verbose output.
@ISA = qw(Exporter);
@EXPORT = qw(readConfig);

sub readConfig {
    my ($filename) = @_;
    open (FILE, $filename) or die "Can't open $filename\n$!";
    my %config;
    while (<FILE>) {
	chomp;
	if (/^\#/) { next;} #comment
	if (/=/) {
	    my ($key, $value) = split (/=/, $_, 2);
	    $config{$key} = $value;
	}
    }
    close FILE;
    return (%config);
}


1; #EOM
