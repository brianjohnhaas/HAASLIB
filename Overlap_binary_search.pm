#!/usr/local/bin/perl
package main;
our $SEE;


package Overlap_binary_search;
use strict;


# Provide a list of elements in the form:
# @ordered_element_list = ( [lend, rend], [lend, rend], ...)
# where lend <= rend and elements are ordered by their midpoint.


sub new {
    my ($packagename, @ordered_element_list) = @_;
    
    my $self = { ordered_list => [@ordered_element_list] };
    
    bless ($self, $packagename);
    return ($self);
}




## Performs a binary search on the ordered element list to identify an element overlapping the 
## given set of lend,rend coordinates.
# returns the index of the overlapping element within the list, or -1 if none found.

sub find_overlapping_element {
    my ($self, $lend, $rend) = @_;
    
    my $midpt = ($lend+$rend)/2;

    my $ordered_element_list_aref = $self->{ordered_list};
    my $last_index = $#$ordered_element_list_aref;
    
    my $leftBound = 0;
    my $rightBound = $last_index;
    
    
    while ($leftBound <= $rightBound) {
	my $currIndex = int( ($leftBound + $rightBound) /2 );
	my $currElement = $ordered_element_list_aref->[$currIndex];
	
	my ($currLend, $currRend) = @$currElement;
	my $curr_mid = ($currLend + $currRend) /2;

	print "comparing ($lend,$rend) to ($currLend,$currRend)\n" if $SEE;


	if ($lend <= $currRend && $rend >= $currLend) { #overlap
	    return ($currIndex);
	}

	if ($midpt < $curr_mid ) {
	    $rightBound = $currIndex - 1;
	} else {
	    $leftBound = $currIndex + 1;
	}
    }

    return (-1); #not found
}


1; #EOM

