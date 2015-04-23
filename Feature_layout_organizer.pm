#!/usr/local/bin/perl

package main;
our ($DEBUG);

package Feature_layout_organizer;

use strict;

## Each feature must be a hashref with keys (lend, rend) where lend <= rend


sub organize_features {
    my @features = @_;
    
    ## First, sort the features so that the longer ones appear first:
    @features = reverse sort { ($a->{rend} - $a->{lend}) 
				   <=>
				   ($b->{rend} - $b->{lend}) } @features;
    


    
    my @rows = ([]);
    # each row contains an array-ref of sorted non-overlapping features
    
    while (@features) {
	my $curr_feature = shift @features;
	my $added_flag = 0;
	for (my $i=0; $i <= $#rows; $i++) {
	    my $curr_row = $rows[$i];
	    if (! &overlaps_existing_features($curr_row, $curr_feature)) {
		# add feature to current row:
		push (@$curr_row, $curr_feature);
		@$curr_row = sort {$a->{lend}<=>$b->{lend}} @$curr_row; #maintain sortedness.
		$added_flag = 1;
		last;
	    }
	}
	unless ($added_flag) {
	    #must add another row:
	    push (@rows, [$curr_feature]);
	}
    }
    return (@rows);
}



sub overlaps_existing_features {
    my ($curr_row, $curr_feature) = @_;
    
    ## check for empty row:
    unless (@$curr_row) {
	return (0); #nothing yet to overlap.
    }

    ## Do a binary search on the row to see if features overlap:
    my $curr_lend = $curr_feature->{lend};
    my $curr_rend = $curr_feature->{rend};
    
    my $left_index = 0;
    my $right_index = $#{$curr_row};
    
    my $done = 0;
    while (!$done) {
	if ($left_index > $right_index) { 
	    $done = 1;
	} else {
	    my $search_index = int (($left_index + $right_index)/2 + 0.5);
	    my $check_feature = $curr_row->[$search_index];
	    my ($check_lend, $check_rend) = ($check_feature->{lend}, $check_feature->{rend});
	    if ($check_lend <= $curr_rend && $check_rend >= $curr_lend) {
		# overlaps:
		return (1);
	    } else {
		# no overlap, reset search range and keep looking.
		if ($curr_lend < $check_lend) {
		    $right_index = $search_index - 1;
		} else {
		    $left_index = $search_index + 1;
		}
	    }
	}
    }
    return (0); # no apparent conflict
}





1; #EOM
