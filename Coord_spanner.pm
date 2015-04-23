package Coord_spanner;

use strict;
use warnings;
use Carp;

## provide a list of coordinates like so: ([end5, end3], [end5, end3], ...)
## returns exclusive regions of overlap: ([lend, rend], [lend, rend], ...)
sub get_coord_spans {
    my @coord_sets = @_;
    
    unless (@coord_sets) { confess "need coordsets as param"; }

    my @coords;
    ## convert end5, end3 to lend,rend w/o strand information.
    foreach my $coordset (@coord_sets) {
        my ($lend, $rend) = sort {$a<=>$b} @$coordset;
        push (@coords, [$lend, $rend]);
    }
    
    # put coordsets in ascending order
    @coords = sort {$a->[0]<=>$b->[0]} @coords;

    
    ## determine non-overlapping maximal spans:
    my @spans = shift @coords; # prime it w/ the first pair
    
    while (@coords) {
        my $prev_coordset = $spans[$#spans];
        my $next_coordset = shift @coords;

        my ($prev_lend, $prev_rend) = @$prev_coordset;
        my ($next_lend, $next_rend) = @$next_coordset;

        if ($next_lend <= $prev_rend) {
            if ($next_rend > $prev_lend) {
                $prev_coordset->[1] = $next_rend; # reset right boundary
            }
        }
        else {
            # don't overlap, so start new span
            push (@spans, $next_coordset);
        }
    }

    return (@spans);
}


####
sub sum_coord_spans {
    my @coordsets = @_;
    
    my $sum_length = 0;
    
    foreach my $coordset (@coordsets) {
        my ($lend, $rend) = @$coordset;
        
        my $length = abs ($rend - $lend) + 1;

        $sum_length += $length;
    }

    return ($sum_length);
}



1; #EOM

    
    
    

    
