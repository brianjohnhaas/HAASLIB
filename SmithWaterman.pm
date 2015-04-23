package main;
our $SEE;

package SmithWaterman;
use strict;
use warnings;


## Defaults:

my $MATCH_DEFAULT = 2;
my $MISMATCH_DEFAULT = -3;
my $GAP_DEFAULT = -10;


## Constants:
my $UP = 1;
my $DIAG = 2;
my $LEFT = 3;


sub new {
    my $packagename = shift;

    my ($seqA_ref, $seqB_ref) = @_;
    unless (ref $seqA_ref eq "SCALAR" && ref $seqB_ref eq "SCALAR") {
        die "Error!!!! need two sequence references\n";
    }

    my $self = {
        seqA => $seqA_ref,
        seqA_len => length($$seqA_ref),
        
        seqB => $seqB_ref,
        seqB_len => length($$seqB_ref),
        
        match_score => $MATCH_DEFAULT,
        mismatch_score => $MISMATCH_DEFAULT,
        gap_score => $GAP_DEFAULT,

        scoring_matrix => [],
        max_alignment_score => 0,
        
        # alignment coordinates (i => seqA, j => seqB)
        start_i => -1,
        start_j => -1,
        end_i => -1,
        end_j => -1,
            


    };

    bless ($self, $packagename);

    $self->_init();
        
    return ($self);

}


sub _init {
    my $self = shift;

    my $x = $self->{seqA_len};
    my $y = $self->{seqB_len};

    my $scores = $self->{scoring_matrix};
    
    ## init scores
    for (my $i=0; $i <= $x; $i++) {
        for (my $j = 0; $j <= $y; $j++) {
            $scores->[$i][$j]= { score => 0,
                               prev => 0 };
        }
    }
}
    

sub align {
    my $self = shift;
    
    my $x = $self->{seqA_len};
    my $y = $self->{seqB_len};
    
    my $seqA = $self->{seqA};
    my $seqB = $self->{seqB};

    my $MATCH = $self->{match_score};
    my $MISMATCH = $self->{mismatch_score};
    my $GAP = $self->{gap_score};

    my $scores = $self->{scoring_matrix};
    

    my $max_i = 0;
    my $max_j = 0;
    my $total_max_score = 0;
    
    ## score smith waterman style
    for (my $i = 1; $i <= $x; $i++) {
        
        my $charA = substr($$seqA, $i-1, 1);
        
        for (my $j = 1; $j <= $y; $j++) {
            
            my $max_score = 0;
            my $max_dir = 0;
            
            my $charB = substr($$seqB, $j-1, 1);
            
            print "Comparing: $charA-$charB, ($i,$j)\n" if $SEE;
            
            if ($charA eq $charB) {
                ## match situation:
                my $score = $scores->[$i-1][$j-1]->{score} + $MATCH;
                if ($score > $max_score) {
                    $max_score = $score;
                    $max_dir = $DIAG;
                }
            } else {
                ## check mismatch situation
                my $score = $scores->[$i-1][$j-1]->{score} + $MISMATCH;
                if ($score > $max_score) {
                    $max_score = $score;
                    $max_dir = $DIAG;
                }
            }
            
            
            ## check gap-i (UP)
            my $score = $scores->[$i][$j-1]->{score} + $GAP;
            if ($score > $max_score) {
                $max_score = $score;
                $max_dir = $UP;
            }
            
            ## check gap-j (LEFT)
            $score = $scores->[$i-1][$j]->{score} + $GAP;
            if ($score > $max_score) {
                $max_score = $score;
                $max_dir = $LEFT;
            }
            
            
            print "max score: $max_score, $max_dir\n" if $SEE;
            
            if ($max_score > 0) {
                $scores->[$i][$j]->{score} = $max_score;
                $scores->[$i][$j]->{prev} = $max_dir;
                
                # track best score ever
                if ($max_score > $total_max_score) {
                    $total_max_score = $max_score;
                    $max_i = $i;
                    $max_j = $j;
                }
                
            }
            
        }
    }
    
    my $alignment_i = "";
    my $alignment_j = "";
    
    $self->{max_alignment_score} = $total_max_score;
    if ($total_max_score > 0) {

        # store terminal coordinates:
        $self->{end_i} = $max_i;
        $self->{end_j} = $max_j;
        
        my $start_i = $max_i;
        my $start_j = $max_j;
        # do traceback:
        my $prev = $scores->[$max_i][$max_j]->{prev};
        while ($prev != 0) {
            
            print "Tracing Alignment: (score: $total_max_score) $max_i,$max_j\n" if $SEE;
            
            $start_i = $max_i;
            $start_j = $max_j;

            my $char_i = substr($$seqA, $max_i-1, 1);
            my $char_j = substr($$seqB, $max_j-1, 1);
            
            if ($prev == $DIAG) {
                #aligned chars
                $alignment_i = $char_i . $alignment_i;
                $alignment_j = $char_j . $alignment_j;
                $max_i--;
                $max_j--;
            }
            
            elsif ($prev == $UP) {
                # gap in i
                $alignment_i = "-" . $alignment_i;
                $alignment_j = $char_j . $alignment_j;
                $max_j--;
            }
            
            elsif ($prev == $LEFT) {
                # gap in j
                $alignment_i = $char_i . $alignment_i;
                $alignment_j = "-" . $alignment_j;
                $max_i--;
            }
            
            $prev = $scores->[$max_i][$max_j]->{prev};
            
        }
        
        $self->{start_i} => $start_i;
        $self->{start_j} = $start_j;
    }
    
    return ($alignment_i, $alignment_j);
}
    


1; #EOM
