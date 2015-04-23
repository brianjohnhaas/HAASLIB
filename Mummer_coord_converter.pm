#!/usr/local/bin/perl

=Description

    This module performs coordinate transformations from one assembly to another
    based on the results from a Mummer comparison.

    #!/usr/local/perl

    use Mummer_coord_converter;
 
    my $converter = new Mummer_coord_converter();
    $converter->parse_Mummer_outputfiles('mum_forward', 'mum_reverse');
    $converter->data_dumper();
    my $num = 1;
    my $converted_num = $converter->transform_coordinate_1_to_2($num);
    print "$num converted to $converted_num\n";

=cut

our $SEE = 0; #set in main.

package Mummer_coord_converter;

use strict;



sub new {
    # no input parameters
    my $self = {
	mum_objs => []
    };
    bless ($self);
    return ($self);
}


sub parse_Mummer_outputfiles {
    my $self = shift;
    my ($forward_run, $reverse_run) = @_;
    
    ## parse forward run file
    open (MUM, "$forward_run") or die "Can't open $forward_run\n";
    while (<MUM>) {
	if (/^\#/) {
	    next;
	}
	if (/\s*(\d+)\s+(\d+)\s+(\d+)\s*/) {
	    print if $main::SEE;
	    my ($coord_seq1, $coord_seq2, $match_length) = ($1, $2, $3);
	    my $coord_seq1a = $coord_seq1;
	    my $coord_seq1b = $coord_seq1 + $match_length - 1;
	    my $coord_seq2a = $coord_seq2;
	    my $coord_seq2b = $coord_seq2 + $match_length - 1;
	    my $orientation = '+';
	    my $mum_obj = Mum_obj::create_Mum_obj ($coord_seq1a, $coord_seq1b, 
						   $coord_seq2a, $coord_seq2b, 
						   $match_length, $orientation);
	    $self->add_mum_obj($mum_obj);
	}
    }
    close MUM;

    ## parse reverse run file
    ## first pass, gather length of seq 1 (which was reverse complemented)
    my $seq1_length;
    open (MUM, "$reverse_run") or die "Can't open $reverse_run\n";
    while (<MUM>) {
	if (/^\#\s(\S+)\s(\d+)\s(\S+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)/) {
	    my $seq1_name = $1;
	    $seq1_length = $2;
	    my $seq2_name = $3;
	    my $min_mum_length_setting = $4;
	    my $num_mums = $5;
	    my $total_mums_length_in_nts = $5;
	    last;
	}
    }
    close MUM;
    unless ($seq1_length) { 
	die "ERROR: I couldn't parse the seq1_length from the mum_r file\n";
    }
    ## second pass, retrieve the mum info
    open (MUM, "$reverse_run") or die "Can't open $reverse_run\n";
    while (<MUM>) {
	if (/^\#/) {
	    next;
	}
	if (/\s*(\d+)\s+(\d+)\s+(\d+)\s*/) {
	    print if $main::SEE;
	    my ($coord_seq1, $coord_seq2, $match_length) = ($1, $2, $3);
	   
	    my $coord_seq1a = $coord_seq1;
	    my $coord_seq1b = $coord_seq1 + $match_length - 1;
	    ## swap positions after revcomping them
	    ($coord_seq1a, $coord_seq1b) = reverse (&revcomp_coord($coord_seq1a, $seq1_length),
						    &revcomp_coord($coord_seq1b, $seq1_length));


	   
	    my $coord_seq2a = $coord_seq2;
	    my $coord_seq2b = $coord_seq2 + $match_length - 1;
	    my $orientation = '-';
	    my $mum_obj = Mum_obj::create_Mum_obj ($coord_seq1a, $coord_seq1b, 
						   $coord_seq2a, $coord_seq2b, 
						   $match_length, $orientation);
	    $self->add_mum_obj($mum_obj);
	}
    }
    close MUM;
    
    $self->sort_MUMs();
}

sub sort_MUMs {
    my $self = shift;
    my $mum_objs_list = $self->{mum_objs};
    @$mum_objs_list = reverse sort {$a->{match_length}<=>$b->{match_length}} @$mum_objs_list;
}






sub add_mum_obj {
    my $self = shift;
    my $mum_obj = shift;
    my $index = $#{$self->{mum_objs}};
    $index++;
    #print "Index: $index\n";
    $self->{mum_objs}->[$index] = $mum_obj;
}

sub transform_coordinate_1_to_2 {
    my $self = shift;
    my $coord = shift;
    my ($mum_obj, $c1, $c2) = $self->get_spanning_mum_obj($coord);
    if ($mum_obj) {
	my $strand = $mum_obj->{orientation};
	my $rel_coord2 = $mum_obj->{coord_seq2a};
	if ($strand eq '+') {
	    my $delta = $coord - $c1;
	    return ($delta + $rel_coord2);
	} else {
	    my $delta = $c2 - $coord;
	    return ($delta + $rel_coord2);
	}
    }
}

sub get_spanning_mum_obj {
    my $self = shift;
    my $coord = shift;
    my $mum_collection = $self->{mum_objs};
    for (my $i=0; $i <=$#{$mum_collection}; $i++) {
	my $c1 = $mum_collection->[$i]->{coord_seq1a};
	my $c2 = $mum_collection->[$i]->{coord_seq1b};
	#print "$i Searching coords: $coord >= $c1 && $coord <= $c2\n";
	if ($coord >= $c1 && $coord <= $c2) {
	    #print "FOUND IT\n";
	    return ($mum_collection->[$i], $c1, $c2);
	}
    }
    #print "DAMNIT!!!!\n";
    return (0);
}


sub data_dumper {
    my ($self) = shift;
    my @mum_objs = @{$self->{mum_objs}};
    foreach my $mum_obj (@mum_objs) {
	print $mum_obj->toString();
    }
}

sub get_all_mums {
    my $self = shift;
    return (@{$self->{mum_objs}});
}


sub get_longest_mum {
    my $self = shift;
    my @mum_objs = @{$self->{mum_objs}};
    my $longest_length = 0; #initialize
    my $longest_mum = undef(); #initialize
    foreach my $mum (@mum_objs) {
	print $mum->toString() if $main::SEE;
	my $match_length = $mum->get_match_length();
	if ($match_length > $longest_length) {
	    $longest_length = $match_length;
	    $longest_mum = $mum;
	}
    }
    return ($longest_mum);
}



sub revcomp_coord {
    my ($coord, $seq_length) = @_;
    return ($seq_length - $coord + 1);
}





##################################
## Mum_obj stores mummer data

package Mum_obj;
use strict;

sub new {
    my $self = {
	coord_seq1a=>0,
	coord_seq1b=>0,
	coord_seq2a=>0,
	coord_seq2b=>0,
	match_length=>0,
	orientation=> 0
	};
    bless ($self);
    return ($self);
}

sub create_Mum_obj {
    my ($coord_seq1a, $coord_seq1b, $coord_seq2a, $coord_seq2b, 
	$match_length, $orientation) = @_;
    print "INCOMING: ($coord_seq1a, $coord_seq1b, $coord_seq2a, $coord_seq2b, $match_length, $orientation\n" if $SEE;
    my $self = new();
    $self->{coord_seq1a} = $coord_seq1a;
    $self->{coord_seq1b} = $coord_seq1b;
    $self->{coord_seq2a} = $coord_seq2a;
    $self->{coord_seq2b} = $coord_seq2b;
    $self->{match_length} = $match_length;
    $self->{orientation} = $orientation;
    return ($self);
}

sub get_coords_seq1 {
    my $self = shift;
    return ($self->{coord_seq1a}, $self->{coord_seq1b});
}

sub get_coords_seq2 {
    my $self = shift;
    return ($self->{coord_seq2a}, $self->{coord_seq2b});
}

sub get_match_length {
    my $self = shift;
    return ($self->{match_length});
}

sub get_orientation {
    my $self = shift;
    return ($self->{orientation});
}

sub toString {
    my $self = shift;
    return "MUM Data: $self->{coord_seq1a}\t$self->{coord_seq1b}\t$self->{coord_seq2a}\t$self->{coord_seq2b}\t$self->{match_length}\t$self->{orientation}\n";
}



1; #end of module

