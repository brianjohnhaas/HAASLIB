#!/usr/local/bin/perl

package Mummer_remap;

use strict;
use warnings;
use Carp;

use Mummer_delta_parser;


our $VERBOSE = 0;
our $DEBUG = 0;

# constants
our $QUERY_TO_REFERENCE = 1;
our $REFERENCE_TO_QUERY = 2;


sub new {
    my $packagename = shift;
    my $delta_file = shift;

    my $delta_parser = new Mummer_delta_parser($delta_file);


    my $self = {
        accessions_to_alignments => {},
        delta_parser => $delta_parser,
        MSP_list => [],
        
    };
    
    bless ($self, $packagename);

    
    $self->_init();

    return ($self);

}

####
sub _init {
    my $self = shift;
    
    my $delta_parser = $self->get_delta_parser();
    
    my $accessions_to_alignments_href = $self->{accessions_to_alignments};

    my @alignments = $delta_parser->get_alignments();
    
    foreach my $alignment (@alignments) {
        my $reference_accession = $alignment->get_reference_accession();
        my $query_accession = $alignment->get_query_accession();

        my $alignment_list_aref = $accessions_to_alignments_href->{$reference_accession}->{$query_accession};
        unless (ref $alignment_list_aref) {
            $alignment_list_aref = $accessions_to_alignments_href->{$reference_accession}->{$query_accession} = [];
        }

        push (@$alignment_list_aref, $alignment);
        
    }
    
    
}



####
sub get_delta_parser {
    my $self = shift;
    return ($self->{delta_parser});
}

####
sub get_accession_to_alignments_mappings {
    my $self = shift;
    
    return ($self->{accessions_to_alignments});
}


####
sub _clear_MSP_list {
    my $self = shift;
    @{$self->{MSP_list}} = ();
    
    return;
}

####
sub _store_MSP_list {
    my $self = shift;
    my @msps = @_;
        
    @{$self->{MSP_list}} = @msps;
    
    return;
}

####
sub get_MSPs {
    my $self = shift;
    return (@{$self->{MSP_list}});
}

####
sub transform_coordinates {
    my $self = shift;
    my ($type, $ids_href, @coordinates) = @_;
    
    unless ($type == $REFERENCE_TO_QUERY || $type == $QUERY_TO_REFERENCE) {
        confess "fatal, don't understand 'type' function argument.  ";
    }
    
    my @alignments = $self->_get_alignments_from_ids($ids_href);
    
    my @converted_coordinates;
    my %msp_tracker; # track those msps used for coordinate conversions, and store them for later access

    ## clear the current msp list
    $self->_clear_MSP_list();
    
 
  COORDINATES:
    foreach my $coordinate (@coordinates) {
        my $alignment = $self->_find_longest_overlapping_coordinate_set($type, $coordinate, \@alignments);
        unless (ref $alignment) {
            push (@converted_coordinates, -1); # no alignment contains this coordinate in the span
            next COORDINATES;
        }
                
        ## got overlapping alignment
        ## find the msp containing the coordinate;
        my @msps = $alignment->get_MSPs();
        my $msp = $self->_find_longest_overlapping_coordinate_set($type, $coordinate, \@msps);
        unless (ref $msp) {
            # no msp contains this coordinate
            push (@converted_coordinates, -1);
            next COORDINATES;
        }
        
        $msp_tracker{$msp} = $msp;
        
        if ($Mummer_remap::VERBOSE) {
            print "Longest alignment MSP found spanning coordinate $coordinate:\n" 
                . $msp->toString() . "\n";
        }
        
        my ($reference_lend, $reference_rend) = $msp->get_reference_coords();
        my ($query_lend, $query_rend) = $msp->get_query_coords();
        my $orientation = $msp->get_orientation();
        
        my $adj_coordinate;
        if ($type == $REFERENCE_TO_QUERY) {
            
            my $delta = $coordinate - $reference_lend;
            
            if ($orientation eq '+') {
                $adj_coordinate = $query_lend + $delta;
                print "CoordConvert: case A\n" if $DEBUG;
            }
            else {
                # revcomp match
                $adj_coordinate = $query_rend - $delta;
                print "CoordConvert: case B\n" if $DEBUG;
            }
           
        }
        else {
            ## QUERY_TO_REFERENCE
                        
            if ($orientation eq '+') {
                my $delta = $coordinate - $query_lend;
                $adj_coordinate = $reference_lend + $delta;
                print "CoordConvert: case C\n" if $DEBUG;
            }
            else {
                # revcomp'd 
                my $delta = $query_rend - $coordinate;
                $adj_coordinate = $reference_lend + $delta;
                print "CoordConvert: case D\n" if $DEBUG;
            }
        }
        
        push (@converted_coordinates, $adj_coordinate);
    }
    
    $self->_store_MSP_list(values %msp_tracker);  #trick to get unique entries (using keys doesn't work due to string conversion)
    
    ## check calling context to make sure we do the right thing in our return value(s)
    if (scalar (@converted_coordinates) > 1) {
        if (wantarray) {
            return (@converted_coordinates);
        }
        else {
            confess "fatal, have list of converted coordinates but function called in scalar context ";
        }
    }
    else {
        if (wantarray) {
            return (@converted_coordinates);
        }
        else {
            return (scalar $converted_coordinates[0]);
        }
    }
}


###
sub _find_longest_overlapping_coordinate_set {
    my $self = shift;
    my ($type, $coordinate, $coord_sets_aref) = @_;
    
    my $longest_coordset = undef;
    my $longest_length = 0;

    foreach my $coordset (@$coord_sets_aref) {
        my $got_overlap_flag = 0;
        my $length;
        if ($type == $REFERENCE_TO_QUERY) {
            my ($reference_lend, $reference_rend) = $coordset->get_reference_coords();
            $length = $coordset->get_reference_coord_span_length();
            if ($reference_lend <= $coordinate && $coordinate <= $reference_rend) {
                ## overlap found
                $got_overlap_flag = 1;
            }
        }
        elsif ($type == $QUERY_TO_REFERENCE) {
            my ($query_lend, $query_rend) = $coordset->get_query_coords();
            $length = $coordset->get_query_coord_span_length();
            if ($query_lend <= $coordinate && $coordinate <= $query_rend) {
                ## overlap found
                $got_overlap_flag = 1;
            }
        }
        else {
            confess "fatal bug, type not recognized.";
        }
        
        if ($got_overlap_flag) {
            if ($length > $longest_length) {
                $longest_length = $length;
                $longest_coordset = $coordset;
            }
        }
    }

    return ($longest_coordset);
}


####
sub _get_alignments_from_ids {
    my $self = shift;
    my ($ids_href) = @_;
    
    my $reference_accession = $ids_href->{reference_accession} or confess "must specify reference_accession as key in hash_ref ";
    my $query_accession = $ids_href->{query_accession} or confess "must specify query_accession as key in hash_ref ";
    
    my $accessions_to_alignments_href = $self->get_accession_to_alignments_mappings();
    my $alignment_list_aref = $accessions_to_alignments_href->{$reference_accession}->{$query_accession};
    
    unless (ref $alignment_list_aref) {
        die "Error, couldn't find alignment based on ref: $reference_accession, query: $query_accession ";
    }
    
    return (@$alignment_list_aref);
    
    
}



1; #EOM
    
