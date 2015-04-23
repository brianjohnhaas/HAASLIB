#!/usr/local/bin/perl

package Mummer_delta_parser;

use strict;
use warnings;
use Carp;

our $VERBOSE = 0; ## for verbose 

####
sub new {
    my $packagename = shift;
    
    my ($delta_filename) = @_;
    
    unless (-s $delta_filename) {
        confess "Error, cannot locate $delta_filename or file size is zero ";
    }

    my $self = {
        delta_filename => $delta_filename,
        match_list => [], # container for all Mummer_match_objs
    };

    bless ($self, $packagename);

    $self->_parse_delta_file();

    return ($self);
}

####
sub get_delta_filename {
    my $self = shift;
    return ($self->{delta_filename});
}

####
sub get_alignments {
    my $self = shift;
    return (@{$self->{match_list}});
}


####
sub _parse_delta_file {
    my $self = shift;
    
    my $delta_file = $self->get_delta_filename();
    open (my $fh, $delta_file) or confess "Error, cannot open $delta_file";

    my $line = "";
    until ($line =~ /^>/) {
        $line = <$fh>; #prime the line reader
        chomp $line;
    }
    
  DELTA_PARSING:
    while ($line) {
        if ($line =~ /^>(.*)$/) {
            ## initiating alignment reading:
            my $seq_pair_defline = $1;
                        
            ## read alignment until get to zero, terminating the alignment
            my $alignment_string = "-1";
          ALIGNMENT_PARSING:
            while (defined ($alignment_string) && $alignment_string !~ /^>/) {
                $alignment_string = $line = <$fh>; ## also prime our $line
                if (!defined($alignment_string)) {
                    #done parsing file
                    last DELTA_PARSING;
                }
                chomp $alignment_string;
                if ($alignment_string =~ /^>/) {
                    # started reading a new alignment, done with this alignment
                    last ALIGNMENT_PARSING;
                }
                
                
                ## parse the indels
                
                my $indel_text = "";
                my $indel_string = ""; 
                while (defined($indel_string) && $indel_string ne "0") {
                    $indel_string = <$fh>;
                    chomp $indel_string;
                    
                    if ($indel_string ne '0') {
                        $indel_text .= "$indel_string,";
                    }
                }
                chop $indel_text; #remove trailing comma
                
                ## add alignment
                $self->_add_alignment_entry($seq_pair_defline, $alignment_string, $indel_text);
          
            } # end of individual alignment
            
        } 
        else {
            warn "line [$line] read and not parsed.\n";
        }
        
        unless ($line =~ /^>/) {
            ## keep reading
            $line = <$fh>;
            chomp $line;
        }
        
    }
}


####
sub _add_alignment_entry {
    my $self = shift;
    my ($seq_pair_defline, $alignment_string, $indel_string) = @_;
    my $mummer_alignment = Mummer_delta_parser::Mummer_alignment->new($seq_pair_defline, $alignment_string, $indel_string);

    my $match_list_aref = $self->{match_list};
    push (@$match_list_aref, $mummer_alignment);
 
    return;
}


###################
##################
package Mummer_delta_parser::Coordinate_set_pair;
# to be a base class for anything that needs accessions, coordinates, and orientations
use strict;
use warnings;
use Carp;

####
sub new {
    my $packagename = shift;
    my ($reference_accession, $reference_lend, $reference_rend, 
        $query_accession, $query_lend, $query_rend,
        $orientation) = @_;
 
    # make sure coords are properly ordered lend <= rend
    ($reference_lend, $reference_rend) = sort {$a<=>$b} ($reference_lend, $reference_rend);
    ($query_lend, $query_rend) = sort {$a<=>$b} ($query_lend, $query_rend);

    my $reference_coord_span_length = abs ($reference_rend - $reference_lend) + 1;
    my $query_coord_span_length = abs ($query_rend - $query_lend) + 1;

    
    my $self = {
        reference_accession => $reference_accession,
        reference_lend => $reference_lend,
        reference_rend => $reference_rend,
        reference_coord_span_length => $reference_coord_span_length,

        query_accession => $query_accession,
        query_lend => $query_lend,
        query_rend => $query_rend,
        query_coord_span_length => $query_coord_span_length,
        
        orientation => $orientation,
    };
    
    bless ($self, $packagename);
    
    $self->_verify_numeric($reference_lend, $reference_rend, $query_lend, $query_rend);
    $self->_verify_orientation($orientation);
    
    return ($self);
}

####
sub get_reference_accession {
    my $self = shift;
    return ($self->{reference_accession});
}

####
sub get_reference_coords {
    my $self = shift;
    my @coords = ($self->{reference_lend}, $self->{reference_rend});
    return (@coords);
}

####
sub get_reference_coord_span_length {
    my $self = shift;
    return ($self->{reference_coord_span_length});
}

####
sub get_query_accession {
    my $self = shift;
    return ($self->{query_accession});
}

####
sub get_query_coords {
    my $self = shift;
    my @coords = ($self->{query_lend}, $self->{query_rend});
    return (@coords);
}

####
sub get_query_coord_span_length {
    my $self = shift;
    return ($self->{query_coord_span_length});
}


####
sub get_orientation {
    my $self = shift;
    return ($self->{orientation});
}


####
sub toString {
    my $self = shift;
    
    my $text = "Ref: " . $self->get_reference_accession()
        . ", Query: " . $self->get_query_accession() . "\n"
        . join ("-", $self->get_reference_coords()) . ", "
        . join ("-", $self->get_query_coords()) 
        . " [" . $self->get_orientation() . "]";
    
    return ($text);
    
}


####
sub toAlignment {
    my $self = shift;
    my ($reference_seqref, $query_seqref) = @_;
    
    unless (ref $reference_seqref && ref $query_seqref) {
        confess "Fatal, need perl-references to sequence data for the refseq and queryseq ";
    }

    my ($refseq_lend, $refseq_rend) = $self->get_reference_coords();
    my ($query_lend, $query_rend) = $self->get_query_coords();
    
    my $refseq_subseq = substr($$reference_seqref, $refseq_lend - 1, $refseq_rend - $refseq_lend + 1);
    my $query_subseq = substr($$query_seqref, $query_lend -1 , $query_rend - $query_lend + 1);
    
    my $orient = $self->get_orientation();
    
    if ($orient eq '-') {
        $query_subseq = $self->_reverse_complement($query_subseq);
    }
        
    my $align_text = "";

    my $align_line_length = 50;
    
    for (my $i = 0; $i < length($refseq_subseq); $i+= $align_line_length) {
        my $refseq_align = lc substr($refseq_subseq, $i, $align_line_length);
        my $query_align = lc substr($query_subseq, $i, $align_line_length);
        
        my @refseq_align_chars = split (//, $refseq_align);
        my @query_align_chars = split (//, $query_align);
        for (my $j = 0; $j <= $#refseq_align_chars; $j++) {
            if ($refseq_align_chars[$j] eq $query_align_chars[$j]) {
                $refseq_align_chars[$j] = $query_align_chars[$j] = uc $query_align_chars[$j]; # make upper case
            }
        }
        
        $refseq_align = join ("", @refseq_align_chars);
        $query_align = join ("", @query_align_chars);
        
        my $refseq_coord = $refseq_lend + $i;
        my $query_coord = ($orient eq '+') ? $query_lend + $i : $query_rend - $i;
        
        $align_text .= sprintf ("%12d $refseq_align\n%12d $query_align\n\n", $refseq_coord, $query_coord);      
        
    }
    
    return ($align_text);
    
}
    


####
sub _reverse_complement {
    my $self = shift;
    my ($sequence) = @_;
    my $rc = reverse ($sequence);
    
    $rc =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
    
    return($rc);
}                  


####
sub _verify_numeric {
    my $self = shift;
    foreach my $number (@_) {
        unless ($number =~ /^\d+$/) {
            confess "fatal, coordinate [$number] is not fully numeric or not an integer ";
        }
    }
    return;
}

####
sub _verify_orientation {
    my $self = shift;
    my ($orientation) = @_;
    
    if ($orientation eq '+' || $orientation eq '-') {
        ## good!
        return;
    }
    else {
        confess "fatal, orientation [$orientation] is not allowed.  Only +|-\n";
    }
}


#########################
#########################

package Mummer_delta_parser::Mummer_alignment;
use strict;
use warnings;
use Carp;
use base qw (Mummer_delta_parser::Coordinate_set_pair);

####
sub new {
    my $packagename = shift;
    my ($seqpair_line, $alignment_line, $indel_text) = @_;
    
    my ($reference, $query, $reference_length, $query_length) = split (/\s+/, $seqpair_line);
    my ($ref_start, $ref_end, $query_start, $query_end, @rest) = split (/\s+/, $alignment_line);
    
    my $orientation = ($query_start < $query_end) ? '+' : '-';
    
    
    my $self = $packagename->SUPER::new($reference, $ref_start, $ref_end,
                                        $query, $query_start, $query_end,
                                        $orientation);
    
    ## additional class-specific attributes
    
    $self->{refseq_match_length} =  abs ($ref_end - $ref_start) + 1;
    $self->{query_match_length} = abs ($query_end - $query_start) + 1;
    $self->{indels} = $indel_text;
    $self->{MSP_list} = []; # list of ungapped matches in this larger gapped alignment
    
    
    if ($Mummer_delta_parser::VERBOSE) {
        print "Got seqpairline: $seqpair_line\n"
            . "alignment_line: $alignment_line\n"
            . "and indels:\n$indel_text\n\n";
    }
    
    bless ($self, $packagename);
    
    $self->_parse_MSPs();

    if ($Mummer_delta_parser::VERBOSE) {
        print $self->toString();
    }

    return ($self);
}

sub _add_MSP {
    my $self = shift;
    my ($a_lend, $a_rend, $b_lend, $b_rend, $orient) = @_;
    
    my $msp = Mummer_delta_parser::Coordinate_set_pair->new($self->get_reference_accession(), $a_lend, $a_rend, 
                                                            $self->get_query_accession(), $b_lend, $b_rend, 
                                                            $orient);
    
    my $msp_list = $self->{MSP_list};
    push (@$msp_list, $msp);
    
    return;
}

####
sub get_MSPs {
    my $self = shift;
    return (@{$self->{MSP_list}});
}

####
sub get_indels {
    my $self = shift;
    my $indel_text = $self->{indels};
    unless ($indel_text =~ /\d/) {
        return (); # no indels to report
    }
 
    my @indels = split (/,/, $indel_text);
    
    unless (@indels) {
        confess "fatal, should have indels but none found"; #defensive progmg
    }

    return (@indels);
}

#### private
sub _parse_MSPs {
    my $self = shift;
    
    ## walk the indels
    my @indels = $self->get_indels();
    
    ## get some alignment properties:
    my $ref_match_length = $self->{refseq_match_length};
    my $query_match_length = $self->{query_match_length};
    my ($query_lend, $query_rend) = $self->get_query_coords();

    my ($refseq_start, $refseq_end) = $self->get_reference_coords(); 
    my $orientation = $self->get_orientation();
    
    my @msps; # coords reference to this alignment substring (starting at 1)

    if (@indels) {
        ## determine MSPs between the indels
        # init cursor pos:
        my $prev_ref_cursor = 0;
        my $prev_query_cursor = 0;
        
        while (@indels) {
            my $indel = shift @indels;
            my $delta = abs ($indel);
            $delta--;
            
            if ($delta > 0) {
                my $new_ref_cursor = $prev_ref_cursor + $delta;
                my $new_query_cursor = $prev_query_cursor + $delta;
                
                
                push ( @msps, 
                       [ $prev_ref_cursor + 1 . "-$new_ref_cursor", 
                         $prev_query_cursor + 1 . "-$new_query_cursor" ] 
                       );
                
                $prev_ref_cursor = $new_ref_cursor;
                $prev_query_cursor = $new_query_cursor;
            }
            
            if ($indel > 0) {
                ## insert in query
                $prev_ref_cursor++;
                # print "-insert in query\n";
                
            }
            else {
                ## insert in ref
                $prev_query_cursor++;
                # print "-insert in reference.\n";
            }
        }
        
        ## add last MSP beyond the last indel
                        
        if ($ref_match_length != $prev_ref_cursor && $query_match_length != $prev_query_cursor) {
            push (@msps, 
                  [ $prev_ref_cursor + 1 . "-$ref_match_length", 
                    $prev_query_cursor + 1 . "-$query_match_length"]
                  );
        }
    }
    else {
        ## no indels
        push (@msps, ["1-$ref_match_length", "1-$query_match_length"]); # the entire alignment is an MSP!
    }

    ## convert msp coordinates so they reflect the complete sequence's coordinates, and alignment orientation

    foreach my $msp (@msps) {
        my ($ref_msp_text, $query_msp_text) = @$msp;
        my ($r_start, $r_end) = split (/-/, $ref_msp_text);
        my ($q_start, $q_end) = split (/-/, $query_msp_text);

        ## adjust ref msp coords:
        my $r_start_adj = $refseq_start + $r_start - 1;
        my $r_end_adj = $refseq_start + $r_end - 1;

        ## adjust query msp coords:
        my ($q_start_adj, $q_end_adj);
        if ($orientation eq '+') {
            $q_start_adj = $query_lend + $q_start - 1;
            $q_end_adj = $query_lend + $q_end - 1;
        }
        else {
            # revcomp'd
            $q_start_adj = $query_rend - $q_start + 1;
            $q_end_adj = $query_rend - $q_end + 1;
            
            if ($q_start_adj <= 0 || $q_end_adj <= 0) {
                confess "Error, ended up with negative coordinates for query after adjustment!\n"
                    . "q_start: $q_start, q_end: $q_end\n"
                    . "query_start: $query_lend, query_end: $query_rend\n";
            }
            
            ## swap them so coords are ascending, orient var specifies strand, not coord vals:
            ($q_start_adj, $q_end_adj) = ($q_end_adj, $q_start_adj);
        }

        $self->_add_MSP($r_start_adj, $r_end_adj, $q_start_adj, $q_end_adj, $orientation);

    }
    
    
}

####
sub toString {
    my $self = shift;
    my $text = "Mummer_delta_parser::Mummer_alignment\n"
        . "(Ref: " 
        . $self->get_reference_accession() 
        . " vs. Query: " 
        . $self->get_query_accession() . ")\n"
        . join ("-", $self->get_reference_coords()) . ", "
        . join ("-", $self->get_query_coords()) . " " 
        . "[" . $self->get_orientation() . "]\nMSPs:\n";
        
    foreach my $msp ($self->get_MSPs()) {
        my $msp_text = $msp->toString();
        $msp_text =~ s/\n/\n\t/g;
        $text .= "\tmsp:\t$msp_text\n";
    }
    
    return ($text);
}

####
sub toAlignment {
    my $self = shift;
    my ($reference_seqref, $query_seqref) = @_;

    my $alignment_text = "";
    foreach my $msp ($self->get_MSPs()) {
        $alignment_text .= $msp->toString() . "\n";
        $alignment_text .= "\n" . $msp->toAlignment($reference_seqref, $query_seqref);
    }

    return ($alignment_text);
}
   


1; #EOM
