#!/usr/local/bin/perl

package ORF_finder;

use strict;
use Nuc_translator;


sub new {
    shift;

    ## This object stores the longest ORF identified.
    my $obj = { pep_seq =>0,
		nt_seq => 0,
		length => 0, #length of nt_seq
		end5 => 0,
		end3 => 0,
		all_ORFS=>[], #container holds all ORFs found in order of decreasing length. Use orfs() method to retrieve them.
		
		# flags to manipulate orf-finding behaviour:
		ALLOW_5PRIME_PARTIALS=>0, #allow for lacking start codon in logest orf.
		ALLOW_3PRIME_PARTIALS=>0, #allow for lacking stop codon in longest orf.
		FORWARD_STRAND=>1, #default set to true (analyze forward strand)
		REVERSE_STRAND=>1, #default set to true.
		ALLOW_NON_MET_STARTS=>0 #allow for non-methionine start codons.
		};
    bless ($obj);
    return ($obj);
}

## can include partial orfs at end of sequence.
sub allow_partials {
    my $self = shift;
    $self->{ALLOW_5PRIME_PARTIALS} = 1;
    $self->{ALLOW_3PRIME_PARTIALS} = 1;
}

sub allow_5prime_partials {
    my $self = shift;
    $self->{ALLOW_5PRIME_PARTIALS} = 1;
}

sub allow_3prime_partials {
    my $self = shift;
    $self->{ALLOW_3PRIME_PARTIALS} = 1;
}

sub forward_strand_only {
    my $self = shift;
    $self->{REVERSE_STRAND} = 0;
}

sub reverse_strand_only {
    my $self = shift;
    $self->{FORWARD_STRAND} = 0;
}

sub allow_non_met_starts {
    my $self = shift;
    $self->{ALLOW_NON_MET_STARTS} = 1;
}


sub get_longest_orf {
    my $self = shift;
    my $input_sequence = shift;
    
    unless (ref $self) {
	$self = new Longest_orf();
    }
    unless ($input_sequence) {
	print STDERR "I require a cDNA nucleotide sequence as my only parameter\n";
	return;
    }
    unless (length ($input_sequence) >= 3) {
	print STDERR "Sequence must code for at least a codon. Your seq_length is too short\n";
	return;
    }
    my @orfList = $self->capture_all_ORFs($input_sequence);
    $self->{all_ORFS} = \@orfList;
    return ($self);
}


sub capture_all_ORFs {
    
    my $self = shift;
    my $input_sequence = shift;
    unless ($input_sequence) {
	print STDERR "I require a cDNA nucleotide sequence as my only parameter\n";
	return;
    }
    unless (length ($input_sequence) >= 3) {
	print STDERR "Sequence must code for at least a codon. Your seq_length is too short\n";
	return;
    }
    
    $input_sequence = lc ($input_sequence);
    
    my (@starts, @stops, @orfs);
    if ($self->{FORWARD_STRAND}) {
	## analyse forward position
	@stops = $self->identify_putative_stops($input_sequence);
	@starts = $self->identify_putative_starts($input_sequence,\@stops);
	@orfs = &get_orfs (\@starts, \@stops, $input_sequence, '+');
    }
    
    if ($self->{REVERSE_STRAND}) {
	## reverse complement sequence and do again
	$input_sequence = &revcomp ($input_sequence);
	@stops = $self->identify_putative_stops($input_sequence);
	@starts = $self->identify_putative_starts($input_sequence, \@stops);
	push (@orfs,  &get_orfs (\@starts, \@stops, $input_sequence, '-'));
    }

    if (@orfs) {
	## set in order of decreasing length
	@orfs = reverse sort {$a->{length} <=> $b->{length}} @orfs;
	
	my $longest_orf = $orfs[0];
	my $start = $longest_orf->{start};
	my $stop = $longest_orf->{stop};
	my $seq = $longest_orf->{sequence};
	my $length = length($seq);
	my $protein = &get_protein ($seq);
	$self->{end5} = $start;  ## now coord is seq_based instead of array based.
	$self->{end3} = $stop;
	$self->{length} = $length;
	$self->{nt_seq} = $seq;
	$self->{pep_seq} = $protein;
	return (@orfs);
    }
}

sub orfs {
    my $self = shift;
    return (@{$self->{all_ORFS}});
}

#####################
# supporting methods
#####################

sub get_end5_end3 {
    my $self = shift;
    return ($self->{end5}, $self->{end3});
}

sub get_peptide_sequence {
    my $self = shift;
    return ($self->{pep_seq});
}

sub get_nucleotide_sequence {
    my $self = shift;
    return ($self->{nt_seq});
}


sub toString {
    my $self = shift;
    my ($end5, $end3) = $self->get_end5_end3();
    my $protein = $self->get_peptide_sequence();
    my $nt_seq = $self->get_nucleotide_sequence();
    my $ret_string = "Coords: $end5, $end3\n" 
	. "Protein: $protein\n"
	    . "Nucleotides: $nt_seq\n";
    return ($ret_string);
}


#################################

#Private methods:


sub get_orfs {
    my ($starts_ref, $stops_ref, $seq, $direction) = @_;
    my %last_delete_pos = ( 0=>-1,
			    1=>-1,
			    2=>-1); #store position of last chosen stop codon in spec reading frame.
    my @orfs;
    my $seq_length = length ($seq);
    foreach my $start_pos (@{$starts_ref}) {
	my $start_pos_frame = $start_pos % 3;
	foreach my $stop_pos (@{$stops_ref}) {
	    if ( ($stop_pos > $start_pos)   && #end3 > end5
		 ( ($stop_pos - $start_pos) % 3 == 0) #must be in-frame
		 && ($start_pos > $last_delete_pos{$start_pos_frame})) #only count each stop once.
	    {
		$last_delete_pos{$start_pos_frame} = $stop_pos;
		my ($start_pos_adj, $stop_pos_adj) = ( ($start_pos+1), ($stop_pos+1+2));
		#print "Startposadj: $start_pos_adj\tStopPosadj: $stop_pos_adj\n";
		# sequence based position rather than array-based
		
		my ($start, $stop) = ($direction eq '+') ? ($start_pos_adj, $stop_pos_adj) 
		    : (&revcomp_coord($start_pos_adj, $seq_length), &revcomp_coord($stop_pos_adj, $seq_length));

		#print "Start: $start\tStop: $stop\n";
		my $orfSeq =  substr ($seq, $start_pos, ($stop_pos - $start_pos + 3)); #include the stop codon too.
		my $protein = &get_protein($orfSeq);
		my $orf = { sequence => $orfSeq,
			    protein => $protein,
			    start=>$start,
			    stop=>$stop,
			    length=>length($orfSeq),
			    orient=>$direction
			    };
		push (@orfs, $orf);
		last;
	    }
	}
    }
    return (@orfs);
}


sub identify_putative_starts {
    my ($self, $seq, $stops_aref) = @_;
    my %starts;
    my %stops;
    foreach my $stop (@$stops_aref) {
	$stops{$stop} = 1;
    }

    if ($self->{ALLOW_5PRIME_PARTIALS} || $self->{ALLOW_NON_MET_STARTS}) {
	$starts{1} = 1 unless $stops{1};
	$starts{2} = 1 unless $stops{2};
	$starts{3} = 1 unless $stops{3};
    }
    
    if (!$self->{ALLOW_NON_MET_STARTS}) { #Look for ATG start codons.
	my $start_pos = index ($seq, "atg");
	while ($start_pos != -1) {
	    $starts{$start_pos} = 1;
	    #print "Start: $start_pos\n";
	    $start_pos = index ($seq, "atg", ($start_pos + 1));
	}
    } else {
	# find all residues just subsequent to a stop codon, in-frame:
	foreach my $stop (@$stops_aref) {
	    my $candidate_non_met_start = $stop +3;
	    unless ($stops{$candidate_non_met_start}) {
		$starts{$candidate_non_met_start} = 1;
	    }
	}
    }
    my @starts = sort {$a<=>$b} keys %starts;
    return (@starts);
}


sub identify_putative_stops {
    my ($self, $seq) = @_;
    my %stops;
    if ($self->{ALLOW_3PRIME_PARTIALS}) {
	## count terminal 3 nts as possible ORF terminators.
	my $seq_length = length ($seq);
	$stops{$seq_length} = 1;
	$seq_length--;
	$stops{$seq_length} = 1;
	$seq_length--;
	$stops{$seq_length} = 1;
    }
    foreach my $stop_codon ("taa", "tga", "tag") {
	my $stop_pos = index ($seq, $stop_codon);
	while ($stop_pos != -1) {
	    $stops{$stop_pos} = 1;
	    #print "Stop: $stop_pos\n";
	    $stop_pos = index ($seq, $stop_codon, ($stop_pos + 1)); #include the stop codon too.
	}
    }
    my @stops = sort {$a<=>$b} keys %stops;
    return (@stops);
}


sub revcomp {
    my ($seq) = @_;
    my $reversed_seq = reverse ($seq);
    $reversed_seq =~ tr/ACGTacgtyrkm/TGCAtgcarymk/;
    return ($reversed_seq);
}


sub revcomp_coord {
    my ($coord, $seq_length) = @_;
    return ($seq_length - $coord + 1);
}




   
1;


