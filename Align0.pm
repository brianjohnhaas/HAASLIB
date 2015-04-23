#!/usr/bin/env perl

package Align0;

use strict;
use warnings;
use Carp;

sub new {
	my $packagename = shift;
	
	my ($seqA, $seqB) = @_;

	unless ($seqA && $seqB) {
		confess "Error, need two sequence strings as parameters";
	}

	my $self = {
		
		seqA => $seqA,
		seqB => $seqB,

		alignA => undef,
		alignB => undef,
		
		num_matches => undef,
		
		num_mismatches => undef,
		mismatch_positions => undef, # becomes arrayref

		num_indels => undef,
		indel_positions => undef, # becomes arrayref

        left_bound => undef,
        right_bound => undef,

		
	};
	
	bless ($self, $packagename);
	
	$self->_run_align0($seqA, $seqB);


	return($self);
}

####
sub get_seqA_align_positions {
    my $self = shift;
    
    my $left_bound = $self->{left_bound};
    my $right_bound = $self->{right_bound};
    
    my @seqA_align_chars = split(//, $self->{alignA});

    my $char_counter = 0;
    my $pos_counter = 0;
    
    my ($left_pos, $right_pos);
    
    foreach my $char (@seqA_align_chars) {
        $char_counter++;
        if ($char =~ /\w/) {
            $pos_counter++;
        }

        if ($char_counter == $left_bound+1) {
            $left_pos = $pos_counter;
        }

        if ($char_counter == $right_bound+1) {
            $right_pos = $pos_counter;
            last;
        }
    }

    return($left_pos, $right_pos);
}



####
sub _run_align0 {
	my $self = shift;
	my ($seqA, $seqB) = @_;
	
	my $tmp_fileA = "$$.tmp.A";
	my $tmp_fileB = "$$.tmp.B";
	my $tmp_output = "$$.tmp.out";
	

	&_write_tmp_fasta("seqA", $seqA, $tmp_fileA);
	&_write_tmp_fasta("seqB", $seqB, $tmp_fileB);
	
	my $cmd = "align0 $tmp_fileA $tmp_fileB > $tmp_output 2>/dev/null";
	system($cmd);
	
	if (! -s $tmp_output) {
		confess "Error, cmd $cmd died with ret $?";
	}
	
    # system("cat $tmp_output");
	
	## parse output file
	my ($alignA, $alignB) = &_parse_alignment_output($tmp_output);

	if (length($alignA) != length($alignB)) {
		confess "Error, alignment lengths nonidentical:\nSeqA: $alignA\nSeqB:$alignB\n\n";
	}

	$self->{alignA} = $alignA;
	$self->{alignB} = $alignB;
	
	unlink($tmp_fileA, $tmp_fileB, $tmp_output);
	
	my ($match, $mismatch, $indel) = $self->_count_aligned_positions($alignA, $alignB);
	
	#print "$match\t$mismatch\t$indel\n";

	
	return;
	
}


####
sub _write_tmp_fasta {
	my ($acc, $seq, $filename) = @_;

	open (my $ofh, ">$filename") or die "Error, cannot write file $filename";
	print $ofh ">$acc\n$seq\n";
	close $ofh;

	return;
}


####
sub _parse_alignment_output {
	my ($file) = @_;

	my ($alignA, $alignB);

	open (my $fh, $file) or die "Error, cannot read file $file";
	while (<$fh>) {
		chomp;
		if (/^seq[AB] /) {
			my $seq_acc = $1;
			my @x = split (/\t/);
			my ($acc, $align) = split(/\s+/);
			if ($align =~ /[^A-Za-z\.\-\s]/) {
				next; # not an alignment string
			}
			
			
			if ($acc eq "seqA") {
				$alignA .= $align;
			}
			elsif ($acc eq "seqB") {
				$alignB .= $align;
			}
			else {
				confess "Error, parsed alignment line but not seqA or seqB: $_";
			}
		}
	}
	close $fh;
	
	$alignA =~ s/\s//g;
	$alignB =~ s/\s//g;

	return($alignA, $alignB);
}

####
sub _count_aligned_positions {
	my $self = shift;
	my ($alignA, $alignB) = @_;

	my @A_chars = split(//, uc $alignA);
	my @B_chars = split(//, uc $alignB);

	## end gaps don't count
	
	my $left_align_index = &_define_left_align_index(\@A_chars, \@B_chars);
	my $right_align_index = &_define_right_align_index(\@A_chars, \@B_chars);

    $self->{left_bound} = $left_align_index;
    $self->{right_bound} = $right_align_index;
    

	my @mismatch_positions;
	my @indel_positions;
		
	my ($num_match, $num_mismatch, $num_indel) = (0,0,0);

	for (my $i = $left_align_index; $i <= $right_align_index; $i++) {
		
		my $charA = $A_chars[$i];
		my $charB = $B_chars[$i];
		
		if ($charA eq "N" || $charB eq "N") {
			# N's don't count.
			next;
		}
		elsif ($charA eq $charB) {
			$num_match++;
		}
		elsif ($charA =~ /[\.\-]/ || $charB =~ /[\.\-]/) {
			$num_indel++;
			push (@indel_positions, $i);
		}
		else {
			$num_mismatch++; # only option left.
			push (@mismatch_positions, $i);
		}
	}

	$self->{num_matches} = $num_match;
	$self->{num_mismatches} = $num_mismatch;
	$self->{mismatch_positions} = \@mismatch_positions;
	$self->{num_indels} = $num_indel;
	$self->{indel_positions} = \@indel_positions;
	
	return($num_match, $num_mismatch, $num_indel);
	
}


####
sub _define_left_align_index {
	my ($A_chars_aref, $B_chars_aref) = @_;

	my $foundA = 0;
	my $foundB = 0;

	for (my $i = 0; $i <= $#$A_chars_aref; $i++) {
		my $charA = $A_chars_aref->[$i];
		my $charB = $B_chars_aref->[$i];
		
		if ($charA =~ /\w/) {
			$foundA++;
		}
		if ($charB =~ /\w/) {
			$foundB++;
		}

		if ($foundA && $foundB) {
			return($i);
		}
	}

	confess "Error, no left bound found in alignment.  Shouldn't happen";
}

####
sub _define_right_align_index {
	my ($A_chars_aref, $B_chars_aref) = @_;

	my $foundA = 0;
	my $foundB = 0;
	
	for (my $i = $#$A_chars_aref; $i >= 0; $i--) {
		
		my $charA = $A_chars_aref->[$i];
		my $charB = $B_chars_aref->[$i];
		
		if ($charA =~ /\w/) {
			$foundA++;
		}
		if ($charB =~ /\w/) {
			$foundB++;
		}

		if ($foundA && $foundB) {
			return($i);
		}
	}

	confess "Error, no left bound found in alignment.  Shouldn't happen";
}



1; #EOM
