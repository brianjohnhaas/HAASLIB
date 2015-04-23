#!/usr/local/bin/perl

package main;
our $SEE;


package Alternative_splice_comparer;
use Gene_obj;
use strict;

sub new {
    my $packagename = shift;
    my $self = {};
    bless ($self, $packagename);
    return ($self);
}


####
sub compare_isoforms {
    my $self = shift;
    my ($gene1, $gene2) = @_;
    ## Looking for:
    #    -unspliced introns
    #    -conventional alt-splice isoforms
    #    -transcriptional start or polyadenylation site within intron
    #    -exon skipping
    #    -alternate exons

    ## Look for Unspliced Introns
    my $struct  = { unspliced_introns => 0,
		    conventional_alt_splice => 0,
		    start_or_end_within_intron => 0,
		    exon_skipping => 0,
		    alternate_exons => 0 };
    
    
    my @unspliced_introns = ($self->find_unspliced_introns($gene1, $gene2), $self->find_unspliced_introns($gene2, $gene1));
    if (@unspliced_introns) {
	print "*** Unspliced introns \n";
	$struct->{unspliced_introns} = 1;
    }

    ## Look for the conventional alt splice isoforms (diff donors/acceptors for introns).
    my @alternate_acceptors_n_donors = $self->find_conventional_alt_splice_isoforms($gene1, $gene2);
    if (@alternate_acceptors_n_donors) {
	print "*** Conventional Alt splice (donor and/or acceptor)\n";
	$struct->{conventional_alt_splice} = 1;
    }
    

    my @intron_starts_or_ends = ($self->find_starts_and_ends_within_introns ($gene1, $gene2), $self->find_starts_and_ends_within_introns ($gene2, $gene1));
    if (@intron_starts_or_ends) {
	print "*** Transcriptional start or polyadenylation site within an intron.\n";
	$struct->{start_or_end_within_intron} = 1;
    }

    ## Look for exon skipping events.
    my @exon_skips = ($self->find_exon_skipping_events($gene1, $gene2), $self->find_exon_skipping_events($gene2, $gene1));
    if (@exon_skips) {
	print "*** Exon skipping event detected.\n";
	$struct->{exon_skipping} = 1;
    }

    ## Look or Alternate exons
    my @alternate_exons = $self->find_alternate_exons($gene1, $gene2);
    if (@alternate_exons) {
	print "*** Found alternate exons\n";
	$struct->{alternate_exons} = 1;
    }
    
    return ($struct);

}


# private
sub enumerate_exons_of_gene {
    my $gene_obj = shift;
    # put everything in forward coordinate axis:
    my %exon_coords;
    my @exons = $gene_obj->get_exons();
    foreach my $exon (@exons) {
	my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();
	$exon_coords{$lend} = $rend;
    }
    return (%exon_coords);
}


#private
####
sub enumerate_introns_of_gene {
    my $gene_obj = shift;
    ## Put everything in forward strand coordinate axis.
    my %introns;
    my @exons = sort {$a->{end5}<=>$b->{end5}} $gene_obj->get_exons();
    for (my $i = 0; $i < $#exons; $i++) {
	my ($exon1_lend, $exon1_rend) = sort {$a<=>$b} $exons[$i]->get_coords();
	my ($exon2_lend, $exon2_rend) = sort {$a<=>$b} $exons[$i+1]->get_coords();
	my ($intron_end5, $intron_end3) = ($exon1_rend + 1, $exon2_lend - 1);
	$introns{$intron_end5} = $intron_end3;
    }
    return (%introns);
}


=over 4

=item find_unspliced_introns()

B<Description:> Find unspliced introns in gene_1 when compared to gene_2

B<Parameters:> $gene1, $gene2

B<Returns:> @unspliced_introns

@unspliced_introns is a list of coordinate pairs representing the unspliced introns found in gene 1 when compared to gene2

@unspliced_introns = ([end5,end3], ...)

=back

=cut

####
sub find_unspliced_introns {
    my $self = shift;
    my ($gene1, $gene2) = @_;
    ## Look for unspliced intron found in gene1 when compared to gene2
    my %gene1_intron_coords = &enumerate_introns_of_gene ($gene1);
    my %gene2_exon_coords = &enumerate_exons_of_gene($gene2);
    
    my @unspliced_introns;
    foreach my $intron_lend (keys %gene1_intron_coords) {
	my $intron_rend = $gene1_intron_coords{$intron_lend};
	
	foreach my $exon_lend (keys %gene2_exon_coords) {
	    my $exon_rend = $gene2_exon_coords{$exon_lend};
	    
	    if ($intron_lend > $exon_lend && $intron_rend < $exon_rend) { #unspliced intron found
		push (@unspliced_introns, [$intron_lend, $intron_rend]);
	    }
	}
    }
    return (@unspliced_introns);
}



=over 4

=item find_conventional_alt_splice_isoforms()

B<Description:> Looks for different donor and acceptor sites within overlapping introns of genes

B<Parameters:> $gene1, $gene2

B<Returns:> @coords

@coords is a list of the intron coordinates which differ in gene 1 when compared to gene 2

@coords = (intron_lend, ...)

=back

=cut


####
sub find_conventional_alt_splice_isoforms {
    my $self = shift;
    my ($gene1, $gene2) = @_;
    print "## Looking for conventional alt splice isoforms (diff donors, acceptors)\n" if $SEE;
    my %introns_gene1 = &enumerate_introns_of_gene ($gene1);
    my %introns_gene2 = &enumerate_introns_of_gene ($gene2);
    
    # algorithm
    #   -intron must match one and only other intron and have diff coordinates

    my @alternate_acceptors_n_donors; # holds coordinates for all gene1 diff boundaries.
    foreach my $intron1_lend (keys %introns_gene1) {
	my $intron1_rend = $introns_gene1{$intron1_lend};
	my $num_matches = 0;
	my ($match_lend, $match_rend);
	foreach my $intron2_lend (keys %introns_gene2) {
	    my $intron2_rend = $introns_gene2{$intron2_lend};
	    print "comparing ($intron1_lend, $intron1_rend) to ($intron2_lend, $intron2_rend)\n" if $SEE;
	    if ($intron1_lend < $intron2_rend && $intron1_rend > $intron2_lend) {
		$num_matches++;
		print "match: $num_matches\n" if $SEE;
		($match_lend, $match_rend) = ($intron2_lend, $intron2_rend);
	    }
	}
	if ($num_matches == 1) {
	    ## Make sure no exon is found within either intron
	    if ($match_lend != $intron1_lend || $match_rend != $intron1_rend) { #intersting boudnary diff.
		my $exon_incapsulated = 0;
		foreach my $exon ($gene1->get_exons(), $gene2->get_exons()) {
		    my ($exon_lend, $exon_rend) = sort {$a<=>$b} $exon->get_coords();
		    if ( ($exon_lend > $intron1_lend && $exon_rend < $intron1_rend) ||
			 ($exon_lend > $match_lend && $exon_rend < $match_rend)) { #exon incapsulated within one of the introns
			$exon_incapsulated = 1;
			print "Exon incapsulated in the suspect intron\n" if $SEE;
			last;
		    }
		}
		unless ($exon_incapsulated) {
		    print "entering match data.\n" if $SEE;
		    if ($match_lend != $intron1_lend) {
			push (@alternate_acceptors_n_donors, $intron1_lend);
			print "diff lend ($match_lend vs. $intron1_lend)\n" if $SEE;
		    }
		    if ($match_rend != $intron1_rend) {
			push (@alternate_acceptors_n_donors, $intron1_rend);
			print "diff rend ($match_rend vs. $intron1_rend)\n" if $SEE;
		    }
		}
	    }
	}
    }
    return (@alternate_acceptors_n_donors);
}




=over 4

=item find_exon_skipping_events()

B<Description:> Finds an exon of gene_1 which reside within an intron of gene_2

B<Parameters:> gene1, gene2

B<Returns:> @skipped_exons

@skipped_exons = ([exon_lend,exon_rend], ...) 

=back

=cut


####
sub find_exon_skipping_events {
    my $self = shift;
    my ($gene1, $gene2) = @_;
    
    # Algorithm:
    #   -find an internal exon of gene 1 that resides within an intron of gene 2.  Flanking exons must be anchored to the other isoform
    my %gene1_exons = &enumerate_exons_of_gene($gene1);
    my %gene2_introns = &enumerate_introns_of_gene($gene2);
    
    my @potential_skipped_exons;
    foreach my $exon1_lend (keys %gene1_exons) {
	my $exon1_rend = $gene1_exons{$exon1_lend};
	
	## See if within intron of second gene
	foreach my $intron2_lend (keys %gene2_introns) {
	    my $intron2_rend = $gene2_introns{$intron2_lend};
	    
	    if ($exon1_lend > $intron2_lend && $exon1_rend < $intron2_rend) { #exon incapsulated in intron
		push (@potential_skipped_exons, [$exon1_lend, $exon1_rend]);
	    }
	}
    }
    
    ## Verify flanking exons are anchorable:
    my @skipped_exons;
    if (@potential_skipped_exons) {
	foreach my $potential_skipped_exon (@potential_skipped_exons) {
	    my ($exon_lend, $exon_rend) = @$potential_skipped_exon;
	    
	    ## Try to anchor left exon
	    my $anchor_left_exon = 0;
	    foreach my $exon1 ($gene1->get_exons()) {
		my ($exon1_lend, $exon1_rend) = sort {$a<=>$b} $exon1->get_coords();
		unless ($exon1_rend < $exon_lend) { next;}
		foreach my $exon2 ($gene2->get_exons()) {
		    my ($exon2_lend, $exon2_rend) = sort {$a<=>$b} $exon2->get_coords();
		    unless ($exon2_rend < $exon_lend) { next;}
		    
		    if ($exon1_lend < $exon2_rend && $exon1_rend > $exon2_lend) { #anchorable
			$anchor_left_exon = 1;
			last;
		    }
		}
		if ($anchor_left_exon) { last;}
	    }
	    unless ($anchor_left_exon) { next;}
	    
	    ## Try to anchor the right exon
	    my $anchor_right_exon = 0;
	    foreach my $exon1 ($gene1->get_exons()) {
		my ($exon1_lend, $exon1_rend) = sort {$a<=>$b} $exon1->get_coords();
		unless ($exon1_lend > $exon_rend) { next;}
		foreach my $exon2 ($gene2->get_exons()) {
		    my ($exon2_lend, $exon2_rend) = sort {$a<=>$b} $exon2->get_coords();
		    unless ($exon2_lend > $exon_rend) { next;}
		    
		    if ($exon1_lend < $exon2_rend && $exon1_rend > $exon2_lend) { #anchorable
			$anchor_right_exon = 1;
			last;
		    }
		}
		if ($anchor_right_exon) { last;}
	    }
	    if ($anchor_right_exon && $anchor_left_exon) {
		push (@skipped_exons, $potential_skipped_exon);
	    }
	}
	
    }
    return (@skipped_exons);
}







=over 4

=item find_alternate_exons()

B<Description:> Finds terminal exons that are different and non-overlapping, and adjacent to overlapping exons.

B<Parameters:> $gene1, $gene2

B<Returns:> @alternate_exons

@alternate_exons = ([exon_lend,exon_rend], ...)



=back

=cut

sub find_alternate_exons {
    my $self = shift;
    my ($gene1, $gene2) = @_;
    my @alternate_exons; # store coords of alternate exons
    # Algorithm:
    #    -Looking at terminal exons, should have non-overlapping exons prior to the first overlapping exon
    
    ## Look from front to back:
    my @gene1_exons = sort {$a->{end5}<=>$b->{end5}} $gene1->get_exons();
    my @gene2_exons = sort {$a->{end5}<=>$b->{end5}} $gene2->get_exons();
    
    for (my $i = 0; $i <= $#gene1_exons; $i++) {
	my ($exon1_lend, $exon1_rend) = sort {$a<=>$b} $gene1_exons[$i]->get_coords();
	my $overlapping_j = undef();
	for (my $j = 0; $j <= $#gene2_exons; $j++) {
	    my ($exon2_lend, $exon2_rend) = sort {$a<=>$b} $gene2_exons[$j]->get_coords();
	    
	    ## check for overlap
	    if ($exon1_lend < $exon2_rend && $exon1_rend > $exon2_lend) {
		$overlapping_j = $j;
		last;
	    }
	}
	if (defined($overlapping_j)) {
	    ## See if i and j are not first:
	    if ($i != 0 && $overlapping_j != 0) {
		for (my $x=0; $x < $i; $x++) {
		    my $exon = $gene1_exons[$x];
		    my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();
		    push (@alternate_exons, [$lend,$rend]);
		}
		for (my $x = 0; $x < $overlapping_j; $x++) {
		    my $exon = $gene2_exons[$overlapping_j];
		    my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();
		    push (@alternate_exons, [$lend,$rend]);
		}
	    }
	    last;
	}
    }
    
    ## Look from back to front:
     for (my $i = $#gene1_exons; $i >= 0; $i--) {
	my ($exon1_lend, $exon1_rend) = sort {$a<=>$b} $gene1_exons[$i]->get_coords();
	my $overlapping_j = undef();
	for (my $j = $#gene2_exons; $j >= 0; $j--) {
	    my ($exon2_lend, $exon2_rend) = sort {$a<=>$b} $gene2_exons[$j]->get_coords();
	    
	    ## check for overlap
	    if ($exon1_lend < $exon2_rend && $exon1_rend > $exon2_lend) {
		$overlapping_j = $j;
		last;
	    }
	}
	if (defined($overlapping_j)) {
	    ## See if i and j are not last:
	    if ($i != $#gene1_exons && $overlapping_j != $#gene2_exons) {
		for (my $x=$#gene1_exons; $x > $i; $x--) {
		    my $exon = $gene1_exons[$x];
		    my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();
		    push (@alternate_exons, [$lend,$rend]);
		}
		for (my $x = $#gene2_exons; $x > $overlapping_j; $x--) {
		    my $exon = $gene2_exons[$overlapping_j];
		    my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();
		    push (@alternate_exons, [$lend,$rend]);
		}
	    }
	    last;
	}
    }

    return (@alternate_exons);
}


=over 4

=item find_starts_and_ends_within_introns()

B<Description:> The first and last exons of gene_1 are compared to the introns of gene_2. 

B<Parameters:> $gene1, $gene2

B<Returns:> @coords

@coords contains the coordinates of either the very end5 or very end3 of terminal exons which fall into introns of gene_2

=back

=cut


####
sub find_starts_and_ends_within_introns {
    my $self = shift;
    my ($gene1, $gene2) = @_;
    
    my @intronic_coords; #store coord of offending exon.
    # Algorithm:
    #   -first and last exon of gene1 is compared to introns of gene2
    my @gene1_exons = $gene1->get_exons();
    my %gene2_introns = &enumerate_introns_of_gene($gene2);
    my %gene2_exons = &enumerate_exons_of_gene($gene2);
    my $first_exon = $gene1_exons[0];
    my ($end5, $end3) = $first_exon->get_coords();
    foreach my $intron_lend (keys %gene2_introns) {
	my $intron_rend = $gene2_introns{$intron_lend};
	if ($end5 >= $intron_lend && $end5 <= $intron_rend) { #endpoint encapsulated by intron.
	    ## make sure exon overlaps another exon
	    foreach my $exon_lend (keys %gene2_exons) {
		my $exon_rend = $gene2_exons{$exon_lend};
		my ($lend, $rend) = sort {$a<=>$b} ($end5, $end3);
		if ($rend > $exon_lend && $lend < $exon_rend) { #overlap
		    push (@intronic_coords, $end5);
		    last;
		}
	    }
	    last;
	}
    }
    
    ## Now try last exon
    my $last_exon = $gene1_exons[$#gene1_exons];
    my ($end5, $end3) = $last_exon->get_coords();
    foreach my $intron_lend (keys %gene2_introns) {
	my $intron_rend = $gene2_introns{$intron_lend};
	if ($end3 >= $intron_lend && $end3 <= $intron_rend) { #endpoint encapsulated by intron.
	    ## Make sure exon overlaps another exon
	    foreach my $exon_lend (keys %gene2_exons) {
		my $exon_rend = $gene2_exons{$exon_lend};
		my ($lend, $rend) = sort {$a<=>$b} ($end5, $end3);
		if ($rend > $exon_lend && $lend < $exon_rend) {
		    push (@intronic_coords, $end3);
		    last;
		}
	    }
	    last;
	}
    }
    return (@intronic_coords);
}



=over 4

=item compare_exons()

B<Description:> Compares all CDS exons between genes 1 and 2, returns number of identical CDS exons and total number of CDS exons between the two genes.

B<Parameters:> $gene1, $gene2

B<Returns:> ($num_identical_CDS_exons, $num_total_CDS_exons)


=back

=cut



####
sub compare_exons {
    my $self = shift;
    my ($gene1, $gene2) = @_;
    print "gene1_strand: $gene1->{strand}\n";
    my $clone_1 = $gene1->clone_gene();
    print "clone1_strand: " . $clone_1->{strand} . "\n";
    my $clone_2 = $gene2->clone_gene();
    $clone_1->trim_UTRs();
    print "clone1_strand, utrs trimmed: " . $clone_1->{strand} . "\n";
    $clone_2->trim_UTRs();
   
    my @exons_1 = $clone_1->get_exons();
    
    my @exons_2 = $clone_2->get_exons();
    
    my @identity_list = ();
    my @all_exons = sort {$a->{end5}<=>$b->{end5}} (@exons_1, @exons_2);
    
    for (my $i=0; $i <= $#all_exons-1; $i++) {
	
	my $curr_exon = $all_exons[$i];
	my $next_exon = $all_exons[$i+1];
	
	my ($curr_exon_end5, $curr_exon_end3) = $curr_exon->get_coords();
	
	my ($next_exon_end5, $next_exon_end3) = $next_exon->get_coords();

	if ($curr_exon_end5 == $next_exon_end5 && $curr_exon_end3 == $next_exon_end3) {
	    $identity_list[$i] = 1;
	    $identity_list[$i+1] = 1;
	    $i++; #if A = B, then go onto comparing C to D, not B to C.
	}
    }

    my ($num_identical_exons, $total_num_exons) = (0,0);
    for (my $i=0; $i <= $#all_exons; $i++) {
	if ($identity_list[$i]) {
	    $num_identical_exons++;
	}
	$total_num_exons++;
    }
    
    return ($num_identical_exons, $total_num_exons);
}

	


=over 4

=item start_or_stop_within_intron()

B<Description:> Compares the start codon and stop codon position of gene1 to the introns of gene2.

B<Parameters:> $gene1, $gene2

B<Returns:> ($start_within_intron, $stop_within_intron)

return values are 0|1 meaning true|false for each return parameter.

=back

=cut


sub start_or_stop_codon_within_intron {
    my $self = shift;
    my ($gene1, $gene2) = @_;
 

    ## Look for annotated start codon or stop codon within intron:
    my ($annotated_start_within_intron, $annotated_stop_within_intron) = (0,0);

    my ($start_codon, $stop_codon) = $gene1->get_model_span();
    my (@alignment_segments) = sort {$a->{end5}<=>$b->{end5}} $gene2->get_exons();
    if ($#alignment_segments > 0) { #multiple segments:
	for (my $i=1; $i <= $#alignment_segments; $i++) {
	    my $prev_seg = $alignment_segments[$i-1];
	    my ($prev_lend, $prev_rend) = sort {$a<=>$b} $prev_seg->get_coords();
	    my $curr_seg = $alignment_segments[$i];
	    my ($curr_lend, $curr_rend) = sort {$a<=>$b} $curr_seg->get_coords();
	    my ($intron_lend, $intron_rend) = ($prev_rend+1, $curr_lend-1);
	    if ($start_codon >= $intron_lend && $start_codon <= $intron_rend) {
		$annotated_start_within_intron = 1;
	    }
	    if ($stop_codon >= $intron_lend && $stop_codon <= $intron_rend) {
		$annotated_stop_within_intron = 1;
	    }
	}
    }
    return ($annotated_start_within_intron, $annotated_stop_within_intron);
}













1;
	
    
