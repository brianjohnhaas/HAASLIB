

## No reason to instantiate this class.  Use as static class.

package Gene_obj_splitter;
use strict;
use Gene_obj;


## cutCoord ends up in the left_half of the gene regardless of orientation.

sub split_gene {
    my ($gene_obj, $cutCoord) = @_;

    print "## Splitting gene model: " . $gene_obj->{Model_feat_name} . "\n";
    
    my $gene_orientation = $gene_obj->get_orientation();
    
    my ($lend, $rend) = sort {$a<=>$b} $gene_obj->get_coords();
        
    my $left_half_gene_obj = $gene_obj->clone_gene();
    $left_half_gene_obj->delete_isoforms();
    $left_half_gene_obj->erase_gene_structure();

    my $right_half_gene_obj = $gene_obj->clone_gene();
    $right_half_gene_obj->delete_isoforms();
    $right_half_gene_obj->erase_gene_structure();

    my @exons = $gene_obj->get_exons();
    
    foreach my $exon (@exons) {
	my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();
	if ($rend <= $cutCoord) {
	    $left_half_gene_obj->add_mRNA_exon_obj($exon->clone_exon());
	} elsif ($lend > $cutCoord) {
	    $right_half_gene_obj->add_mRNA_exon_obj($exon->clone_exon());
	} else {
	    ## cutCoord must be within exon coordinate range:
	    unless ($cutCoord >= $lend && $cutCoord <= $rend) {
		die "ERROR, cutCoord ($cutCoord) should be within exon [$lend-$rend] but isn't.!\n";
	    }
	    
	    my ($leftExon, $rightExon) = &split_exon($exon, $cutCoord, $gene_orientation);
	    $left_half_gene_obj->add_mRNA_exon_obj($leftExon);
	    $right_half_gene_obj->add_mRNA_exon_obj($rightExon);
	}
    }

    if (my @isoforms = $gene_obj->get_additional_isoforms()) {
	foreach my $isoform (@isoforms) {
	    my ($leftGene, $rightGene) = split_gene($isoform, $cutCoord);
	    if (ref $leftGene) {
		$left_half_gene_obj->add_isoform($leftGene);
	    }
	    if (ref $rightGene) {
		$right_half_gene_obj->add_isoform($rightGene);
	    }
	}
    }

    
    #### Check gene obj content:
    if ( scalar($left_half_gene_obj->get_exons()) == 0) { #no exons
	# see if isoforms exist:
	print "left_half_gene_obj has no exons.\n";
	if (scalar($left_half_gene_obj->get_additional_isoforms()) == 0) {
	    undef($left_half_gene_obj); #nothing to report
	    print "\tno isoforms either.  Undef'ing.\n";
	} else {
	    # use an isoform to template the gene:
	    print "\tgot isoforms.  Templating an isoform.\n";
	    my @isoforms = $left_half_gene_obj->get_additional_isoforms();
	    my $templateIsoform = shift @isoforms;
	    foreach my $iso (@isoforms) {
		$templateIsoform->add_isoform($iso);
	    }
	    $left_half_gene_obj = $templateIsoform;
	}
    }

    

    ## do the same for right_gene_obj:
    if ( scalar($right_half_gene_obj->get_exons()) == 0) { #no exons
	# see if isoforms exist:
	print "right_half_gene_obj has no exons.\n";
	if (scalar($right_half_gene_obj->get_additional_isoforms()) == 0) {
	    undef($right_half_gene_obj); #nothing to report
	    print"\tno isoforms either. Undef'ing.\n";
	} else {
	    # use an isoform to template the gene:
	    print "\tgot isoforms.  Templating an isoform.\n";
	    my @isoforms = $right_half_gene_obj->get_additional_isoforms();
	    my $templateIsoform = shift @isoforms;
	    foreach my $iso (@isoforms) {
		$templateIsoform->add_isoform($iso);
	    }
	    $right_half_gene_obj = $templateIsoform;
	}
    }
    
    
    return ($left_half_gene_obj, $right_half_gene_obj);
}




## private
sub split_exon {
    my ($exon, $cutCoord, $gene_orientation) = @_;
    
    my ($end5, $end3) = $exon->get_coords();
    my ($exon_lend, $exon_rend) = sort {$a<=>$b} ($end5, $end3); 
    unless ($cutCoord >= $exon_lend && $cutCoord < $exon_rend) {
	die "$cutCoord not within range of exon: $exon_lend-$exon_rend\n";
    }
    
    my ($left_exon, $right_exon) = ($exon->clone_exon(), $exon->clone_exon());
    
    ## extract cds info:
    my ($cds, $cds_end5, $cds_end3, $cds_lend, $cds_rend);
    if ($cds = $exon->get_CDS_obj()) {
	($cds_end5, $cds_end3) = $cds->get_coords();
	($cds_lend, $cds_rend) = sort {$a<=>$b} ($cds_end5, $cds_end3);
    }


    ## process left exon:
    my ($left_end5, $left_end3);
    if ($gene_orientation eq "+") {
	($left_end5, $left_end3) = ($end5, $cutCoord);
    } else {
	($left_end5, $left_end3) = ($cutCoord, $end3);
    }
    $left_exon->set_coords($left_end5, $left_end3);
   
    if ($cds) {
	if ($cutCoord >= $cds_lend && $cutCoord < $cds_rend) {
	    # cutCoord contained within cds.
	    
	    my ($left_cds_end5, $left_cds_end3);
	    if ($gene_orientation eq "+") {
		($left_cds_end5, $left_cds_end3) = ($cds_end5, $cutCoord);
	    } else {
		($left_cds_end5, $left_cds_end3) = ($cutCoord, $cds_end3);
	    }
	    $left_exon->get_CDS_obj()->set_coords($left_cds_end5, $left_cds_end3);
	} else {
	    ## not with cds.  Delete the cds from this exon:
	    $left_exon->remove_CDS_exon();
	}
    }
    

    ## Process right exon:
    my $posAfterCut = $cutCoord + 1;
    my ($right_end5, $right_end3);
    if ($gene_orientation eq "+") {
	($right_end5, $right_end3) = ($posAfterCut, $end3);
    } else {
	($right_end5, $right_end3) = ($end5, $posAfterCut);
    }
    $right_exon->set_coords($right_end5, $right_end3);
   
    if ($cds) {
	if ($posAfterCut >= $cds_lend && $posAfterCut <= $cds_rend) {
	    my ($right_cds_end5, $right_cds_end3);
	    if ($gene_orientation eq "+") {
		($right_cds_end5, $right_cds_end3) = ($posAfterCut, $cds_end3);
	    } else {
		($right_cds_end5, $right_cds_end3) = ($cds_end5, $posAfterCut);
	    }
	    $right_exon->get_CDS_obj()->set_coords($right_cds_end5, $right_cds_end3);
	} else {
	    # not within cds. Delete it:
	    $right_exon->remove_CDS_exon();
	}
    }

    return ($left_exon, $right_exon);
}

    
1; #EOM




	    

