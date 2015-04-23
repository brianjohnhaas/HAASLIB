#!/usr/local/bin/perl

package main;

our $SEE;

package Mummer_gene_converter;

use Mummer_coord_converter;
use Gene_obj;
use strict;

=description
    Given a Gene_obj object and a Mummer_coord_converter object, it will return a new Gene_obj with coordinates based on the new assembly ONLY IF THE ENTIRE GENE CAN BE TRANSPOSED.  The returned Gene_obj is not guaranteed to validate.  Gene validation should be performed separately.

NO NEED TO INSTANTIATE.

=cut

{#end of description
     ;
}
    
    
####
sub transpose_gene {
    my ($converter, $gene) = @_;
    unless (ref ($converter) && ref ($gene)) {
	die "I require a Mummer_coord_converter object and a Gene_obj object!\n";
    }
    my $new_gene = $gene->clone_gene();
    $new_gene->delete_isoforms();
    $new_gene->erase_gene_structure();
    
    my @exons = $gene->get_exons();
    my $gene_empty = 1;
    foreach my $exon (@exons) {
	my @exon_coords = $exon->get_mRNA_exon_end5_end3();
	my $exon_length = abs ($exon_coords[1] - $exon_coords[0]) + 1;
	my $transformed_exon_end5 = $converter->transform_coordinate_1_to_2 ($exon_coords[0]);
	my $transformed_exon_end3 = $converter->transform_coordinate_1_to_2 ($exon_coords[1]);
	
	my $transformed_exon_length = abs($transformed_exon_end5 - $transformed_exon_end3) + 1;
	print "Exon: ($exon_coords[0], $exon_coords[1]) ----> ($transformed_exon_end5, $transformed_exon_end3)\n" if $SEE;
	
	unless ($transformed_exon_end5 > 0 && $transformed_exon_end3 > 0 && $exon_length == $transformed_exon_length) {
	    next;
	}
	my $transformed_exon = new mRNA_exon_obj ($transformed_exon_end5, $transformed_exon_end3);
	my $cds_exon = $exon->get_CDS_obj();
	if ($cds_exon) {
	    my @cds_coords = $cds_exon->get_CDS_end5_end3();
	    my $cds_length = abs($cds_coords[1] - $cds_coords[0]) + 1;
	    my $transformed_cds_end5 = $converter->transform_coordinate_1_to_2 ($cds_coords[0]);
	    my $transformed_cds_end3 = $converter->transform_coordinate_1_to_2 ($cds_coords[1]);
	    my $transformed_cds_length = abs ($transformed_cds_end5 - $transformed_cds_end3) + 1;
	    print "CDS: ($cds_coords[0], $cds_coords[1]) ----> ($transformed_cds_end5, $transformed_cds_end3)\n" if $SEE;
	    if ($transformed_cds_end5 > 0 && $transformed_cds_end3 > 0 && $cds_length == $transformed_cds_length) {
		my $transformed_cds = new CDS_exon_obj($transformed_cds_end5, $transformed_cds_end3);
		$transformed_exon->set_CDS_exon_obj($transformed_cds);
	    }
	}
	$new_gene->add_mRNA_exon_obj($transformed_exon);
	$gene_empty = 0;
	
    }
    if ((!$gene_empty) && ($gene->number_of_exons() == $new_gene->number_of_exons())) { #avoid adding more partials for now.
	$new_gene->refine_gene_object();
	if ($SEE) {
	    print "Incoming Gene *************\n" . $gene->toString() . "############# Converted gene\n";
	    print $new_gene->toString();
	}
	
	if (my @isoforms = $gene->get_additional_isoforms()) {
	    foreach my $isoform (@isoforms) {
		my $new_isoform = transpose_gene($converter, $isoform);
		unless ($new_isoform) {
		    print "Error. Couldn't transpose isoform " . $isoform->toString();
		    return (undef);
		}
		$new_gene->add_isoform($new_isoform);
	    }
	    $new_gene->refine_gene_object();
	}
	return ($new_gene);
    } else {
	return (undef());
    }
}

1;  #end of module.

	
