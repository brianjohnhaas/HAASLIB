package main;

our $SEE;


package Sub_to_final_gene_obj_converter;

use strict;
use Sub_to_final_coord_transform;
use Gene_obj;

## No reason to instantiate.
## Use as a static class

sub convert_gene_coordinates {
    
    print "Converting gene coordinates.\n" if $SEE;
    
    my ($gene, $bac_asmbl_id, $chromo_asmbl_id, $converter, $convert_chromo_to_bac_flag) = @_;
    ## convert_chromo_to_bac must be set, otherwise, converts from bac to chromo coordinates.


    unless (ref ($converter) && ref ($gene)) {
        die "I require a gene object asn sub_to_final_coord_transform object";
    }
    my $new_gene = $gene->clone_gene();
    $new_gene->erase_gene_structure();
    my @exons = $gene->get_exons();
    my $gene_empty = 1;
    foreach my $exon (@exons) {
        my @exon_coords = $exon->get_mRNA_exon_end5_end3();
        
	my ($transformed_exon_end5, $transformed_exon_end3);

	if ($convert_chromo_to_bac_flag) {
	    my ($transformed_exon_end5_struct, $transformed_exon_end3_struct) = $converter->convert_chromo_to_bac($chromo_asmbl_id, @exon_coords);
	    if ($transformed_exon_end5_struct->{bac_asmbl_id} == $bac_asmbl_id) {
		$transformed_exon_end5 = $transformed_exon_end5_struct->{coord};
	    }
	    if ($transformed_exon_end3_struct->{bac_asmbl_id} == $bac_asmbl_id) {
		$transformed_exon_end3 = $transformed_exon_end3_struct->{coord};
	    }

	} else {
	    $transformed_exon_end5 = $converter->convert_bac_to_chromo ($bac_asmbl_id, $chromo_asmbl_id, $exon_coords[0]);
	    $transformed_exon_end3 = $converter->convert_bac_to_chromo($bac_asmbl_id, $chromo_asmbl_id, $exon_coords[1]);
	}
	
	print "$exon_coords[0] --> $transformed_exon_end5\n" 
	    . "$exon_coords[1] --> $transformed_exon_end3\n" if $SEE;
	
	unless ($transformed_exon_end5 > 0 && $transformed_exon_end3 > 0) {
            next;
        }
        my $transformed_exon = new mRNA_exon_obj ($transformed_exon_end5, $transformed_exon_end3);
        my $cds_exon = $exon->get_CDS_obj();
        if ($cds_exon) {
            my @cds_coords = $cds_exon->get_CDS_end5_end3();
           
	    my ($transformed_cds_end5, $transformed_cds_end3);
	    
	    if ($convert_chromo_to_bac_flag) {
		my ($transformed_cds_end5_struct, $transformed_cds_end3_struct) = $converter->convert_chromo_to_bac($chromo_asmbl_id, @cds_coords);
		if ($transformed_cds_end5_struct->{bac_asmbl_id} == $bac_asmbl_id) {
		    $transformed_cds_end5 = $transformed_cds_end5_struct->{coord};
		}
		if ($transformed_cds_end3_struct->{bac_asmbl_id} == $bac_asmbl_id) {
		    $transformed_cds_end3 = $transformed_cds_end3_struct->{coord};
		}
		
	    } else {
		$transformed_cds_end5 = $converter->convert_bac_to_chromo ($bac_asmbl_id, $chromo_asmbl_id, $cds_coords[0]);
		$transformed_cds_end3 = $converter->convert_bac_to_chromo ($bac_asmbl_id, $chromo_asmbl_id, $cds_coords[1]);
	    }

	    if ($transformed_cds_end5 > 0 && $transformed_cds_end3 > 0) {
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
                my $new_isoform = convert_gene_coordinates($isoform, $bac_asmbl_id, $chromo_asmbl_id, $converter, $convert_chromo_to_bac_flag);
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


