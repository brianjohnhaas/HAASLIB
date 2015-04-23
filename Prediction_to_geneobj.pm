#!/usr/local/bin/perl

## There's no reason to create objects from this class.  Use methods as 
## static methods.  Simply returns Gene_obj's based on model feat_names.


package Prediction_to_geneobj;

use strict;
use DBI;
use lib ($ENV{EGC_SCRIPTS});
use Egc_library;
use Gene_obj;

# get gene based on model feat_name
sub get_gene {
    my ($dbproc, $feat_name) = @_;
    my $gene_obj = new Gene_obj();
     ## set up gene parameters:
    $gene_obj->{is_pseudogene} = 0;
    $gene_obj->{com_name} = "predicted protein";
    my $query = "select ev_type from phys_ev where feat_name = \"$feat_name\"\n";
    my $pred_derived = &first_result_sql ($dbproc, $query);
    if ($pred_derived) {
	$gene_obj->{pub_comment} = "predicted by $pred_derived";
    }
    
    ## build the gene structure
    my $query = "select a_exon.end5, a_exon.end3, a_exon.feat_name\n"
	. "from asm_feature a_exon, feat_link f_m_e\n"
	. "where f_m_e.parent_feat = \"$feat_name\" and f_m_e.child_feat = a_exon.feat_name\n"
	    . "and a_exon.feat_type = \"exon\"\n";
    
    my @results = &do_sql ($dbproc, $query);
    foreach my $result (@results) {
	my ($exon_end5, $exon_end3, $exon_feat_name) = split (/\t/, $result);
	my $mRNA_exon = new mRNA_exon_obj($exon_end5, $exon_end3);
	$mRNA_exon->set_feat_name($exon_feat_name);
	$gene_obj->add_mRNA_exon_obj($mRNA_exon);
	my $cds_exon = new CDS_exon_obj($exon_end5, $exon_end3);
	$mRNA_exon->set_CDS_exon_obj($cds_exon);
    }
    
    ## set the asmbl_id 
    my $query = "select asmbl_id from asm_feature where feat_name = \"$feat_name\"\n";
    my $asmbl_id = &first_result_sql ($dbproc, $query);
    $gene_obj->{asmbl_id} = $asmbl_id;

        
    ## refine and return gene_obj
    $gene_obj->refine_gene_object();
    return ($gene_obj);
}

1; #endof package

