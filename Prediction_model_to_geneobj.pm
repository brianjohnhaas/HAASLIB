#!/usr/local/bin/perl

## There's no reason to create objects from this class.  Use methods as 
## static methods.  Simply returns Gene_obj's based on model feat_names.


package Prediction_model_to_geneobj;

use strict;
use DBI;
use Egc_library;
use Gene_obj;

# get gene based on model feat_name
sub get_prediction {
    my ($dbproc, $ev_type, $feat_name) = @_;
    unless ($feat_name =~ /^\d+\.m\d+/) {
	die "Can't obtain gene obj for feat_name: ($feat_name)\n";
    }
    

    my $gene_obj = new Gene_obj();
    my $tu_feat_name = $feat_name;
    $tu_feat_name =~ s/\.m/\.tPRED/; ## Make different from a regular TU feat_name so it can be easily distinguished.
    
    ## set up gene parameters:
    $gene_obj->{TU_feat_name} = $tu_feat_name;
    $gene_obj->{Model_feat_name} = $feat_name;
    $gene_obj->{is_pseudogene} = 0;
    $gene_obj->{com_name} = "$ev_type model $feat_name"; 
        
    ## build the gene structure
    my $has_gene_structure = 0;
    my $query = "select a_exon.end5, a_exon.end3, a_exon.feat_name \n"
	. "from asm_feature a_exon, feat_link f_m_e, phys_ev p_exon \n"
	. "where f_m_e.parent_feat = \"$feat_name\" and f_m_e.child_feat = a_exon.feat_name\n"
	. "and a_exon.feat_name = p_exon.feat_name and p_exon.ev_type = \"$ev_type\"\n and a_exon.feat_type = \"exon\"\n";

    
    my @results = &do_sql ($dbproc, $query);
    foreach my $result (@results) {
	my ($exon_end5, $exon_end3, $exon_feat_name) = split (/\t/, $result);
	my $mRNA_exon = new mRNA_exon_obj($exon_end5, $exon_end3);
	$mRNA_exon->set_feat_name($exon_feat_name);
	$gene_obj->add_mRNA_exon_obj($mRNA_exon);
	my $cds_exon = new CDS_exon_obj($exon_end5, $exon_end3);
	$mRNA_exon->set_CDS_exon_obj($cds_exon);
	$has_gene_structure = 1;
    }
    
    if ($has_gene_structure) {

	## set the asmbl_id 
	my $query = "select asmbl_id from asm_feature where feat_name = \"$feat_name\"\n";
	my $asmbl_id = &first_result_sql ($dbproc, $query);
	$gene_obj->{asmbl_id} = $asmbl_id;
	
	## refine and return gene_obj
	$gene_obj->refine_gene_object();
	return ($gene_obj);
    } else {
	return (); # No gene structure
    }
}


## Get all gene objects on asmbl_id
sub get_predicted_models_on_assembly {
    my ($dbproc, $ev_type, $asmbl_id) = @_;
    
    my @predicted_gene_objs;

    my $query = "select a.feat_name from asm_feature a, phys_ev p where a.asmbl_id = $asmbl_id and a.feat_name = p.feat_name and p.ev_type = \"$ev_type\" and a.feat_type = 'model' ";
    my @prediction_feat_names = &do_sql($dbproc, $query);
    foreach my $model_feat_name (@prediction_feat_names) {
	my $gene_obj = get_prediction($dbproc, $ev_type, $model_feat_name);
	push (@predicted_gene_objs, $gene_obj);
    }
        
    return (@predicted_gene_objs);
}

1; #EOM
