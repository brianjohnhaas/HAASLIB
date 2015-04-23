#!/usr/local/bin/perl

## There's no reason to create objects from this class.  Use methods as 
## static methods.  Simply returns Gene_obj's based on model feat_names.


package DBmodel_to_geneobj;

use strict;
use DBI;
use lib ($ENV{EGC_SCRIPTS});
use Egc_library;
use Gene_obj;
use Gene_ontology;

# get gene based on model feat_name
sub get_gene {
    my ($dbproc, $feat_name) = @_;
    unless ($feat_name =~ /^\d+\.m\d+/) {
	die "Can't obtain gene obj for feat_name: ($feat_name)\n";
    }
    

    my $gene_obj = new Gene_obj();
    my $query = "select i.feat_name, i.locus, i.pub_locus, i.is_pseudogene, i.pub_comment, i.comment  from ident i, feat_link f where i.feat_name = f.parent_feat and f.child_feat = \"$feat_name\"\n";
    
    my @results = &do_sql_2D($dbproc, $query);
    my @x;
    my $top_result = shift @results;
    if (ref $top_result) {
	@x = @$top_result;
    }
    my ($tu_feat_name, $locus, $pub_locus, $is_pseudogene, $pub_comment, $comment) = @x;
    
    ## Get ident_xref info:
    
    # mappings of xref_type to Gene_obj attribute:
    my %xrefMap = ( "ec number" => {  1 => "ec_num",
				      2 => "secondary_ec_numbers" },
		    "gene name" => { 1 => "gene_name",
				     2 => "secondary_gene_names" },
		    
		    "gene symbol" => { 1 => "gene_sym",
				       2 => "secondary_gene_symbols" },

		    "product name" => { 1 => "com_name",
					2 => "secondary_product_names" }
		    );


    my $query = "select ident_val, xref_type, relrank from ident_xref where feat_name = \"$tu_feat_name\"";
    my @results = &do_sql_2D($dbproc, $query);
    
    foreach my $result (@results) {
	my ($ident_val, $xref_type, $relrank) = @$result;
	my $gene_att = $xrefMap{$xref_type};
	
	if (ref $gene_att) {
	    $gene_att = $gene_att->{$relrank};
	}
	unless ($gene_att) {
	    #print STDERR "Currently do not recognize xref_type of $xref_type, rank $relrank from ident_xref\n";
	    next;
	}
	
	if ($relrank == 1) {
	    $gene_obj->{$gene_att} = $ident_val;
	} else {
	    ## add to list of secondary attributes.
	    push (@{$gene_obj->{$gene_att}}, $ident_val);
	}
    }
    
    ## get com_name curated info:
    my $query = "select curated from asm_feature where feat_name = \"$tu_feat_name\"";
    my $com_name_curated = &first_result_sql($dbproc, $query);

    ## get model structure curation info:
    my $query = "select curated from asm_feature where feat_name = \"$feat_name\"";
    my $curated_gene_structure = &first_result_sql($dbproc, $query);
    
    
    ## set up gene parameters:
    $gene_obj->{TU_feat_name} = $tu_feat_name;
    $gene_obj->{Model_feat_name} = $feat_name;
    $gene_obj->{locus} = $locus;
    $gene_obj->{pub_locus} = $pub_locus;
    $gene_obj->{is_pseudogene} = $is_pseudogene;
    $gene_obj->{pub_comment} = $pub_comment;
    $gene_obj->{comment} = $comment;
    $gene_obj->{curated_com_name} = $com_name_curated;
    $gene_obj->{curated_gene_structure} = $curated_gene_structure;
    

    ## See if has model_pub_locus
    my $query = "select pub_locus, locus from ident where feat_name = \"$feat_name\"";
    my $result = &first_result_sql($dbproc, $query);
    my ($model_pub_locus, $model_locus) = split (/\t/, $result);
    if ($model_pub_locus) {
	$gene_obj->{model_pub_locus} = $model_pub_locus;
    }
    if ($model_locus) {
	$gene_obj->{model_locus} = $model_locus;
    }
    
    my $has_gene_structure = 0;
    
    ## build the gene structure
    $query = "select a_exon.end5, a_exon.end3, a_exon.feat_name \n"
	. "from asm_feature a_exon, feat_link f_m_e, phys_ev p_exon \n"
	. "where f_m_e.parent_feat = \"$feat_name\" and f_m_e.child_feat = a_exon.feat_name\n"
	. "and a_exon.feat_name = p_exon.feat_name and p_exon.ev_type = \"working\"\n and a_exon.feat_type = \"exon\"\n";
        
    my @results = &do_sql ($dbproc, $query);
    foreach my $result (@results) {
	my ($exon_end5, $exon_end3, $exon_feat_name) = split (/\t/, $result);
	my $mRNA_exon = new mRNA_exon_obj($exon_end5, $exon_end3);
	$mRNA_exon->set_feat_name($exon_feat_name);
	$gene_obj->add_mRNA_exon_obj($mRNA_exon);
	$has_gene_structure = 1;
	
	## check for CDS
	my $query = "select a_cds.end5, a_cds.end3, a_cds.feat_name from asm_feature a_cds, feat_link f_e_c, phys_ev p_cds \n"
	    . " where f_e_c.parent_feat = \"$exon_feat_name\" and f_e_c.child_feat = a_cds.feat_name and a_cds.feat_type = 'CDS' \n"
	    . " and a_cds.feat_name = p_cds.feat_name and p_cds.ev_type = 'working' ";
	my $result = &first_result_sql($dbproc, $query);
	if ($result) {
	    my ($cds_end5, $cds_end3, $cds_feat_name) = split (/\t/, $result);
	    
	    if ($cds_end5 && $cds_end3) {
		my $cds_exon = new CDS_exon_obj($cds_end5, $cds_end3);
		$cds_exon->set_feat_name($cds_feat_name);
		$mRNA_exon->set_CDS_exon_obj($cds_exon);
	    }
	}
    }
    
    if ($has_gene_structure) {

	## set the asmbl_id 
	my $query = "select asmbl_id from asm_feature where feat_name = \"$feat_name\"\n";
	my $asmbl_id = &first_result_sql ($dbproc, $query);
	$gene_obj->{asmbl_id} = $asmbl_id;
	
	## populate gene_synonym list
	my @gene_synonyms = gather_all_gene_synonyms($dbproc, $feat_name);
	$gene_obj->{gene_synonyms} = \@gene_synonyms;
	
	## Set the partial flags if partial:
	my $query = "select score, score2 from ORF_attribute where feat_name = \"$feat_name\" and att_type = \"is_partial\"";
	my $result = &first_result_sql($dbproc, $query);
	if ($result) {
	    my ($is_5prime_partial, $is_3prime_partial) = split (/\t/, $result);
	    if ($is_5prime_partial) {
		$gene_obj->set_5prime_partial(1);
	    }
	    if ($is_3prime_partial) {
		$gene_obj->set_3prime_partial(1);
	    }
	}
	
	
	## Get the gene ontology assignments:
	my @gene_ontology_objs = &retrieve_GO_assignments($dbproc, $tu_feat_name);
	if (@gene_ontology_objs) {
	    $gene_obj->add_gene_ontology_objs (@gene_ontology_objs);
	}
	
	
	## refine and return gene_obj
	$gene_obj->refine_gene_object();
	return ($gene_obj);
    } else {
	return (); # No gene structure for requested model.
    }
}


## return list of gene models based on TU
sub get_gene_models_via_TU {
    my ($dbproc, $tu_feat_name, $alt_splice_rep_flag) = @_;
    ## if alt_splice_rep_flag, then only a single gene_obj will be returned with additional isoforms set within the additional_isoform list of the single gene_obj


    my $query = "select a.feat_name from asm_feature a, feat_link f, phys_ev p where f.parent_feat = \"$tu_feat_name\" and f.child_feat = a.feat_name and a.feat_type = \"model\" and a.feat_name = p.feat_name and p.ev_type = \"working\"\n";
    my @model_feats = &do_sql ($dbproc, $query);
    my @genemodels;
    foreach my $model_feat (@model_feats) {
	push (@genemodels, &get_gene ($dbproc, $model_feat));
    }
    
    if ($alt_splice_rep_flag && $#genemodels > 0) {
	my @copy = @genemodels;
	my $template_obj = shift @copy;
	foreach my $gene_obj (@copy) {
	    $template_obj->add_isoform($gene_obj);
	}
	$template_obj->refine_gene_object();
	@genemodels = ($template_obj);
    }
    
    return (@genemodels);
}



## Get all gene objects on asmbl_id
sub get_all_genes_on_assembly {
    my ($dbproc, $asmbl_id, $TU_based_only_flag, $First_model_only_flag, $alt_splice_rep_flag) = @_;
    ## set the TU_based_only_flag if you want to exclude all RNA genes.
    ## First_model_only_flag: if set, only the first isoform of multiple models will be included.


    ## get gene list
    my $query = "select distinct a.feat_name from ident i, asm_feature a where i.feat_name = a.feat_name and a.asmbl_id = $asmbl_id and a.feat_type = \"TU\" order by a.end5 \n";
    my @TUs = &do_sql ($dbproc, $query);
    my @complete_gene_list;
    foreach my $TU (@TUs) {
	my @genes = &get_gene_models_via_TU($dbproc, $TU, $alt_splice_rep_flag);
	if (@genes) {
	    if ($First_model_only_flag) {
		push (@complete_gene_list, $genes[0]);
	    } else {
		push (@complete_gene_list, @genes);
	    }
	}
    }
    
    unless ($TU_based_only_flag) {
	my @RNA_genes = &get_RNA_genes($dbproc, $asmbl_id);
	if (@RNA_genes) {
	    push (@complete_gene_list, @RNA_genes);
	}
    }
    return (@complete_gene_list);
}


## get tRNA and rRNA genes 

## For now, keeping it simple... RNA genes have only one coordinate set.
sub get_RNA_genes {
    my ($dbproc, $asmbl_id) = @_;
    my @rna_genes;

    foreach my $feat_type ("tRNA", "rRNA", "snoRNA", "snRNA", "slRNA", "ncRNA", "sRNA", "miRNA") {
	my $query = "select a.feat_name, a.feat_type, a.end5, a.end3, coalesce(i.com_name, o.score), coalesce (i.pub_locus, o.score5), i.locus from asm_feature a, ORF_attribute o, ident i where a.asmbl_id = $asmbl_id and a.feat_name *= o.feat_name and a.feat_type = \"$feat_type\" and a.feat_name *= i.feat_name and o.att_type = \"$feat_type\"";
	my @results = &do_sql ($dbproc, $query, $;);
	foreach my $result (@results) {
	    my ($feat_name, $feat_type, $end5, $end3, $name, $pub_locus, $locus) = split (/$;/, $result);
	    my $gene_obj = new Gene_obj();
	    my $mRNA_exon =  new mRNA_exon_obj($end5, $end3);
	    $gene_obj->add_mRNA_exon_obj($mRNA_exon);
	    $gene_obj->{TU_feat_name} = $gene_obj->{Model_feat_name} = $feat_name;
	    $gene_obj->{pub_locus} = $pub_locus;
	    $gene_obj->{locus} = $locus;
	    $gene_obj->{com_name} = $name;
	    $gene_obj->{gene_type} = lc $feat_type;
	    $gene_obj->refine_gene_object();
	    push (@rna_genes, $gene_obj);
	}
    }
    return (@rna_genes);
}



## Gene to XML format.
sub TU_to_XML {
    my ($dbproc, $tu_feat_name) = @_;
    my $xml = "<?xml version=\"1.0\"?>\n";
    $xml .= <<_EODOCTYPE;
    <!DOCTYPE GENEXML [
		       <!ELEMENT GENEXML (TU?)>
		       <!ELEMENT TU (FEAT_NAME, PUB_LOCUS, LOCUS, ALT_LOCUS, COM_NAME, PUB_COMMENT, IS_PSEUDOGENE, MODEL+)>
		       <!ELEMENT MODEL (FEAT_NAME, EXON+)>
		       <!ELEMENT EXON (FEAT_NAME, COORDS, CDS?)>
		       <!ELEMENT CDS (FEAT_NAME, COORDS)>
		       <!ELEMENT FEAT_NAME (\#PCDATA)>
		       <!ELEMENT PUB_LOCUS (\#PCDATA)>
		       <!ELEMENT LOCUS (\#PCDATA)>
		       <!ELEMENT ALT_LOCUS (\#PCDATA)>
		       <!ELEMENT COM_NAME (\#PCDATA)>
		       <!ELEMENT PUB_COMMENT (\#PCDATA)>
		       <!ELEMENT IS_PSEUDOGENE (\#PCDATA)>
		       <!ELEMENT COORDS (\#PCDATA)> ]>
_EODOCTYPE
			   
    $xml .= "\n<GENEXML>\n";
    my @geneobjs = get_gene_models_via_TU($dbproc, $tu_feat_name);
    if (@geneobjs) {
	## create xml
	my $template = $geneobjs[0];
	$xml .= "<TU>";
	$xml .= "<FEAT_NAME>$template->{TU_feat_name}</FEAT_NAME>";
	$xml .= "<PUB_LOCUS>$template->{pub_locus}</PUB_LOCUS>";
	$xml .= "<LOCUS>$template->{locus}</LOCUS>";
	$xml .= "<ALT_LOCUS>$template->{alt_locus}</ALT_LOCUS>";
	$xml .= "<COM_NAME>$template->{com_name}</COM_NAME>";
	$xml .= "<PUB_COMMENT>$template->{pub_comment}</PUB_COMMENT>";
	$xml .= "<IS_PSEUDOGENE>$template->{is_pseudogene}</IS_PSEUDOGENE>";
	
	## Now, illustrate gene structures: only need exon and cds info for each model
	foreach my $model (@geneobjs) {
	    $xml .= "<MODEL>";
	    $xml .= "<FEAT_NAME>$model->{Model_feat_name}</FEAT_NAME>";
	    my @exons = $model->get_exons();
	    foreach my $exon (@exons) {
		my @coords = $exon->get_coords();
		$xml .= "<EXON>\n";
		$xml .= "<FEAT_NAME>$exon->{feat_name}</FEAT_NAME>";
		$xml .= "<COORDS>$coords[0]-$coords[1]</COORDS>";
		## cds?
		my $cdsobj = $exon->get_CDS_obj();
		if (ref $cdsobj) {
		    my @coords = $cdsobj->get_coords();
		    $xml .= "<CDS>";
		    $xml .= "<FEAT_NAME>$cdsobj->{feat_name}</FEAT_NAME>";
		    $xml .= "<COORDS>$coords[0]-$coords[1]</COORDS>";
		    $xml .= "</CDS>";
		}
		$xml .= "</EXON>";
	    }
	    $xml .= "</MODEL>";
	}
	$xml .= "</TU>";
    }
    $xml .= "</GENEXML>";
    $xml =~ s/\"/\'/g;
    return ($xml);
}


####
sub retrieve_GO_assignments {
    my ($dbproc, $tu_feat_name) = @_;
    my $query = "select grl.id, grl.go_id, grl.assigned_by, grl.qualifier, gt.name, gt.type from go_role_link grl, common..go_term gt where grl.feat_name = \"$tu_feat_name\" and grl.go_id = gt.go_id";
    my @results = &do_sql_2D($dbproc, $query);
    if (@results) {
	## Got go assignments:
	my @gene_ontology_objs;
	foreach my $result (@results) {
	    my ($role_link_id, $go_id, $assignedby, $qualifier, $description, $type) = @$result;
	    my $gene_ontology_obj = new Gene_ontology($go_id, $description, $type);
	    $gene_ontology_obj->{qualifier} = $qualifier if ($qualifier =~ /\w/);
	    $gene_ontology_obj->{assignby} = $assignedby;
	    
	    ## get all evidences:
	    my $query = "select distinct ge.ev_code, ge.evidence, ge.with_ev from go_evidence ge where ge.role_link_id = $role_link_id";
	    my @results = &do_sql_2D($dbproc, $query);
	    foreach my $result (@results) {
		my ($ev_code, $evidence, $with_ev) = @$result;
		$gene_ontology_obj->add_evidence($ev_code, $evidence, $with_ev);
	    }
	    push (@gene_ontology_objs, $gene_ontology_obj);
	}
	return (@gene_ontology_objs);
    } else {
	return (); #no assignments
    }
}










1; #endof package DBmodel_to_geneobj
