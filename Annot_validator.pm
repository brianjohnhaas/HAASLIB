#!/usr/local/bin/perl

our ($SEE, $DEBUG);

package Annot_validator;

## There's no reason to create an instance of this class.  Class method will be used to validate
## a gene model.

## Requirements for annot validator:
#   -gene components must have same orientation as the TU.
#   -gene component coordinates must behave as expected; ie. CDS must reside within the exon.
#   -protein translations must exist for non-pseudogenes and must be complete ORFs.
#   -consensus splice sites must exist at internal exon boundaries.

require Exporter;
use lib ($ENV{EGC_SCRIPTS}, $ENV{EUK_MODULES});
use Egc_library;
use strict;
use DBI;
use Data::Dumper;


## Current validations:
# -A product name is available in the ident_xref table
# -The TU has a gene structure (models, exons) 
# -consensus GT/GC-AG donor-acceptor splice sites.
# -models, exons, and CDS's have the same strand orientation as the TU
# -models, exons, and CDS's are encapsulated within the coordinates of the TU
# -a working model has only one TU
# -an exon can have at most one CDS
# -exons cannot overlap each other
# -introns are at least 10 bp long.
# -protein sequence starts with M, ends with *, and lacks intervening stops.
# -CDS sequence starts with ATG and ends with a stop codon {TAG,TGA,TAA}.
# -feature coordinates are within the length of the contig.
# - a non-pseudogene encodes at least one CDS exon.
##



use vars qw (@ISA @EXPORT); ## set in script using this module for verbose output.
@ISA = qw(Exporter);
@EXPORT = qw(validate_gene);

#global 
my $asmbl_seq_ref;
my $seq_length;

## validate gene 
sub validate_gene {
    my ($dbproc, $asmbl_id, $TU, $asmbl_seq_reference) = @_;
    $asmbl_seq_ref = $asmbl_seq_reference;
    $seq_length = length($$asmbl_seq_ref);
    my ($error_text); #hold errors for gene analysis
    my $err_ref = \$error_text; #pass this ref to subs.
    &validate_ident_info ($dbproc, $TU, $err_ref); #should be an entry in the ident table
    &validate_gene_structure ($dbproc, $TU, $asmbl_id, $err_ref); ## do all the serious work here.
    if ($error_text) {
	return ("ASMBL_ID: $asmbl_id\tTU: $TU\n" . $error_text . "\n\n");
    } else {
	return (undef());
    }
}



################################################################################################################3


## SUPPORTING METHODS ##

####
sub validate_ident_info {
    my ($dbproc, $TU, $err_ref) = @_;
    #this can be expounded upon later.
    my $query = "select feat_name, ident_val from ident_xref where feat_name = \"$TU\" and xref_type = 'product name' and relrank = 1 \n\n";
    my $result = &first_result_sql ($dbproc, $query);
    unless ($result) {
	&add_error($err_ref, "IDENT_XREF: apparent, no info exists here for $TU\n");
	return;
    }
    my ($feat_name, $com_name) = split (/\t/, $result);
    unless ($com_name =~ /\w/) {
	&add_error($err_ref, "IDENT_XREF: no com_name assignment for gene $TU\n");
	return;
    }
}

####
sub validate_gene_structure {
    my ($dbproc, $TU, $asmbl_id, $err_ref) = @_;
    ## here, we analyse the gene structure, splice site data, and ORF composition.
    my $query = "select a_TU.feat_name, a_TU.end5, a_TU.end3,\n"
	. "a_Model.feat_name, a_Model.end5, a_Model.end3,\n"
	. "a_Exon.feat_name, a_Exon.end5, a_Exon.end3,\n"
	. "a_CDS.feat_name, a_CDS.end5, a_CDS.end3\n"
	. "from asm_feature a_TU, asm_feature a_Model, asm_feature a_Exon, asm_feature a_CDS,\n"
	. "feat_link f_t_m, feat_link f_m_e, feat_link f_e_c, \n"
	. "phys_ev p_Model, phys_ev p_Exon, phys_ev p_CDS\n"
	. "where a_TU.feat_name = \"$TU\" and a_TU.feat_type = \"TU\"\n"
	. "and a_TU.feat_name = f_t_m.parent_feat and f_t_m.child_feat = a_Model.feat_name\n"
	. "and a_Model.feat_name = p_Model.feat_name and p_Model.ev_type = \"working\"\n"
        . "and a_Model.feat_type = \"model\"\n"
	. "and a_Model.feat_name = f_m_e.parent_feat and f_m_e.child_feat = a_Exon.feat_name\n"
        . "and a_Exon.feat_name = p_Exon.feat_name and p_Exon.ev_type = \"working\"\n"
	. "and a_Exon.feat_type = \"exon\"\n"
	. "and a_Exon.feat_name *= f_e_c.parent_feat and f_e_c.child_feat *= a_CDS.feat_name\n"
        . "and a_CDS.feat_name *= p_CDS.feat_name and p_CDS.ev_type = \"working\"\n"
	. "and a_CDS.feat_type = \"CDS\"\n";

    print $query if $SEE;

    my @results = &do_sql ($dbproc, $query);
    unless (@results) { #no gene structure here
	
	&add_error ($err_ref, "No Gene Structure available for $TU\n");
	return;
    }
    
    
    ## Create data structure
    my (%TU_struct, %Model_struct, %Exon_struct, %CDS_struct);
    my %seen; #don't analyse the same data twice.
    foreach my $result (@results) {
	my ($TU_feat_name, $TU_end5, $TU_end3, 
	    $Model_feat_name, $Model_end5, $Model_end3, 
	    $Exon_feat_name, $Exon_end5, $Exon_end3, 
	    $CDS_feat_name, $CDS_end5, $CDS_end3) = split (/\t/, $result);
	$TU_struct{$TU_feat_name}->{end5} = $TU_end5;
	$TU_struct{$TU_feat_name}->{end3} = $TU_end3;
	if ($Model_feat_name) {
	    if (!$seen{$Model_feat_name}) {
		push (@{$TU_struct{$TU_feat_name}->{Models}}, $Model_feat_name); #implicit array ref
		$seen{$Model_feat_name} = 1;
	    }
	    $Model_struct{$Model_feat_name}->{end5} = $Model_end5;
	    $Model_struct{$Model_feat_name}->{end3} = $Model_end3;
	    
	    ## compound exons with Model indices
	    if ($Exon_feat_name) {
		if (!$seen{$Exon_feat_name}) {
		    push (@{$Model_struct{$Model_feat_name}->{Exons}}, $Exon_feat_name); #implicit array ref
		    $seen{$Exon_feat_name} = 1;
		}
		$Exon_struct{$Exon_feat_name}->{end5} = $Exon_end5;
		$Exon_struct{$Exon_feat_name}->{end3} = $Exon_end3;
		if ($CDS_feat_name) {
		    push (@{$Exon_struct{$Exon_feat_name}->{CDS}}, $CDS_feat_name);
		    $CDS_struct{$CDS_feat_name}->{end5} = $CDS_end5;
		    $CDS_struct{$CDS_feat_name}->{end3} = $CDS_end3;
		}
	    }
	}
    }

    ## Traverse data structure, identify problems.
    foreach my $TU_feat_name (keys %TU_struct) {
	my $TU_end5 = $TU_struct{$TU_feat_name}->{end5};
	my $TU_end3 = $TU_struct{$TU_feat_name}->{end3};
	&check_feat_within_contig($TU_feat_name, $err_ref, $TU_end5, $TU_end3);

	my $gene_orient = &get_orient ($TU_end5, $TU_end3);
	my $is_pseudogene = &get_pseudogene_status ($dbproc, $TU_feat_name);
	print "TU: $TU_feat_name\t$gene_orient\t$TU_end5\t$TU_end3\n" if $SEE;
	my $number_of_models = ($#{$TU_struct{$TU_feat_name}->{Models}} + 1);
	my @model_coord_structs;
	my @model_feat_names; #tracking
	foreach my $model_feat_name (@{$TU_struct{$TU_feat_name}->{Models}}) {
	    push (@model_feat_names, $model_feat_name);
	    &confirm_one_TU($dbproc, $model_feat_name, $err_ref);
	    my $model_end5 = $Model_struct{$model_feat_name}->{end5};
	    my $model_end3 = $Model_struct{$model_feat_name}->{end3};
	    &check_feat_within_contig($model_feat_name, $err_ref, $model_end5, $model_end3);

	    push (@model_coord_structs, { feat_name=>$model_feat_name,
					  end5 => $model_end5,
					  end3 => $model_end3 });
	    
	    my $model_orient = &confirm_orientation ($gene_orient, $model_end5, $model_end3, "model", $err_ref);
	    &confirm_encapsulation ($TU_end5, $TU_end3, $model_end5, $model_end3, "TU/model", $err_ref);
	    print "Model: $model_feat_name\t$model_orient\t$model_end5\t$model_end3\n" if $SEE;
	    my @exons = sort {$Exon_struct{$a}->{end5}<=>$Exon_struct{$b}->{end5}} @{$Model_struct{$model_feat_name}->{Exons}};
	    @exons = reverse @exons if ($gene_orient eq '-');
	    my $num_exons = $#exons;
	    print "Model has " . ($num_exons+1) . " number of exons.\n" if $SEE;
	    my $seen_CDS = 0;
	    for (my $i = 0; $i <= $num_exons; $i++) {
		my $exon_feat_name = $exons[$i];
		my $exon_end5 = $Exon_struct{$exon_feat_name}->{end5};
		my $exon_end3 = $Exon_struct{$exon_feat_name}->{end3};
		&check_feat_within_contig ($exon_feat_name, $err_ref,$exon_end5, $exon_end3);
		my $exon_orient = &confirm_orientation($gene_orient, $exon_end5, $exon_end3, "exon", $err_ref);
		print "Exon: $exon_feat_name\t$exon_orient\t$exon_end5\t$exon_end3\n" if $SEE;
		&confirm_encapsulation($TU_end5, $TU_end3, $exon_end5, $exon_end3, "TU/exon", $err_ref);
		my $exon_type = &get_exon_type ($i, $num_exons);
		&confirm_splice_sites ($exon_end5, $exon_end3, $exon_orient, $exon_type, 
				       $exon_feat_name, $asmbl_seq_ref, $err_ref) if (!$is_pseudogene); 
		my $CDS_feat_name_array_ref = $Exon_struct{$exon_feat_name}->{CDS};
		&confirm_number_of_elements( ($#{$CDS_feat_name_array_ref} + 1), $err_ref, $exon_feat_name, "CDSs");
		if (my $CDS_feat_name =  $CDS_feat_name_array_ref->[0]) {
		    my $CDS_end5 = $CDS_struct{$CDS_feat_name}->{end5};
		    my $CDS_end3 = $CDS_struct{$CDS_feat_name}->{end3};
		    &check_feat_within_contig($CDS_feat_name, $CDS_end5, $CDS_end3);
		    my $CDS_orient = &confirm_orientation($gene_orient, $CDS_end5, $CDS_end3, "CDS", $err_ref);
		    &confirm_encapsulation($exon_end5, $exon_end3, $CDS_end5, $CDS_end3, "exon:$exon_feat_name/CDS:$CDS_feat_name", $err_ref);
		    print "CDS: $CDS_feat_name\t$CDS_orient\t$CDS_end5\t$CDS_end3\n\n" if $SEE;
		    $seen_CDS = 1;
		}
	    }
	    if ((! $seen_CDS) && !$is_pseudogene) {
		&add_error ($err_ref, "No CDS exon reported for gene.\n");
	    }
	    
	    &check_exons_for_overlap(\%Exon_struct, \@exons, $err_ref);
	    &check_intron_lengths(\%Exon_struct, \@exons, $err_ref);
	    my %partials = &get_partial_info_from_db($dbproc, $model_feat_name);
	    &examine_protein_sequence ($dbproc, $model_feat_name, $err_ref, \%partials) if (!$is_pseudogene);
	    &examine_CDS_sequence ($dbproc, $model_feat_name, $err_ref, $is_pseudogene, \%partials);
	}
	
	## Check multiple isoforms for consistent coordinates, and report alt splicing.
	&check_overlapping_models(\@model_coord_structs, $err_ref);
	
	&check_working_models_no_exons($dbproc, $TU_feat_name, \@model_feat_names, $err_ref);
    }
}



####
sub confirm_splice_sites {
    my ($exon_end5, $exon_end3, $exon_orient, $exon_type, $exon_feat_name, $asmbl_seq_ref, $err_ref) = @_;
    print "exon_type: $exon_type\n" if $SEE;
    if ($exon_type eq "gene") { return;} ## no splicing here.
    my ($coord1, $coord2) = sort {$a<=>$b} ($exon_end5, $exon_end3);
    ## get two coordinate sets corresponding to potential splice sites
    my $splice_1_start = $coord1-2-1;
    my $splice_2_start = $coord2-1+1;
    print "confirming splice sites at "  . ($splice_1_start +1) . " and " . ($splice_2_start + 1) . "\n"if $SEE;
    my $splice_1 = substr ($$asmbl_seq_ref, $splice_1_start, 2);
    my $splice_2 = substr ($$asmbl_seq_ref, $splice_2_start, 2);
    if ($SEE) {
	my $exon_length = $coord2 - $coord1 + 1;
	my $exonseq = lc (substr ($$asmbl_seq_ref, $splice_1_start, $exon_length + 4));
	if ($exon_orient eq "-") { $exonseq = &reverse_complement($exonseq);}
	$exonseq =~ /^(\w{2})/;
	my $ac = uc $1;
	$exonseq =~ s/^\w{2}/$ac/;
	$exonseq =~ /(\w{2})$/;
	my $d = uc $1;
	$exonseq =~ s/\w{2}$/$d/;
	print "ExonSeq: $exonseq\n";
    }
 
    my ($acceptor, $donor) = ($exon_orient eq '+') ? ($splice_1, $splice_2) : (&reverse_complement($splice_2), &reverse_complement($splice_1)); 
    my $check_acceptor = ($acceptor =~ /ag/i);
    my $check_donor = ($donor =~ /gt|gc/i);
    ## associate results of checks with exon type.
    if ($exon_type eq "initial" || $exon_type eq "internal") {
	unless ($check_donor) {
	    &add_error ($err_ref, "non-consensus $donor donor splice site at $exon_feat_name\n");
	}
    }
    if ($exon_type eq "internal" || $exon_type eq "terminal") {
	unless ($check_acceptor) {
	    &add_error ($err_ref, "non-consensus $acceptor acceptor splice site at $exon_feat_name\n");
	}
    }
    return;
}



####
sub get_exon_type {
    my ($curr_exon_num, $num_exons) = @_;

    ## count starts at 0; array based.

    ## types defined
    # gene: single exon gene
    # initial :first exon of a multi-exon gene
    # internal: internal exon of a multi-exon gene
    # terminal: last exon of a multi-exon gene


    if ($curr_exon_num == 0) {
	if ($num_exons == 0) {
	    return "gene";
	} else {
	    return "initial";
	}
    } elsif ($curr_exon_num == $num_exons) {
	return "terminal";
    } else {
	return "internal";
    }
}



####
sub confirm_orientation {
    my ($gene_orient, $coord1, $coord2, $type, $err_ref) = @_;
    my $orient = &get_orient($coord1, $coord2);
    if (($orient) && ($gene_orient ne $orient)) {
	&add_error($err_ref, "$type has opposite orientation to gene ($gene_orient)\n");
    }
    return ($orient);
}


####
sub confirm_encapsulation {
    my ($coord_big1, $coord_big2, $coord_sm1, $coord_sm2, $type, $err_ref) = @_;
    ($coord_big1, $coord_big2) = sort {$a<=>$b} ($coord_big1, $coord_big2);
    ($coord_sm1, $coord_sm2) = sort {$a<=>$b} ($coord_sm1, $coord_sm2);
    unless ( ($coord_sm1 >= $coord_big1) && ($coord_sm2 <= $coord_big2) ) {
	print "Not Encapsulated: ($coord_sm1, $coord_sm2, $type) within ($coord_big1, $coord_big2)\n" if $SEE;
	&add_error($err_ref, "$type, coords not encapsulated\n");
    }
}

####
sub get_orient {
    my ($coord1, $coord2) = @_;
    my $orient;
    if ($coord1 < $coord2) {
	$orient = '+';
    } elsif ($coord1 > $coord2) {
	$orient = '-';
    }
    return ($orient);
}

####
sub get_pseudogene_status {
    my ($dbproc, $TU_feat_name, $err_ref) = @_;
    my $query = "select is_pseudogene from ident where feat_name = \"$TU_feat_name\"\n";
    my $status = &first_result_sql ($dbproc, $query);
    if ($status == 1) {
	return (1);
    } else {
	return(0);
    }
}

####
sub examine_protein_sequence {
    my ($dbproc, $model_feat_name, $err_ref, $partials_ref) = @_;
    my $query = "select protein from asm_feature where feat_name = \"$model_feat_name\"\n";
    my $prot_seq = &first_result_sql ($dbproc, $query);
    unless ($prot_seq) {
	&add_error ($err_ref, "No protein sequence.\n");
	return;
    }

    ## see if it is a gene_synonym
    my @gene_synonyms = &gather_all_gene_synonyms($dbproc, $model_feat_name);
    if (@gene_synonyms) {
	$$err_ref .= "\t" . '***' . " gene_synonym here, I may be a partial gene" . ' ***' . "\n";
	
	if ($model_feat_name =~ /^(\d+)\./) {
	    my $asmbl_id = $1;
	    foreach my $syn (@gene_synonyms) {
		if ($syn =~ /^$asmbl_id\./) {
		    &add_error($err_ref, "Model $model_feat_name has synonym $syn on the same assembly.\n");
		}
	    }
	}
    }
    
    if ($#gene_synonyms >= 1) {
	&add_error ($err_ref, "Model $model_feat_name has multiple gene synonyms.  Could be problematic.  GSs: @gene_synonyms\n");
    }

   
    my $first_char = substr ($prot_seq, 0, 1);
    unless ($first_char =~ /M/i) {
	unless ($partials_ref->{'5prime'}) {
	    &add_error ($err_ref,  "Protein seq doesn't start with Methionine.....[$first_char] instead.\n");
	}
    }
    my $prot_length = length ($prot_seq);
    my $last_char = substr ($prot_seq, ($prot_length - 1), 1);
    unless ($last_char =~ /\*/) {
	unless ($partials_ref->{'3prime'}) {
	    &add_error ($err_ref,  "Protein seq doesn't terminate with a stop codon....[$last_char] instead.\n");
	}
    }

    ## remove stop codon, continue analysis:
    $prot_seq =~ s/\*$//;
    $prot_length = length($prot_seq);
    my $num_stops = &get_number_stops ($prot_seq);
    if ($num_stops > 0) {
	&add_error ($err_ref,  "corrupt protein sequence: [$num_stops] stops in protein sequence\n");
	return;
    }

    ## check protein length
    if ($prot_length < 50) {
	&add_error($err_ref, "Warning: Protein is very short (length = $prot_length)\n");
    }
    
}



####
sub get_number_stops {
    my ($prot_seq) = @_;
    my $stop_num = 0;
    while ($prot_seq =~ /\*/g) {
	$stop_num++;
    } 
    return ($stop_num);
}

####
sub confirm_one_TU {
    my ($dbproc, $model_feat_name, $err_ref) = @_;
    my $query = "select parent_feat from feat_link where child_feat = \"$model_feat_name\" and parent_feat like \"%.t%\"\n";
    my @results = &do_sql ($dbproc, $query);
    if ($#results >= 1) {
	&add_error($err_ref, "Model $model_feat_name has multiple TUs: @results\n");
    }
}

####
sub add_error {
    my ($err_ref, $comment) = @_;
    $comment = "\tERROR:\t$comment";
    $$err_ref .= $comment;
}

####
sub add_warning {
    my ($err_ref, $comment) = @_;
    $comment = "\tWARNING:\t$comment";
    $$err_ref .= $comment;
}


####
sub confirm_number_of_elements {
    ## each exon should have only a single CDS!
    my ($num_elements, $err_ref, $exon_feat_name, $type) = @_;
    if ($num_elements > 1) {
	my $comment = "$exon_feat_name has multiple $type\n";
	&add_error($err_ref, $comment);
    }
}


####
sub examine_CDS_sequence {
    my ($dbproc, $model_feat_name, $err_ref, $is_pseudogene, $partials_ref) = @_;
    
    unless ($is_pseudogene) {
        
        my $query = "select sequence from asm_feature where feat_name = \"$model_feat_name\"\n";
        my $sequence = &first_result_sql ($dbproc, $query);
        if (!$sequence) {
            &add_error ($err_ref, "No CDS sequence could be found for this gene. Possibly no CDS-exons?\n");
        } elsif (!$is_pseudogene) {
            chomp $sequence; #shouldn't need to do this, but just to be safe.
            my $cds_length = length ($sequence);
            my $start_codon = substr ($sequence, 0, 3);
            my $stop_codon = substr ($sequence, $cds_length - 3, 3);
            if ($start_codon !~ /atg/i) {
                unless ($partials_ref->{'5prime'}) {
                    &add_error ($err_ref, "CDS sequence does not begin with ATG; instead [$start_codon]\n");
                }
            }
            if ($stop_codon !~ /tag|tga|taa/i) {
                unless ($partials_ref->{'3prime'}) {
                    &add_error ($err_ref, "CDS sequence does not include a STOP codon; instead [$stop_codon]\n");
                }
            }
        }
    }
}

sub check_exons_for_overlap {
    my ($exon_struct_ref, $exon_list, $err_ref) = @_;
    my $i = 0;
   
    my @exons = @$exon_list;
    my $num_exons = $#exons + 1;

    for ($i = 0; $i < $num_exons - 1; $i++) {
	for (my $j = $i + 1; $j < $num_exons; $j++) {
	    my $feat_name1 = $exons[$i];
	    my $feat_name2 = $exons[$j];
	    
	    my ($a1, $a2) = sort {$a<=>$b} ($exon_struct_ref->{$feat_name1}->{end5}, $exon_struct_ref->{$feat_name1}->{end3});
	    my ($b1, $b2) = sort {$a<=>$b} ($exon_struct_ref->{$feat_name2}->{end5}, $exon_struct_ref->{$feat_name2}->{end3});
	   
	    if ($a1 < $b2 && $a2 > $b1) { #overlap
		&add_error($err_ref, "Exons $feat_name1 and $feat_name2 overlap.\n");
	    }
	}
    }
        
}


sub check_intron_lengths {
    my ($exon_struct_ref, $exon_list, $err_ref) = @_;
    my $i = 0;
    
    my @exons = @$exon_list; #already sorted by end5
    my $num_exons = $#exons + 1;
    
    if ($num_exons <= 1) {
	return();
    }

    my $min_intron_length = 10; #arbitrary (targeting 1 or 2 bp introns).

    for ($i = 1; $i < $num_exons; $i++) {
	my $feat_name1 = $exons[$i-1];
	my $feat_name2 = $exons[$i];
	    
	my ($a1, $a2) = ($exon_struct_ref->{$feat_name1}->{end5}, $exon_struct_ref->{$feat_name1}->{end3});
	my ($b1, $b2) = ($exon_struct_ref->{$feat_name2}->{end5}, $exon_struct_ref->{$feat_name2}->{end3});
	   
	my $intron_length = abs ($b1 - $a2) + 1;
	print "intron length: ($feat_name1, $feat_name2) = $intron_length\n" if $SEE;
	if ($intron_length < $min_intron_length) {
	    &add_error ($err_ref, "short intron between $feat_name1 and $feat_name2 of length $intron_length bp.\n");
	}
    }
}

####
sub check_feat_within_contig {
    my ($feat_name, $err_ref, $end5, $end3) = @_;
    if ($end5 > $seq_length) {
	&add_error($err_ref, "end5 of $feat_name is > seq_length ($end5 < $seq_length).\n");
    }
    if ($end3 > $seq_length) {
	&add_error($err_ref, "end3 of $feat_name is > seq_length ($end3 < $seq_length).\n");
    }
}



####
sub check_overlapping_models {
    my ($model_coords_structs_aref, $err_ref) = @_;
    
    if (scalar (@$model_coords_structs_aref) > 1) {
	
	my %models;
	
	for (my $i=0; $i <$#$model_coords_structs_aref; $i++) {
	    my $struct_i = $model_coords_structs_aref->[$i];
	    
	    my ($feat_name_i, $end5_i, $end3_i) = ($struct_i->{feat_name},
						   $struct_i->{end5},
						   $struct_i->{end3} );
	    my ($orient_i, $lend_i, $rend_i) = ($end5_i < $end3_i) ? ('+', $end5_i, $end3_i) : ('-', $end3_i, $end5_i);

	    $models{$feat_name_i} = 1;

	    for (my $j=$i+1; $j <= $#$model_coords_structs_aref; $j++) {
		my $struct_j = $model_coords_structs_aref->[$j];
		
		my ($feat_name_j, $end5_j, $end3_j) = ($struct_j->{feat_name},
						       $struct_j->{end5},
						       $struct_j->{end3} );
		my ($orient_j, $lend_j, $rend_j) = ($end5_j < $end3_j) ? ('+', $end5_j, $end3_j) : ('-', $end3_j, $end5_j);
		
		if (! ($lend_i < $rend_j && $rend_i > $lend_j)) {
		    &add_error($err_ref, "Isoform models $feat_name_i and $feat_name_j do not overlap.\n");
		}
		
		if ($orient_i ne $orient_j) {
		    &add_error($err_ref, "Isoform models $feat_name_i ($orient_i) and $feat_name_j ($orient_j) have opposite orientations\n");
		}

		$models{$feat_name_j} = 1;
	    }
	}
	
	my @modelList = keys (%models);
	my $model_list = join (", ", @modelList);
	
	&add_warning($err_ref, "multiple models (isoforms) exist: $model_list\n");

    }
}


####
sub check_working_models_no_exons {
    my ($dbproc, $TU_feat_name, $model_list_aref, $err_ref) = @_;
    my %existingModels;
    foreach my $model (@$model_list_aref) {
	$existingModels{$model} = 1;
    }

    ## get list of all working models linked to TU:
    my $query = "select a.feat_name from asm_feature a, phys_ev p, feat_link f where f.parent_feat = \"$TU_feat_name\" and f.child_feat = a.feat_name and a.feat_type = 'model' and a.feat_name = p.feat_name and p.ev_type = 'working'";
    my @models = &do_sql($dbproc, $query);
    foreach my $model (@models) {
	unless ($existingModels{$model}) {
	    &add_error($err_ref, "working model $model lacks gene structure.\n");
	}
    }
}


####
sub get_partial_info_from_db {
    my ($dbproc, $model_feat_name) = @_;
    
    my %partials;

    my $query = "select score, score2 from ORF_attribute where att_type = 'is_partial'";
    my $result = &first_result_sql($dbproc, $query);
    if ($result) {
	my ($prime5, $prime3) = split (/\t/, $result);
	$partials{'5prime'} = $prime5;
	$partials{'3prime'} = $prime3;
    }

    return (%partials);
}


1; #end of Annot_validator.pm



