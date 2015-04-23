#!/usr/local/bin/perl

package Genbankxmlparser;

use strict;
use Gene_obj;
use Xml_parser_facilitator;
use Data::Dumper;

## tree properties
#         tree node: (node_name, node_ref)
#         array_ref: (attribute_hash, (next_node_name, next_node_ref) ......)
#             if a node_name == 0, the node_ref is pcdata

my $SEE; #single file-scoped lexical

sub new {
    shift;
    my $self = {
	tree => 0,     #XML::Parser (Style=>Tree)
	genes => 0,    #Gene_obj's created to store gene info
	sequence => 0, #nucleotide sequence of genbank accession
        locus => 0,    #locus field
	accession => 0, # genbank accession
	version => 0,   #submission version
	title  => 0, #title
	update_date => 0, #date the submission was last updated.
	affiliation => 0, #seq_group information.
	SEE => 0,       #Verbose flag 
	chromosome=>0
    };
    while (@_) {
	my $att = shift;
	my $value = shift;
	if (exists ($self->{$att})) {
	    $self->{$att} = $value;
	}
    }
    unless (ref ($self->{tree}) eq "ARRAY") {
	die "I require an XML::Parser Style=Tree object\n";
    }
    $SEE = $self->{SEE};
    my $tree = $self->{tree};

    ## Specify Nodes for Desired Component Branches
    my $gene_components = "Seq-feat";
    my $sequence_data = "Seq-data";
    my $locus_info = "Textseq-id";
    my $title_info = "Seqdesc";
    my $update_date_info = "Seqdesc_update-date";
    my $affiliation = "Affil";
    my $chromosome = "SubSource_name";

    ## %gene_data holds the branch references for the Node classes specified above
    my %gene_data = ( ## all handlers are referenced here
		      $gene_components => [],
		      $sequence_data => [],
		      $locus_info => [],
		      $title_info =>[],
		      $affiliation =>[],
		      $chromosome =>[]
		      );

    # start tree parsing from the tip
    # all branches are retrieved for post-processing.
    &Xml_parser_facilitator::tree_parser (@{$tree}, \%gene_data);


    ## ------------ Pt. 1 ------------------
    ## Process the gene component branches
    my $x = 1;
    my @allcomponents;
    foreach my $ref (@{$gene_data{$gene_components}}) {
	print "$x\n" if ($SEE);
	$x++;

        # traverse the branch
	my $obj = { "Gene-ref_locus" => 0, #pub_locus
		    "Gene-ref_syn_E" => 0, #locus
		    "Seq-feat_comment" => 0, #pub_comment
		    "Seq-interval_from" => [],
		    "Seq-interval_to" => [], 
		    "RNA-ref_type" => 0, #flag to identify mRNA component (all mRNA exons)
		    "Seq-id_gi" => 0,
		    "Prot-ref_name_E" => 0, #flag to identify CDSinfo (com_names, etc.)
		    "Cdregion" => 0, #flag to identify CDS component (all CDS exons)
		    "Na-strand" => 0, #nucleic acid strand (default = "plus")
		    "Prot-ref_desc" => 0, #additional pub_comment information.
		    "Seq-feat_pseudo" => 0  #indicates pseudogene.
		    };
	&Xml_parser_facilitator::branch_parser ("start", $ref, $obj);
	push (@allcomponents, $obj);
	## print out parsed gene component data
	if ($SEE) {
	    print "--------------------Parsing gene-related info from XML---------------\n";
	    foreach my $key (keys %$obj) {
		if (! ref ($obj->{$key})) {
		    if ($obj->{$key}) {
			print "$key\t$obj->{$key} (not ref)\n";
		    }
		} else {
		    if (ref ($obj->{$key}) eq "ARRAY") {
			print "$key\t@{$obj->{$key}} (array ref)\n";
		    } elsif (ref ($obj->{$key}) eq "HASH") {
			my @stuff = %{$obj->{$key}};
			print "$key\t@stuff (hash ref)\n";
		    }
		}
	    }
	    print "\n";
	}
    }

    # separate all gene components into categories and then link them.
    # use CDS's as the primary keys because they're the most informative and most
    # probable to be present to represent the gene.

    my (@CDSs, @mRNAs, @genes, @CDS_info); #store the separated gene component references.

    ## Newly created components have the following data structure
    #$new_obj = { DIR => <+|->,
    #             F_SPAN => [end5, end3],, #only relative to forward strand
    #             COORDS => { end5 => end3, ....},
    #             locus  => "locus", #only for gene
    #             pub_locus => "pub_locus", #only for gene
    #             com_name => "com_name", #for CDS_info
    #             pub_comment => "pub_comment", #for gene
    #             gi => "gi-A / gi-B ....", #could use to track features.
    #             type => WHAT AM I
    #           }


    my @attributes = qw (DIR F_SPAN COORDS locus com_name pub_comment gi type);


    foreach my $component (@allcomponents) {
	if ($component->{"Cdregion"}) {
	    push (@CDSs, &create_exon_component($component));
	}  elsif ($component->{"Prot-ref_name_E"} || $component->{"Prot-ref_desc"}) {
	    push (@CDS_info, &create_CDSinfo_component($component));
	}  elsif ($component->{"RNA-ref_type"}) {
	    push (@mRNAs, &create_exon_component($component));
	}  elsif ($component->{"Gene-ref_locus"} || $component->{"Gene-ref_syn_E"}) { # gene
	    push (@genes, &create_gene_component($component));
	}
    }

    # Here, just print out the obj info to make sure I'm doing this right.


    if ($SEE) {
	print "\n----------------Separate Data Structs------------------\n";
	my @objects;
	push (@objects, "GENES", \@genes, "mRNAs", \@mRNAs, "CDSs", \@CDSs, "CDS_info", \@CDS_info);

	foreach my $object (@objects) {
	    if (!ref($object)) {
		print "\n\nPARSING $object !!!!!!!!!!!!!!!\n";
		next;
	    }
	    foreach my $obj (@$object) {
		print "\n";
		foreach my $att (@attributes) {
		    if (exists ($obj->{$att})) {
			if (ref ($obj->{$att}) eq "ARRAY") {
			    print "$att\t@{$obj->{$att}}\n";
			} elsif (ref ($obj->{$att}) eq "HASH") {
			    my @temp = %{$obj->{$att}};
			    print "$att\t@temp\n";
			} else {
			    print "$att\t$obj->{$att}\n";
			}
		    }
		}
	    }
	}
    }
    ## create Gene_obj based on CDS coordinates
    ## look for overlapping gene and mRNA components
    ## store info in Gene_obj's
    ## after identifying overlapping object, add "DELETE" flag to obj.

    my @gene_storage;


    foreach my $CDS_obj (@CDSs) {
	my $curr_Gene_obj = new Gene_obj ();
	my $CDS_dir = $CDS_obj->{"DIR"}; #if a gene is found, use gene direction (more reliable)
	my $CDS_gi = $CDS_obj->{"gi"};
	$curr_Gene_obj->{"strand"} = $CDS_dir;
	$curr_Gene_obj->{"is_pseudogene"} = ($CDS_obj->{"is_pseudogene"});
	## pubcomment info stored in many diff places.  Append or create.
	if ($curr_Gene_obj->{"pub_comment"}) {
	    unless ($CDS_obj->{"pub_comment"} == 0) {
		$curr_Gene_obj->{"pub_comment"} .= $CDS_obj->{"pub_comment"};
	    }
	} else {
	    $curr_Gene_obj->{"pub_comment"} = $CDS_obj->{"pub_comment"};
	}
	my ($CDS_span_left, $CDS_span_right) = @{$CDS_obj->{"F_SPAN"}};
	$curr_Gene_obj->{"CDS_coords"} = $CDS_obj->{"COORDS"};
	
	## find overlapping gene
	my $found = &find_overlapping_element ($CDS_span_left, $CDS_span_right, $CDS_dir, $CDS_gi, \@genes);
	if ($found) {
	    print "**** Found a gene, direction is $found->{DIR}\n" if $SEE;
	    $curr_Gene_obj->{"locus"} = $found->{"locus"};
	    $curr_Gene_obj->{"pub_locus"} = $found->{"pub_locus"};
	    if ($curr_Gene_obj->{"pub_comment"}) {
		unless ($found->{"pub_comment"} == 0) {
		    $curr_Gene_obj->{"pub_comment"} .= $found->{"pub_comment"};
		}
	    } else {
		$curr_Gene_obj->{"pub_comment"} = $found->{"pub_comment"};
	    }
	    $curr_Gene_obj->{"strand"} = $found->{"DIR"};
	}
	
	## find overlapping mRNA
	my $found = &find_overlapping_element ($CDS_span_left, $CDS_span_right, $CDS_dir, $CDS_gi, \@mRNAs);
	if ($found) {
	    $curr_Gene_obj->{"mRNA_coords"} = $found->{"COORDS"};
	} else {
	    # use CDS coords for mRNA if mRNA isn't specified.
	    $curr_Gene_obj->{"mRNA_coords"} = $CDS_obj->{"COORDS"};
	}
	
	## gather CDS_info
	my  $found = &find_overlapping_element ($CDS_span_left, $CDS_span_right, $CDS_dir, $CDS_gi, \@CDS_info);
	if ($found) {
	    $curr_Gene_obj->{"com_name"} = $found->{"com_name"};
	    $curr_Gene_obj->{"gi"} = $found->{"gi"};
	    if ($curr_Gene_obj->{"pub_comment"}) {
		unless ($found->{"pub_comment"} == 0) {
		    $curr_Gene_obj->{"pub_comment"} .= $found->{"pub_comment"};
		}
	    } else {
		$curr_Gene_obj->{"pub_comment"} = $found->{"pub_comment"};
	    }
	}
	
	## add gene to gene_storage
	$curr_Gene_obj->refine_gene_object();
	push (@gene_storage, $curr_Gene_obj);
	
    }
    
    if ($SEE) {
	print "\n------------Gene Objects----------------\n";
	foreach my $gene (@gene_storage) {
	    print "\n\nGENE\n";
	    print $gene->toString();
	}
    }
    
    ## store Gene_obj references as an instance member component.
    $self->{genes} = \@gene_storage;


    ####################################################
    #-------Pt 2: Parse Sequence Data ------------------
    ####################################################

    my $x = 1;
    foreach my $ref (@{$gene_data{$sequence_data}}) {
	print "\nFound Seq-data component $x\n" if ($SEE);
	$x++;

        # traverse the branch
	my $obj = { "IUPACna" => 0
		    };
	&Xml_parser_facilitator::branch_parser ("start", $ref, $obj);
	
	if ($obj->{IUPACna}) {
	    $self->{sequence} = $obj->{IUPACna};
	    last;
	}
    }
    
    #######################################################
    #-------Pt 3: Parse Locus Information -----------------
    #######################################################

    my $obj = { "Textseq-id_name" => 0,      #locus for genbank record
		"Textseq-id_accession" => 0, #genbank accession for record
		"Textseq-id_version" => 0    #apparently version info
		};
    
    ## Info should be in first entry.
    foreach my $ref (@{$gene_data{$locus_info}}) {
	
	&Xml_parser_facilitator::branch_parser ("start", $ref, $obj);
	
	my ($locus, $accession, $version) = ($obj->{"Textseq-id_name"},
					     $obj->{"Textseq-id_accession"},
					     $obj->{"Textseq-id_version"});
	if ($accession) {
	    $self->{locus} = $locus;
	    $self->{accession} = $accession;
	    $self->{version} = $version;
	    last;
	}
    }

    ###################################################
    #-------Pt 4: Parse title information ------------
    ##################################################
   
    my $obj = { "Seqdesc_title" => 0};
    foreach my $ref (@{$gene_data{$title_info}}) {
	&Xml_parser_facilitator::branch_parser ("start", $ref, $obj);
	my ($title) = $obj->{"Seqdesc_title"};
	if ($title) {
	    $self->{title} = $title;
	    last;
	}
    }

    ###################################################
    #-------Pt 5: Parse update_date  information -----
    ##################################################

   
    # higher node already retrieved within title_info.
    # must retree parse this section to gather the update_date_info node.
    my %date_parser = ( $update_date_info => [] );
    #search thru all title_info nodes for update_date_info node.
    foreach my $date_ref (@{$gene_data{$title_info}}) {
	&Xml_parser_facilitator::tree_parser ("start", $date_ref, \%date_parser);
    }

    my $obj = { "Date-std_year" => 0,
		"Date-std_month" => 0,
		"Date-std_day" => 0
		};
    
    foreach my $ref (@{$date_parser{$update_date_info}}) {
	&Xml_parser_facilitator::branch_parser ("start", $ref, $obj);
	if ($obj->{"Date-std_year"}) {
	    my ($date) = $obj->{"Date-std_month"} . "-" . $obj->{"Date-std_day"} . "-" . $obj->{"Date-std_year"};
	    $self->{update_date} = $date;
	    last;
	}
    }


    ####################################################
    #-------- Pt 6: Parse Affiliate Information -------
    ####################################################
    
    my $obj = { "Affil_std_affil" => 0,
	    "Affil_str" => 0};
    my $affil;
    foreach my $ref (@{$gene_data{$affiliation}}) {
	&Xml_parser_facilitator::branch_parser ("start", $ref, $obj);
	foreach my $key (%$obj) {
	    if ($obj->{$key}) {
		$affil .= $obj->{$key} . "\n\n";
	    }
	}
    }
    if ($affil) {
	$self->{affiliation} = $affil;
    }

    ###################################################
    #----------- Pt 7: Parse Chromosome Information ---
    ###################################################

    my $obj = { "SubSource_name" => 0};
    foreach my $ref (@{$gene_data{$chromosome}}) {
	&Xml_parser_facilitator::branch_parser ("SubSource_name", $ref, $obj);
	if ($obj->{SubSource_name}) {
	    my $chromosome = $obj->{SubSource_name};
	    $self->{chromosome} = $chromosome;
	    last;
	}
    }
    



    #END OF TREE PARSING


    #########################################################################
    # Finess the genes (shift pub_comment to com_name if no com_name specified 
    #########################################################################

    my @genes = &get_genes($self);
    foreach my $gene (@genes) {
	if (!$gene->{'com_name'} && $gene->{'pub_comment'}) {
	    $gene->{'com_name'} = $gene->{'pub_comment'};
	    $gene->{'pub_comment'} = 0;
	}
    }


    ## return Genbankxmlparser object
    bless($self);
    return ($self);
}

sub get_genes {
    my $self = shift;
    return (@{$self->{genes}});
}

sub get_sequence {
    my $self = shift;
    return ($self->{sequence});
}


####
sub find_overlapping_element {
    my ($input_left, $input_right, $input_dir, $input_gi, $collect_ref) = @_;
    foreach my $compare_obj (@{$collect_ref}) {
	if ($compare_obj->{"DELETE"}) { next;} #ignore previously chosen ones.
	my ($temp_left, $temp_right) = @{$compare_obj->{"F_SPAN"}};
	my ($temp_dir, $found_gi, $overlap);
	if ($temp_left && $temp_right) {
	    $temp_dir = $compare_obj->{"DIR"};
	    #test overlap
	    print "$input_left < $temp_right) && ($input_right > $temp_left\t" if ($SEE);
	    if (($input_left < $temp_right) && ($input_right > $temp_left)) {
		$overlap = 1;
		print "YES\n" if ($SEE);
	    } else {
		print "NO\n" if ($SEE);
	    }
	} else {
	    #no coords, so must examine gi numbers for matches.
	    my $temp_gi = $compare_obj->{"gi"};
	    if ($input_gi =~ /\b$temp_gi\b/) {
		$found_gi = 1;
	    }
	}

	## Determine if match found
	if ( (($input_dir eq $temp_dir) && $overlap) || $found_gi) {
	    $compare_obj->{"DELETE"} = 1;
	    return ($compare_obj);
	}
    }
    return (0);
}


####
sub create_CDSinfo_component {
    my ($obj) = @_;
    my $ret = {};
#    my $info = $obj->{"Prot-ref_name_E"};
#    if ($info =~ /^(.*);\s+(\d+)-(\d+)\s*$/) { #try reg-exp grabbing
#	$ret->{"com_name"} = $1;
#	my @coords = ($2, $3);
#	$ret->{"F_SPAN"} = &get_forward_span (@coords);
#	$ret->{"DIR"} = &determine_orientation(@coords);
#    } else {
    $ret->{"com_name"} = $obj->{"Prot-ref_name_E"};
    $ret->{"pub_comment"} = $obj->{"Prot-ref_desc"};
    $ret->{"F_SPAN"} = []; #empty ...no useful coords here.
#    }
    $ret->{"type"} = "CDSinfo";
    $ret->{"gi"} = $obj->{"Seq-id_gi"};
    return ($ret);
}


####
sub create_exon_component {
    my ($obj) = @_;
    my $ret = {};
    my $coord1 = $obj->{"Seq-interval_from"};
    my $coord2 = $obj->{"Seq-interval_to"};
    $ret->{"DIR"} = ($obj->{"Na-strand_atts"}->{"value"} eq "minus") ? '-' : '+';
    $ret->{"is_pseudogene"} = ($obj->{"Seq-feat_pseudo_atts"}->{"value"} eq "true") ? 1:0;
    my @coords;
    for (my $i = 0; $i <= $#{$coord1}; $i++) {
	push (@coords, ++$coord1->[$i], ++$coord2->[$i]);
	if ($ret->{"DIR"} eq '+') {
	    $ret->{"COORDS"}->{$coord1->[$i]} = $coord2->[$i];
	} else {
	    $ret->{"COORDS"}->{$coord2->[$i]} = $coord1->[$i];
	}  
    }
    $ret->{"F_SPAN"} = &get_forward_span (@coords);
    $ret->{"type"} = "EXONS (CDS's or mRNA's)";
    $ret->{"gi"} = $obj->{"Seq-id_gi"};
    $ret->{"pub_comment"} = $obj->{"Seq-feat_comment"};
    return ($ret);
}

####
sub create_gene_component {
    my ($obj) = @_;
    my $ret = {};
    if ($obj->{"Gene-ref_locus"}) {
	$ret->{"pub_locus"} = $obj->{"Gene-ref_locus"};
    }
    if ($obj->{"Gene-ref_syn_E"}) {
	$ret->{"locus"} = $obj->{"Gene-ref_syn_E"};
    } else {
	$ret->{"locus"} = $ret->{"pub_locus"};
    }
    $ret->{"pub_comment"} = $obj->{"Seq-feat_comment"};
    $ret->{"DIR"} = ($obj->{"Na-strand_atts"}->{"value"} eq "minus") ? '-' : '+';
    my $coord1 = $obj->{"Seq-interval_from"}->[0];
    my $coord2 = $obj->{"Seq-interval_to"}->[0];
    if ($ret->{"DIR"} eq '+') {    
	$ret->{"COORDS"}->{$coord1} = ++$coord2;
    } else {
	$ret->{"COORDS"}->{$coord2} = ++$coord1;
    }
    $ret->{"F_SPAN"} = &get_forward_span ($coord1, $coord2);
    $ret->{"type"} = "GENE";
    $ret->{"gi"} = $obj->{"Seq-id_gi"};
    return ($ret);
}


####
sub get_forward_span {
    my (@coords) = @_;
    @coords = sort {$a<=>$b} @coords;
    my @ret = ($coords[0], $coords[$#coords]);
    return (\@ret);
}


####
sub determine_orientation {
    my ($end5, $end3) = @_;
    my $direction;
    if ($end5 < $end3) {
	$direction = "+";
    } elsif ($end5 > $end3) {
	$direction = "-";
    }
    return ($direction);
}

sub toString {
    my ($self) = shift;
    my $output;

    # print submission information
    my ($locus, $accession, $version, $title, $update_date, $affiliation, $chromosome) = ($self->{locus},
											  $self->{accession},
											  $self->{version},
											  $self->{title},
											  $self->{update_date},
											  $self->{affiliation},
											  $self->{chromosome});
    
    $output .= "Locus: $locus\nAccession: $accession\nVersion: $version\nTitle: $title\nUpdate date: $update_date\naffiliation: $affiliation\n"
	. "chromosome: $chromosome\n";

    #print gene information
    my @genes = $self->get_genes();
    my $x = 0;
    foreach my $gene (@genes) {
	$x++;
	$output .=  "\nGENE $x\n";
	$output .=  $gene->toString();
    }

    #print submission's sequence.
    my $seq = $self->get_sequence();    
    #$output .=  "SEQUENCE:\n$seq\n";
    return ($output);
}

1;



















