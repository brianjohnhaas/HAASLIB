#!/usr/local/bin/perl

our ($SEE, $DEBUG);

package Annot_db_loader;

use strict;
use lib ("$ENV{EGC_SCRIPTS}");
use Egc_library;
use Data::Dumper;
use Storable;
use Feat_name_generator;

## if sequence is given in constructor, the sequence is loaded and the asmbl_id is set.
## otherwise, asmbl_id must be set explicitly.

sub new {
    shift;
    my ($dbproc, $sequence_ref_or_asmbl_id) = @_;  ## If reference to a sequence, loads sequence and sets asmbl_id.
    
    my ($asmbl_id, $sequence_ref);
    if (ref $sequence_ref_or_asmbl_id eq "SCALAR") {
        $sequence_ref = $sequence_ref_or_asmbl_id;
    } elsif ($sequence_ref_or_asmbl_id) {
        $asmbl_id = $sequence_ref_or_asmbl_id;
        if (length ($asmbl_id) > 6) {
            die "Sorry, $asmbl_id is too large to be an asmbl_id.  Please check your parameters to this constructor.\n";
        }
    }
    
    my $self = { asmbl_id=>0,
                 dbproc=>0,
                 failed_genes=>[] #store list of gene_obj's which failed the loading process.
                 };
    bless $self;
    if ($sequence_ref && !$asmbl_id) {
        $self->load_seq_get_asmbl_id ($dbproc, $sequence_ref);
    } elsif ($dbproc && $asmbl_id) {
        $self->set_dbproc_asmbl_id ($dbproc, $asmbl_id);	    
    }
    return ($self);
}


sub reset_failed_gene_list {
    my $self = shift;
    $self->{failed_genes} = [];
}

sub get_failed_genes {
    ## Return a list of genes that failed
    my $self = shift;
    return (@{$self->{failed_genes}});
}

sub store_failed_gene {
    my $self = shift;
    my $gene = shift;
    if (ref $gene) {
        push (@{$self->{failed_genes}}, $gene);
        return (1); #succeeded.
    } else {
        return(0); #failed
    }
} 

sub load_seq_get_asmbl_id {
    my $self = shift;
    my $dbproc = shift;
    my $seqref = shift;
    my $asmbl_id;
    if ($DEBUG) {
        $asmbl_id = "ASMBL_ID";
    } else {
        $asmbl_id = &load_seq_assembly($dbproc, $seqref);
    }
    if ($asmbl_id) {
        print STDERR "Created New Assembly: $asmbl_id\n";
    } else {
        die "Couldn't insert asmbl into assembly\n";
    }
    ## set asmbl_id, dbproc
    $self->set_dbproc_asmbl_id ($dbproc, $asmbl_id);
}

sub set_dbproc_asmbl_id {
    my $self = shift;
    my $dbproc = shift;
    my $asmbl_id = shift;
    $self->{asmbl_id} = $asmbl_id;
    $self->{dbproc} = $dbproc;
}

sub get_asmbl_id {
    my $self = shift;
    return ($self->{asmbl_id});
}


sub load_genes {
    print STDERR "Loading Genes\n";
    my $self = shift;
    my (@genes) = @_;
    my $asmbl_id = $self->{asmbl_id};
    my $dbproc = $self->{dbproc};
    
    my @nonProteinCodingGenes; ## process later.
    
    $dbproc->{AutoCommit} = 0; #set autocommit off so transactions can be used.
    
    my $x = 0;
    my (@TU_feat_names); # to be returned as a list of each gene loaded into the db.

    ## prepare for feat_name generation:

    my $next_TU_feat_name = &getNextName ($dbproc, $asmbl_id, "TU");
    my $next_model_feat_name = &getNextName($dbproc, $asmbl_id, "model");
    my $next_exon_feat_name = &getNextName($dbproc, $asmbl_id, "exon");
    my $next_cds_feat_name = &getNextName($dbproc, $asmbl_id, "CDS");

    if ($SEE) {
        print "Next names:\n"
            . "\tTU: $next_TU_feat_name\n"
            . "\tmodel: $next_model_feat_name\n"
            . "\texon: $next_exon_feat_name\n"
            . "\tCDS: $next_cds_feat_name\n\n";
    }
    
    my $feat_name_generator = new Feat_name_generator($next_TU_feat_name, $next_model_feat_name, $next_exon_feat_name, $next_cds_feat_name);
    
    
    foreach my $gene (@genes) {
        ## initialize the feat_names which will be set upon successful db loading:
        $gene->{TU_feat_name} = undef;
        $gene->{Model_feat_name} = undef;
        
        unless ($gene->{gene_type} eq "protein-coding") {
            ## only loading protein coding (or pseudogene) genes right now.
            push (@nonProteinCodingGenes, $gene);
            next;
        } 
        if ($SEE) {
            print $gene->toString();
        }
                
        my ($TU_feat_name,$model_feat_name);

        eval {
            $gene->refine_gene_object();
            my @exons = $gene->get_exons();
            
            unless (@exons) {
                die "No exons listed for gene " . $gene->toString();
            }
            
            ### TU RELATED PROCESS
            
            $TU_feat_name = ($DEBUG) ? ('NextTU') : $feat_name_generator->next_TU_feat_name();
            $self->process_functional_annotations($gene, $TU_feat_name);
            
            my ($gene_span, $model_span) = ($gene->{gene_span}, $gene->{model_span});
            my ($curated_com_name, $curated_gene_structure) = ($gene->{curated_com_name},
                                                               $gene->{curated_gene_structure});
            
            ############################################################################################################
            ######## Create Gene Structure #############################################################################
            ############################################################################################################
            
            ## Create TU
            &insert_asm_feature($dbproc, 
                                $TU_feat_name, 
                                $asmbl_id, 
                                "TU", 
                                $gene_span->[0], 
                                $gene_span->[1], 
                                "annot_db_loader", 
                                "annot_db_loader", 
                                1, 0);
            
            ## reset gene's TU_feat_name
            $gene->{TU_feat_name} = $TU_feat_name;
            
            if ($curated_com_name) {
                my $query = "update asm_feature set curated = 1 where feat_name = \"$TU_feat_name\"";
                &RunMod($dbproc, $query);
            }
            
            
            ## MODEL RELATED PROCESS
            $model_feat_name = ($DEBUG) ? ('NextModel') : $feat_name_generator->next_model_feat_name();
            &insert_asm_feature($dbproc, 
                                $model_feat_name, 
                                $asmbl_id, 
                                "model", 
                                $model_span->[0], 
                                $model_span->[1], 
                                "annot_db_loader", 
                                "annot_db_loader", 
                                1, 0);
            
            &insert_phys_ev($dbproc, $model_feat_name, "working", "annot_db_loader");
            &insert_feat_link($dbproc, $model_feat_name, $TU_feat_name, "annot_db_loader","now");
            
            ## reset gene's model_feat_name
            $gene->{Model_feat_name} = $model_feat_name;
            
            if ($curated_gene_structure) {
                my $query = "update asm_feature set curated = 1 where feat_name = \"$model_feat_name\"";
                &RunMod($dbproc, $query);
            }
            
            if (my $model_pub_locus = $gene->{model_pub_locus}) {
                my $query = qq { insert ident (feat_name, com_name, pub_locus) values ("$model_feat_name", "", "$model_pub_locus") };
                &RunMod($dbproc, $query);
            }
            
            
            ## EXON AND CDS RELATED PROCESS
            foreach my $exon (@exons) {
                my $exon_feat_name = ($DEBUG) ? ('NextExon') : $feat_name_generator->next_exon_feat_name();
                my ($end5, $end3) = ($exon->{end5}, $exon->{end3});
                &insert_asm_feature($dbproc, 
                                    $exon_feat_name, 
                                    $asmbl_id, 
                                    "exon", 
                                    $end5, 
                                    $end3, 
                                    "annot_db_loader", 
                                    "annot_db_loader", 
                                    1, 0);
                &insert_phys_ev($dbproc, $exon_feat_name, "working", "annot_db_loader");
                &insert_feat_link($dbproc, $exon_feat_name, $model_feat_name, "annot_db_loader","now");
                
                my $cds = $exon->{CDS_exon_obj};
                if ($cds) {
                    my $cds_feat_name = ($DEBUG) ? ('NextCDS') : $feat_name_generator->next_cds_feat_name();
                    my ($end5, $end3) = ($cds->{end5}, $cds->{end3});
                    &insert_asm_feature($dbproc, 
                                        $cds_feat_name, 
                                        $asmbl_id, 
                                        "CDS", 
                                        $end5, 
                                        $end3, 
                                        "annot_db_loader", 
                                        "annot_db_loader", 
                                        1, 0);
                    &insert_phys_ev($dbproc, $cds_feat_name, "working", "annot_db_loader");
                    &insert_feat_link($dbproc, $cds_feat_name, $exon_feat_name, "annot_db_loader","now");
                }
            }
        };

        ## End transaction
        if (! $@) {
            push (@TU_feat_names, [$TU_feat_name, $model_feat_name]);
            $dbproc->commit;
            ## load alt splice isoforms if available:
            if (my @isoforms = $gene->get_additional_isoforms()) {
                foreach my $isoform (@isoforms) {
                    print STDERR "loading isoform.\n";
                    $self->load_alt_splice_isoform($TU_feat_name, $isoform, $feat_name_generator);
                    ## alt splice loading turns autocommit back on.
                    $dbproc->{AutoCommit} = 0; #set autocommit off again so transactions can be used.
                }
            }
            
            print STDERR "loaded $TU_feat_name, $model_feat_name\n";
            
        } else {
            $dbproc->rollback;
            # reset TU and Model feats to undef:
            $gene->{TU_feat_name} = undef;
            $gene->{Model_feat_name} = undef;
            
            $self->store_failed_gene($gene);
            print STDERR "Failed loading " . $gene->toString() . "\nStored for further use.\n";
        }
        
        
        
    }
    
    ## reset transaction behavior to autocommit:
    $dbproc->{AutoCommit} = 1;
    
    
    ## Write failed genes to a file:
    my @failed_genes = $self->get_failed_genes();
    if (@failed_genes) {
        my $persistent_objs_file = "failed_gene_loads.$$.stobj";
        &Storable::lock_nstore(\@failed_genes, $persistent_objs_file);
        open (FILE, ">failed_gene_loads.$$.txt");
        my $counter = 0;
        foreach my $failed_gene (@failed_genes) {
            $counter++;
            print FILE "Failed gene $counter:\n";
            print FILE $failed_gene->toString() . "\n\n";
        }
        close FILE;
        
        print "****\nThe loading of $counter genes have failed.  Please find their description in file \'failed_gene_loads.$$.txt\', and recoverable Gene_obj's for them can be found in \'$persistent_objs_file\' recoverable using Storable::lock_retrieve()\n****\n";
        
    }
    
    if (@nonProteinCodingGenes) {
        $self->load_non_protein_coding_genes(@nonProteinCodingGenes);
    }
    
}


sub load_clone_info {
    print STDERR "Loading Clone_info Table\n";
    my $self = shift;
    my $dbproc = $self->{dbproc};
    my $asmbl_id = $self->{asmbl_id};
    my ($clone_id, $clone_name, $seq_group, $assignby, $chromo) = @_;
    
    unless ($chromo) { $chromo = 0;}
    
    unless ($clone_id =~ /^\d+$/) {
        $clone_id = 0;
    }
    my $query = "insert clone_info (asmbl_id, clone_id, clone_name, seq_group, orig_annotation, tigr_annotation, status, length, assignby, date, is_public, chromo) values ($asmbl_id, $clone_id, \"$clone_name\", \"$seq_group\", 1,0, \"Annotation\", 0, \"$assignby\", getdate(), 0, $chromo)\n";
    if ($DEBUG) {
        print "query: $query\n";
    } else {
        &do_sql ($dbproc, $query);
    }
    my $query = "update clone_info set length = (select datalength(sequence) from assembly where asmbl_id = $asmbl_id) where asmbl_id = $asmbl_id\n";
    if ($DEBUG) {
        print "query: $query\n";
    } else {
        &do_sql ($dbproc, $query);
    }
}


###############################
## updating a single gene structure:
####

sub update_curr_gene {
    my $self = shift;
    my ($curr_gene, $new_gene) = @_;
    my $dbproc = $self->{dbproc};
    ### delete old exons/cds's
    
    $dbproc->{AutoCommit} = 0; #set autocommit off so transactions can be used.
    
    my $asmbl_id = $curr_gene->{asmbl_id};
    
    eval {
        
        my $curr_model_feat_name = $curr_gene->{Model_feat_name};
        my @curr_gene_exons = $curr_gene->get_exons();
        foreach my $exon (@curr_gene_exons) {
            my $exon_feat_name = $exon->{feat_name};
            if (my $cds_obj = $exon->get_CDS_obj()) {
                my $cds_feat_name = $cds_obj->{feat_name};
                &delete_this_cds($cds_feat_name, $dbproc, $exon_feat_name);
            }
            &delete_this_exon($exon_feat_name, $dbproc, $curr_model_feat_name);
        }
        
        ## Load new features
        
        ## Gene components first:
        foreach my $exon ($new_gene->get_exons()) {
            my $exon_feat_name = &getNextName ($dbproc, $asmbl_id, "exon");
            my ($end5, $end3) = ($exon->{end5}, $exon->{end3});
            &insert_asm_feature($dbproc, 
                                $exon_feat_name, 
                                $asmbl_id, 
                                "exon", 
                                $end5, 
                                $end3, 
                                "annot_db_loader", 
                                "annot_db_loader", 
                                0, 0);
            &insert_phys_ev($dbproc, $exon_feat_name, "working", "annot_db_loader");
            &insert_feat_link($dbproc, $exon_feat_name, $curr_model_feat_name, "annot_db_loader","now");
            
            my $cds = $exon->{CDS_exon_obj};
            if ($cds) {
                my $cds_feat_name = &getNextName ($dbproc, $asmbl_id, "CDS");
                my ($end5, $end3) = ($cds->{end5}, $cds->{end3});
                &insert_asm_feature($dbproc, 
                                    $cds_feat_name, 
                                    $asmbl_id, 
                                    "CDS", 
                                    $end5, 
                                    $end3, 
                                    "annot_db_loader", 
                                    "annot_db_loader", 
                                    0, 0);
                &insert_phys_ev($dbproc, $cds_feat_name, "working", "annot_db_loader");
                &insert_feat_link($dbproc, $cds_feat_name, $exon_feat_name, "annot_db_loader","now");
            }
        }
        ## End transaction
    };
    
    
    if (! $@) {
        $dbproc->commit;
        my $curr_tu_feat_name = $curr_gene->{TU_feat_name};
        &update_gene_attributes($dbproc, $asmbl_id, $curr_tu_feat_name);
        # update obj atts
        $new_gene->{TU_feat_name} = $curr_gene->{TU_feat_name};
        $new_gene->{Model_feat_name} = $curr_gene->{Model_feat_name};
        
        $dbproc->commit;
    } else {
        $dbproc->rollback;
        $self->store_failed_gene($new_gene);
        print STDERR "Failed loading " . $new_gene->toString() . "\nStored for further use.\n";
    }
    ## reset transaction behavior to autocommit:
    $dbproc->{AutoCommit} = 1;
}

sub load_alt_splice_isoform {
    my $self = shift;
    my ($TU_feat_name, $new_gene_obj, $feat_name_generator) = @_; #feat_name_generator is optional
    my $dbproc = $self->{dbproc};
    my $asmbl_id = $self->{asmbl_id};
    
    $dbproc->{AutoCommit} = 0; #set autocommit off so transactions can be used.
        
    eval {
        ## MODEL RELATED PROCESS
        my $model_feat_name = ($DEBUG) ? ('NextModel') : 
            ($feat_name_generator) 
            ? $feat_name_generator->next_model_feat_name()
            : &getNextName ($dbproc, $asmbl_id, "model");
        
        $new_gene_obj->{TU_feat_name} = $TU_feat_name;
        $new_gene_obj->{Model_feat_name} = $model_feat_name;
        my $model_span = $new_gene_obj->{model_span};
        &insert_asm_feature($dbproc, 
                            $model_feat_name, 
                            $asmbl_id, 
                            "model", 
                            $model_span->[0], 
                            $model_span->[1], 
                            "annot_db_loader", 
                            "annot_db_loader", 
                            1, 0);
        &insert_phys_ev($dbproc, $model_feat_name, "working", "annot_db_loader");
        &insert_feat_link($dbproc, $model_feat_name, $TU_feat_name, "annot_db_loader","now");
        
        ## reset gene's model_feat_name
        $new_gene_obj->{Model_feat_name} = $model_feat_name;
        
        if (my $model_pub_locus = $new_gene_obj->{model_pub_locus}) {
            my $query = qq { insert ident (feat_name, com_name, pub_locus) values ("$model_feat_name", "", "$model_pub_locus") };
            &RunMod($dbproc, $query);
        }
        
        
        my @exons = $new_gene_obj->get_exons();
        ## EXON AND CDS RELATED PROCESS
        foreach my $exon (@exons) {
            my $exon_feat_name = ($DEBUG) ? ('NextExon') : ($feat_name_generator) 
                ? $feat_name_generator->next_exon_feat_name()
                : &getNextName ($dbproc, $asmbl_id, "exon");
            
            my ($end5, $end3) = ($exon->{end5}, $exon->{end3});
            &insert_asm_feature($dbproc, 
                                $exon_feat_name, 
                                $asmbl_id, 
                                "exon", 
                                $end5, 
                                $end3, 
                                "annot_db_loader", 
                                "annot_db_loader", 
                                1, 0);
            &insert_phys_ev($dbproc, $exon_feat_name, "working", "annot_db_loader");
            &insert_feat_link($dbproc, $exon_feat_name, $model_feat_name, "annot_db_loader","now");
            $exon->{feat_name} = $exon_feat_name;
            my $cds = $exon->{CDS_exon_obj};
            if ($cds) {
                my $cds_feat_name = ($DEBUG) ? ('NextCDS') : ($feat_name_generator)
                    ? $feat_name_generator->next_cds_feat_name()
                    : &getNextName ($dbproc, $asmbl_id, "CDS");
                
                
                my ($end5, $end3) = ($cds->{end5}, $cds->{end3});
                &insert_asm_feature($dbproc, 
                                    $cds_feat_name, 
                                    $asmbl_id, 
                                    "CDS", 
                                    $end5, 
                                    $end3, 
                                    "annot_db_loader", 
                                    "annot_db_loader", 
                                    1, 0);
                &insert_phys_ev($dbproc, $cds_feat_name, "working", "annot_db_loader");
                &insert_feat_link($dbproc, $cds_feat_name, $exon_feat_name, "annot_db_loader","now");
                $cds->{feat_name} = $cds_feat_name;
            }
        }
        ## End transaction
    };

    
    if (! $@) {
        $dbproc->commit;
        &update_gene_attributes($dbproc, $asmbl_id, $TU_feat_name);
        $dbproc->commit;
        
        
    } else {
        $dbproc->rollback;
        $self->store_failed_gene($new_gene_obj);
        print STDERR "Failed loading " . $new_gene_obj->toString() . "\nStored for further use.\n";
    }
    
    
    ## reset transaction behavior to autocommit:
    $dbproc->{AutoCommit} = 1;
}




#### private methods.
sub load_seq_assembly {
    print STDERR "Loading Assembly Table\n";
    my ($dbproc, $string_ref) = @_;
    my $query = "set textsize 1000000000";
    &RunMod ($dbproc, $query);
    my $query = "insert assembly (sequence, seq#, ed_pn, ed_date) values( \"$$string_ref\", 1, \"0\", getdate())";
    my $asmbl_id = &first_result_sql ($dbproc, $query);
    return ($asmbl_id); # return new asmbl_id
}

sub load_gene_ontology_assignments {
    my ($self, $TU_feat_name, @gene_ontology_objs) = @_;
    my $dbproc = $self->{dbproc};
    my %seen;
    foreach my $gene_ontology_obj (@gene_ontology_objs) {
        
        my $go_id = $gene_ontology_obj->{go_id};
        my $assigned_by = $gene_ontology_obj->{assignby} || 'unknown';
        my $date;
        if ($date = $gene_ontology_obj->{date}) {
            $date = "\"$date\"";
        } else {
            $date = 'getdate()';
        }
        my $qualifier = $gene_ontology_obj->{qualifier};
        if ($qualifier) {
            $qualifier = "\"$qualifier\"";
        } else {
            $qualifier = "NULL";
        }
        
        my $role_link_id;
        
        ## Check to see if an ID already exists for this:
        my $query = "select id from go_role_link where go_id = \"$go_id\" and feat_name = \"$TU_feat_name\" and assigned_by = \"$assigned_by\"";
        $role_link_id = &first_result_sql($dbproc, $query);
        
        unless ($role_link_id) {
            
            my $query = "insert go_role_link (feat_name, go_id, assigned_by, date, qualifier) values (\"$TU_feat_name\", \"$go_id\", \"$assigned_by\", $date, $qualifier)\n";
            &RunMod($dbproc, $query);
            die if $QUERYFAIL;
            
            my $query = 'select @@identity';
            $role_link_id = &first_result_sql($dbproc, $query);
            die unless $role_link_id;
            
        }
        
        my @evidences = $gene_ontology_obj->get_evidence();
        
        foreach my $evidence_ref (@evidences) {
            my $ev_code = $evidence_ref->{ev_code};
            my $evidence = $evidence_ref->{evidence};
            my $with_ev = $evidence_ref->{with_ev};
            
            my $ev_data = join (",", $role_link_id, $ev_code, $evidence, $with_ev);
            
            my $query = "insert go_evidence (role_link_id, ev_code, evidence, with_ev) values ($role_link_id, \"$ev_code\", \"$evidence\", __WITH_EV__)\n";	
            if ($with_ev) {
                $query =~ s/__WITH_EV__/\"$with_ev\"/;
            } else {
                $query =~ s/__WITH_EV__/NULL/;
            }
            
            if ($seen{$ev_data}) { 
                ## already inputed
                next;
            } else {
                $seen{$ev_data}= 1; ## avoid duplicates
            }
            &RunMod($dbproc, $query);
            
        }
    }
}


####
sub load_non_protein_coding_genes {
    my $self = shift;
    my @genes = @_;
    
    ## types of genes to load:
    #            rRNA|snoRNA|snRNA|tRNA
    #   use the {gene_type} attribute to determine how to load the feature.
    #       tRNAs treated specially, non-tRNAs treated identically
    #   for non-tRNAs, use {gene_type} as value of asm_feature.feat_type
    
    foreach my $gene (@genes) {
        my $gene_type = lc $gene->{gene_type};
        if ($gene_type eq "trna") {
            $self->load_tRNAgene($gene);
        } else {
            $self->load_RNAgene($gene);
        }
    }
}



####
sub load_tRNAgene {
    my $self = shift;
    my $gene = shift;
    my $dbproc = $self->{dbproc};
    my $asmbl_id = $self->{asmbl_id};
    
    $dbproc->{AutoCommit} = 0; #set autocommit off so transactions can be used.
    
    my ($pre_trna_feat_name, $trna_feat_name, $rna_exon_feat_name);

    eval {
        my ($end5, $end3) = $gene->get_coords();
                
        ## load pre-trna
        $pre_trna_feat_name = &getNextName($dbproc, $asmbl_id, "pre-tRNA");
        &insert_asm_feature($dbproc, 
                            $pre_trna_feat_name, 
                            $asmbl_id, 
                            "pre-tRNA", 
                            $end5, 
                            $end3, 
                            "annot_db_loader", 
                            "annot_db_loader", 
                            0, 0);
        ## load tRNA
        $trna_feat_name = &getNextName($dbproc, $asmbl_id, "tRNA");
        &insert_asm_feature($dbproc, 
                            $trna_feat_name, 
                            $asmbl_id, 
                            "tRNA", 
                            $end5, 
                            $end3, 
                            "annot_db_loader", 
                            "annot_db_loader", 
                            0, 0);
        ## load rna-exon
        $rna_exon_feat_name = &getNextName($dbproc, $asmbl_id, "rna-exon");
        &insert_asm_feature($dbproc, 
                            $rna_exon_feat_name, 
                            $asmbl_id, 
                            "rna-exon", 
                            $end5, 
                            $end3, 
                            "annot_db_loader", 
                            "annot_db_loader", 
                            0, 0);
        
        ## link pre-tRNA to tRNA:
        &insert_feat_link($dbproc, $trna_feat_name, $pre_trna_feat_name, "annot_db_loader","now");
        
        ## link tRNA to rna-exon:
        &insert_feat_link($dbproc, $rna_exon_feat_name, $trna_feat_name, "annot_db_loader","now");
        
        ## put entry in ident/ident_xref:
        $self->process_functional_annotations($gene, $trna_feat_name);
    };
        
    if (! $@) {
        $dbproc->commit;
        $gene->{TU_feat_name} = $pre_trna_feat_name;
        $gene->{Model_feat_name} = $trna_feat_name;
    } else {
        print STDERR "tRNA gene loading failed: " . $gene->toString();
        $dbproc->rollback;
    }
    
    ## reset autocommit behaviour
    $dbproc->{AutoCommit} = 1;
    
}





####
sub load_RNAgene {
    my $self = shift;
    my $gene = shift;
    my $dbproc = $self->{dbproc};
    my $asmbl_id = $self->{asmbl_id};
    
    $dbproc->{AutoCommit} = 0; #set autocommit off so transactions can be used.
    
    my $feat_name;

    eval {
        my ($end5, $end3) = $gene->get_coords();
        
        my $gene_type = lc $gene->{gene_type};
        $gene_type =~ s/rna/RNA/; 
        ## load asm_feature
        
        $feat_name = &getNextName($dbproc, $asmbl_id, $gene_type);
        
        &insert_asm_feature($dbproc, 
                            $feat_name, 
                            $asmbl_id, 
                            $gene_type, 
                            $end5, 
                            $end3, 
                            "annot_db_loader", 
                            "annot_db_loader", 
                            0, 0);
                
        ## put entry in ident/ident_xref:
        $self->process_functional_annotations($gene, $feat_name);
        
    };
    
    if (! $@) {
        $dbproc->commit;
        $gene->{TU_feat_name} = $feat_name;
    } else {
        print STDERR "tRNA gene loading failed: " . $gene->toString();
        $dbproc->rollback;
    }
    
    ## reset autocommit behaviour
    $dbproc->{AutoCommit} = 1;
    
}


######################################
## ident and ident_xref operations:
######################################

sub process_functional_annotations {
    my $self = shift;
    my $gene = shift;
    my $TU_feat_name = shift;  ## can be an rna feat_name instead.
    
    my $dbproc = $self->{dbproc};
    my $asmbl_id = $self->{asmbl_id};
    
    
    my ($locus, $pub_locus,$com_name, $pub_comment, 
        $is_pseudogene, $alt_locus, 
        $ec_num, $gene_sym, $comment) = ($gene->{locus},
                                         $gene->{pub_locus},
                                         $gene->{com_name},
                                         $gene->{pub_comment},
                                         $gene->{is_pseudogene},
                                         $gene->{alt_locus},
                                         $gene->{ec_num},
                                         $gene->{gene_sym},
                                         $gene->{comment});
    
	if (length($pub_comment) > 255) {
	    $pub_comment = substr ($pub_comment, 0, 255);
	}
	
	unless ($com_name =~ /\w/) {
        
        ## check secondary names:
         if (my @secondary_product_names = $gene->get_secondary_product_names()) {
             $com_name = shift @secondary_product_names;
             # clear secondary product names:
             $gene->{secondary_product_names} = [];
             # set if some left:
             if (@secondary_product_names) {
                 $gene->add_secondary_product_names(@secondary_product_names);
             }
         }
         else {
             $com_name = "!!! ERROR, product name unspecified !!!";
             print stderr "MISSING com_name, so now = $com_name\n";
         }
     }
	
	$is_pseudogene = ($is_pseudogene == 1) ? 1:0;
    
	if ($is_pseudogene && ($com_name !~ /pseudogene/i)) {
	    $com_name .= ", pseudogene";
	}
    
	
	## insert ident info
	
	## Primary functional annotation data, support ident and ident_xref
	
	## Process ComName
	
	my $query = "insert ident (feat_name, locus, com_name, assignby, date, pub_comment, comment, save_history, pub_locus, is_pseudogene, alt_locus, ec\#, gene_sym )\n"
	    . "values (\"$TU_feat_name\", \"$locus\", " . &FormatValue($com_name) . ", \"$0\", getdate(), \"$pub_comment\", \"$comment\",  1, \"$pub_locus\", $is_pseudogene, \"$alt_locus\", \"$ec_num\", \"$gene_sym\")\n";
	
	my $query = "insert ident (feat_name, locus, com_name, assignby, date, pub_comment, comment, save_history, pub_locus, is_pseudogene, alt_locus, ec\#, gene_sym )\n"
	    . "values (?,?,?,?,getdate(),?,?,?,?,?,?,?,?)";
    
	my @values = ($TU_feat_name, $locus, $com_name, $0, $pub_comment, $comment, 1, $pub_locus, $is_pseudogene, $alt_locus, $ec_num, $gene_sym);
	if ($DEBUG) {
	    print "Query: $query\n";
	} else {
	    &do_sql ($dbproc, $query, undef, @values);
    }
	
	## Ident Xref ##
	## Product Names
	&insert_ident_xref($dbproc, $TU_feat_name, $com_name, "product name", "annot_db_loader", "auto", 1);
    
    if (my @secondary_product_names = $gene->get_secondary_product_names()) {
	    foreach my $secondary_product_name (@secondary_product_names) {
            unless ($secondary_product_name =~ /\w/) { next;}
            &insert_ident_xref($dbproc, $TU_feat_name, $secondary_product_name, "product name", "annot_db_loader", "auto", 2);
        }
	}
	
	## Gene names:
	my $gene_name = $gene->{gene_name};
	if ($gene_name =~ /\w/) {
	    &insert_ident_xref($dbproc, $TU_feat_name, $gene_name, "gene name", "annot_db_loader", "auto", 1);
    }
	if (my @secondary_gene_names = $gene->get_secondary_gene_names()) {
	    foreach my $secondary_gene_name (@secondary_gene_names) {
            unless ($secondary_gene_name =~ /\w/) {next;}
            &insert_ident_xref($dbproc, $TU_feat_name, $secondary_gene_name, "gene name", "annot_db_loader", "auto", 2);
        }
	}
	
	
	## EC numbers
	if ($ec_num =~ /\w/) {
	    &insert_ident_xref($dbproc, $TU_feat_name, $ec_num, "ec number", "annot_db_loader", "auto", 1);
    }
	if (my @secondary_ec_numbers = $gene->get_secondary_ec_numbers()) {
	    foreach my $ec_number (@secondary_ec_numbers) {
            unless ($ec_number =~ /\w/) { next;}
            &insert_ident_xref($dbproc, $TU_feat_name, $ec_number, "ec number", "annot_db_loader", "auto", 2);
        }
	}
	
	
	## Gene symbols
	if ($gene_sym =~ /\w/) {
	    &insert_ident_xref($dbproc, $TU_feat_name, $gene_sym, "gene symbol", "annot_db_loader", "auto", 1);
    }
	
    if (my @secondary_gene_symbols = $gene->get_secondary_gene_symbols()) {
	    foreach my $gene_symbol (@secondary_gene_symbols) {
            unless ($gene_symbol =~ /\w/) { next;}
            &insert_ident_xref($dbproc, $TU_feat_name, $gene_symbol, "gene symbol", "annot_db_loader", "auto", 2);
        }
	}
    
}

####
sub load_repeats {
    my $self = shift;
    my @repeats  = @_;
    
    my $dbproc = $self->{dbproc};
    my $asmbl_id = $self->{asmbl_id};
    
    $dbproc->{AutoCommit} = 1;
    
    ## Repeat Structure:
    
    #    my $repeat_struct = { feat_name => $feat_name,
    #			      name => $repeat_type,
    #			      type => 'repeat',
    #			      end5 => $end5,
    #			      end3 => $end3 };
    
    foreach my $repeat (@repeats) {
        my ($desc, $end5, $end3) = ($repeat->{name}, $repeat->{end5}, $repeat->{end3});
        $desc =~ s/\"//g;
        
        my $repeat_feat_name = ($DEBUG) ? ('NextRepeat') : &getNextName ($dbproc, $asmbl_id, "repeat");
        &insert_asm_feature($dbproc, 
                            $repeat_feat_name, 
                            $asmbl_id, 
                            "repeat", 
                            $end5,
                            $end3,
                            "annot_db_loader", 
                            "annot_db_loader", 
                            1, 0);
        
        my $query = "insert ORF_attribute (feat_name, att_type, curated, method, date, assignby, score, score_desc) values (\"$repeat_feat_name\", 'repeat', 0, \"annot_db_loader\", getdate(), \"egc\", \"$desc\", 'repeat description')";
        &RunMod($dbproc, $query);
        
    }
}


1; #end of Annot_db_loader













