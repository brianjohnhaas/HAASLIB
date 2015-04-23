#!/usr/local/bin/perl

package Annot_prediction_loader;

use strict;
use lib ("$ENV{EGC_SCRIPTS}");
use Egc_library;

use vars qw ($DEBUG);
$DEBUG = 0;


## if sequence is given in constructor, the sequence is loaded and the asmbl_id is set.
## otherwise, asmbl_id must be set explicitly.

sub new {
    shift;
    my ($dbproc, $asmbl_id,$ev_type) = @_; #sequence,asmbl_id, update_toggle fields optional
    my $self = {asmbl_id=>0,
                dbproc=>0,
                ev_type=>0};
    if ($dbproc) {
        $self->{dbproc} = $dbproc;
    }
    if ($asmbl_id) {
        $self->{asmbl_id} = $asmbl_id;
    }
    if ($ev_type) {
        $self->{ev_type} = $ev_type;
    }
    
    bless $self;
    return ($self);
}


sub get_asmbl_id {
    my $self = shift;
    return ($self->{asmbl_id});
}


sub load_predictions {
    print STDERR "Loading Predictions\n";
    my $self = shift;
    my (@genes) = @_;
    my $asmbl_id = $self->{asmbl_id};
    my $dbproc = $self->{dbproc};
    my $ev_type = $self->{ev_type};
    my $x = 0;
    my $assignby = "egc";
    my $prog = $0;
    if ($prog =~ /\/([^\/]+)$/) {
        $prog = $1;
    }
    my $method = $prog;
    
    my $next_model_feat_name = &getNextName ($dbproc, $asmbl_id, "model");
    my $next_exon_feat_name = &getNextName ($dbproc, $asmbl_id, "exon");

    my $feat_name_generator = Gene_feat_name_generator->new($next_model_feat_name, $next_exon_feat_name);

    foreach my $gene (@genes) {

        my $model_feat_name;
        
        eval {
            my $gene_type = $gene->{gene_type};
            unless ($gene_type =~ /protein-coding/) { next;}
            
            $gene->trim_UTRs();
            $gene->refine_gene_object();
            my @exons = $gene->get_exons();
            
            ## BEGIN TRANSACTION (Chunks up the inserts, otherwise Susan gets mad and the transaction log fills up.
            $dbproc->{AutoCommit} = 0; #set autocommit off so transactions can be used.
                        
            ## MODEL RELATED PROCESS
            my $model_span = $gene->{model_span};
            $model_feat_name = ($DEBUG) ? ('NextModel') : $feat_name_generator->next_model_feat_name();
            &Egc_library::insert_asm_feature_only($dbproc, 
                                                  $model_feat_name, 
                                                  $asmbl_id, 
                                                  "model", 
                                                  $model_span->[0], 
                                                  $model_span->[1], 
                                                  $assignby, 
                                                  $method, 
                                                  0, 0);
            
            &Egc_library::insert_phys_ev_only($dbproc, $model_feat_name, $ev_type, $assignby);
            ## reset gene's model_feat_name
            $gene->{Model_feat_name} = $model_feat_name;
            
            ## EXON AND CDS RELATED PROCESS
            foreach my $exon (@exons) {
                my $exon_feat_name = ($DEBUG) ? ('NextExon') : $feat_name_generator->next_exon_feat_name();
                my ($end5, $end3) = ($exon->{end5}, $exon->{end3});
                &Egc_library::insert_asm_feature_only($dbproc, 
                                                      $exon_feat_name, 
                                                      $asmbl_id, 
                                                      "exon", 
                                                      $end5, 
                                                      $end3, 
                                                      $assignby, 
                                                      $method, 
                                                      0, 0);
                
                &Egc_library::insert_phys_ev_only($dbproc, $exon_feat_name, $ev_type, $assignby);
                &insert_feat_link($dbproc, $exon_feat_name, $model_feat_name, $assignby,"now");
                
            }
            
        };

        if (! $@) {
            print "Successful Loading of $model_feat_name\n";
            $dbproc->commit();
        } else {
            print "Sorry, loading of $model_feat_name failed. :(   Rolling back\n";
            print $gene->toString();
            $dbproc->rollback();
            
        }
        
        
    }
    
    ## reset transaction behavior to autocommit:
    $dbproc->{AutoCommit} = 1;
    
}


package Gene_feat_name_generator;
use strict;

my $asmbl_id;
my $next_model_num;
my $next_exon_num;
my $padding_size;

sub new {
    my $packagename = shift;
    my ($next_model, $next_exon) = @_;
    
    #print "next model: $next_model, next exon: $next_exon\n";
    
    $next_model =~ /^(\d+)\.m(\d+)$/;
    $asmbl_id = $1;
    my $model_string = $2;
    $padding_size = length("$model_string");
    print "padding_size: $padding_size\n";
    $next_model_num = int($model_string);
    
    $next_exon =~ /\.e(\d+)$/;
    $next_exon_num = int($1);
    
    my $self = {};
    bless ($self, $packagename);
    return($self);
}


sub next_model_feat_name {
    my $self = shift;
    my $model_feat_name = $self->_create_feat_name("m", $next_model_num);
    $next_model_num++;
    return ($model_feat_name);
}

sub next_exon_feat_name {
    my $self = shift;
    my $exon_feat_name = $self->_create_feat_name("e", $next_exon_num);
    $next_exon_num++;
    return ($exon_feat_name);
}


sub _create_feat_name {
    my $self = shift;
    my ($type, $num) = @_;
    
    #print "creating feat_name: padding_size: $padding_size\n";

    my $padding = $padding_size - length("$num");
    if ($padding < 0) {
        $padding = 0;
    }
    
    my $feat_name = "$asmbl_id.$type" . ("0" x $padding) . "$num";
    return ($feat_name);
}












1; #end of Annot_db_loader













