#!/usr/local/bin/perl

=documentation


my $snsp_analyzer = new SnSp_analysis_manager("genscan", "genemarkHMM", "fgenesh");


my %ev_type__to_gene_list = (    'genscan' => [ $geneA ],
                                 'genemarkHMM' => [ $geneB, $geneC], #predicted two genes here where there's really one
                                 'fgenesh' => [ $geneD ]
                                 );


$snsp_analyzer->add_analysis_entry($gold_standard_gene_obj,  \%ev_type_to_gene_list, 1); # last param is for verbose descriptions



=cut




package main;
our $SEE;

package SnSp_analysis;
use strict;
use warnings;
use Carp;

sub new {
    my $packagename = shift;
    

    confess "This class is abstract, never to be instantiated directly";
    
}

sub _init {
    my $self = shift;
    my @prediction_names = @_;

    unless (@prediction_names) {
        confess "Error, need prediction names for initialization.";
    }

    $self->{AP} = 0;
    $self->{AN} = 0;
    $self->{prediction_to_TFPN} = {};
    $self->{prediction_names} = [@prediction_names];
    $self->{number_genes_analyzed} = 0;
    $self->{number_exons_analyzed} = 0;
}


sub get_prediction_names {
    my $self = shift;
    return (@{$self->{prediction_names}});
}


sub get_TFPN {
    my $self = shift;
    my $ev_type = shift;

    return ($self->{prediction_to_TFPN}->{$ev_type});
}


sub summarize_SnSp {
    my $self = shift;
        
    my ($verbose_summary) = @_;
    
    my $number_genes_analyzed = $self->{number_genes_analyzed};
    my $number_exons_analyzed = $self->{number_exons_analyzed};

    my $AP = $self->{AP};
    my $AN = $self->{AN};

    my @prediction_names = $self->get_prediction_names();
    
    foreach my $ev_type (@prediction_names) {
        my $TFPN_obj = $self->get_TFPN($ev_type);
        
        my $TP = $TFPN_obj->{TP};
        my $TN = $TFPN_obj->{TN};
        my $FP = $TFPN_obj->{FP};
        my $FN = $TFPN_obj->{FN};
        
        my $PP = $TFPN_obj->{PP};
        my $PN = $TFPN_obj->{PN};
        
        my $number_predictions = $TFPN_obj->{number_predicted_genes};
        
        my $number_exons_correct = $TFPN_obj->{number_exons_correct};
        my $number_correct_genes = $TFPN_obj->{genes_predicted_correct};
        
        
        # Sn = TP/(TP+FN)
        my $sensitivity_val = $TP / ($TP + $FN);
        $TFPN_obj->{sensitivity} = $sensitivity_val;
        
        # Sp = TP/(TP+FP)
        my $specificity_val = 0;
        if ($TP || $FP) {
            $specificity_val = $TP / ($TP + $FP);
        }
       
        $TFPN_obj->{specificity} = $specificity_val;
        
        ## Correlation coeff
        my $correl_coeff = ($TP * $TN + $FP * $FN) * ( ($PP * $PN * $AP * $AN) ** (-0.5));
        $TFPN_obj->{correlation_coeff} = $correl_coeff;
        
        ## Fscore
        my $Fscore = 0;

        if ($sensitivity_val > 0 && $specificity_val > 0) {
            $Fscore = (2 * $sensitivity_val * $specificity_val) / ($sensitivity_val + $specificity_val);
            $TFPN_obj->{Fscore} = $Fscore;
        }
        
        if ($verbose_summary) {
            printf("EvType $ev_type ($number_predictions preds)\tsensitivity: %.2f\tspecificity: %.2f\tCC: %.2f\tCexons: %d=%.2f\tCpred: %d=%.2f\tFscore: %.2f\n",  
                   $sensitivity_val, $specificity_val, $correl_coeff, 
                   $number_exons_correct, $number_exons_correct/$number_exons_analyzed*100, 
                   $number_correct_genes, $number_correct_genes/$number_genes_analyzed*100,
                   $Fscore*100);
        }
    } 
}





##########################################################
package SnSp_analysis_manager;
use strict;
use warnings;
use Carp;
use base qw (SnSp_analysis);

sub new {
    my $packagename = shift;
    
    my @prediction_names = @_;
    unless (@prediction_names) {
        confess "wrong args";
    }

    my $self = { 
        intergenic_included => 500, #default
        
    };
    
    
    bless ($self, $packagename);
    

    $self->SUPER::_init(@prediction_names);
    $self->_init();

    return ($self);
}



sub _init {
    my $self = shift;
    
    foreach my $ev_type ($self->get_prediction_names()) {
        $self->{prediction_to_TFPN}->{$ev_type} = TFPN->new();
    }
}


sub get_intergenic_included {
    my $self = shift;
    return ($self->{intergenic_included});
}


sub add_analysis_entry {
    my $self = shift;
    my ($template_gene_obj, $others_hashref, $verbose_flag) = @_;
    
    my $snsp_analysis_entry = SnSp_analysis_entry->new($self, $template_gene_obj, $others_hashref);
    
    $snsp_analysis_entry->calc_SnSp();

    print "\nSingle comparison:\n" if $verbose_flag;
    $snsp_analysis_entry->summarize_SnSp($verbose_flag);
    
    
    ## sum current results to total results:
    $self->_append_analysis_results($snsp_analysis_entry);
    my $number_genes_analyzed = $self->{number_genes_analyzed};
    print "\nTotal so far ($number_genes_analyzed in training set):\n" if $verbose_flag;
    $self->summarize_SnSp($verbose_flag);
    

}


sub _append_analysis_results {
    my ($self, $analysis_ref) = @_;
    

    ## append results of individual genefinders.
    foreach my $pred_name ($self->get_prediction_names()) {
        my $global_TFPN_obj = $self->get_TFPN($pred_name);

        my $analysis_TFPN_obj = $analysis_ref->get_TFPN($pred_name);

        foreach my $att qw (TP FP TN FN PP PN 
                            number_exons_correct 
                            genes_predicted_correct 
                            number_predicted_genes) {
            
            $global_TFPN_obj->{$att} += $analysis_TFPN_obj->{$att};
        }
    }

    ## increment template tallies:
    foreach my $att qw (AP AN number_genes_analyzed number_exons_analyzed) {
        $self->{$att} += $analysis_ref->{$att};
    }
}


################################################################################
package TFPN;
use strict;
use warnings;

sub new {
    my $packagename = shift;
 
    my $self = { TP => 0,  # true positive
                 FP => 0,  # false positive
                 TN => 0,  # true negative
                 FN => 0,  # false negative
                 PP => 0,  # predicted positive
                 PN => 0,  # predicted negative

                 number_exons_correct => 0,
                 genes_predicted_correct => 0,

                 number_predicted_genes => 0,

                 sensitivity => 0,
                 specificity => 0,
                 correlation_coeff => 0,
                 Fscore => 0, #  Fscore = 2SnSp/(Sn+Sp)
                 
             };
    
    bless ($self, $packagename);
    return ($self);
}


##################################################################################

package SnSp_analysis_entry;
use strict;
use warnings;
use base qw (SnSp_analysis);
use Carp;
use Gene_obj_comparator;


sub new {
    my $packagename = shift;
    my $SnSp_analysis_manager = shift;
    unless (ref $SnSp_analysis_manager eq "SnSp_analysis_manager") {
        confess "wrong args.";
    }
    
    my $template_gene_obj = shift;
    
    my $others_hashref = shift;

    my $self = {

        template_gene => $template_gene_obj,
        other_predictions_href => $others_hashref,
        
        manager => $SnSp_analysis_manager, 
    };

    bless ($self, $packagename);


    $self->SUPER::_init($SnSp_analysis_manager->get_prediction_names());
    $self->_init();
    
    return ($self);
}


sub get_gene_predictions {
    my $self = shift;
    my $ev_type = shift;

    my $list_ref = $self->{other_predictions_href}->{$ev_type};
    if (ref $list_ref eq "ARRAY") {
        return (@$list_ref);
    } 
    else {
        ## no predictions here
        confess "Error, no prediction for $ev_type.  If this is true, present an empty array ref.\n";
        return ();
    }
}


sub _init {
    my $self = shift;
    
    foreach my $ev_type ($self->{manager}->get_prediction_names()) {
        $self->{prediction_to_TFPN}->{$ev_type} = TFPN->new();
    }
}

sub calc_SnSp {
    my $self = shift;

    $self->{number_genes_analyzed}++;
    
    my $template_gene_obj = $self->{template_gene};
    my $manager = $self->{manager};

    my $intergenic_included = $manager->get_intergenic_included();
    
    
    ## Process the Template Gene:
    my @gold_std_array;
    
    my ($model_lend, $model_rend) = sort {$a<=>$b} $template_gene_obj->get_model_span();
    
    my ($from_pos, $to_pos) = ($model_lend - $intergenic_included -1, $model_rend + $intergenic_included);
    unless ($from_pos > 0) { 
        # shorter upstream end included since close to seq boundary.
        $from_pos = 1;
    }
        
    ## init gold std array
    for (my $i = 0; $i <= $to_pos-$from_pos; $i++) {
        $gold_std_array[$i] = 0;
    }
    
    ## populate gold_std_array

    my %gold_exons;
    
    my @exons = $template_gene_obj->get_exons();
    foreach my $exon (@exons) {
        my @cds_coords = sort {$a<=>$b} $exon->get_CDS_end5_end3();
        if (@cds_coords) {
            $gold_exons { "$cds_coords[0]" . "_" . "$cds_coords[1]" } = 1; 
            $self->{number_exons_analyzed}++;
            for (my $i = $cds_coords[0]; $i <= $cds_coords[1]; $i++) {
                $gold_std_array[$i-$from_pos] = 1;
            }
        }
    }
    
    ## Score the actual positives and actual negatives:
    for (my $i = 0; $i <= $to_pos-$from_pos; $i++) {
        if ($gold_std_array[$i]) {
            $self->{AP}++;
        } else {
            $self->{AN}++;
        }
    }


    ## Process each of the gene structures:
    
    
    foreach my $ev_type ($self->{manager}->get_prediction_names()) {
        my $same_gene_structure_flag = 0;
        
        my %predicted_exons;
        
        my @prediction_array;
        ## init prediction array
        for (my $i=0; $i <= $to_pos-$from_pos; $i++) {
            $prediction_array[$i] = 0;
        }
        
        my @overlapping_predictions = $self->get_gene_predictions($ev_type);

        foreach my $predicted_gene_obj (@overlapping_predictions) {
            
            compare_genes($template_gene_obj, $predicted_gene_obj);
            if (are_CDS_same()) {
                $same_gene_structure_flag = 1;
            }
            
            my @exons = $predicted_gene_obj->get_exons();
            foreach my $exon (@exons) {
                my @cds_coords = sort {$a<=>$b} $exon->get_CDS_end5_end3();
                if (@cds_coords) {
                    
                    $predicted_exons{ "$cds_coords[0]" . "_" . "$cds_coords[1]" } = 1;
                    
                    for (my $i = $cds_coords[0]; $i <= $cds_coords[1]; $i++) {
                        my $pos = $i - $from_pos;
                        if ($pos >= 0 && $pos <= $to_pos - $from_pos) {
                            $prediction_array[$pos] = 1;
                        }
                    }
                }
            }
        }
		
        my $number_exons_predicted_correct = 0;
        foreach my $exon_key (keys %gold_exons) {
            if ($predicted_exons{$exon_key}) {
                $number_exons_predicted_correct++;
            }
        }

        my $TFPN_obj = $self->get_TFPN($ev_type);

        ## Score the Sensitivity and Specificity parameters:
        
        
        $TFPN_obj->{number_predicted_genes} = scalar (@overlapping_predictions);
        $TFPN_obj->{number_exons_correct} = $number_exons_predicted_correct;
        if ($same_gene_structure_flag) {
            #print "$ev_type, same structure.\n";
            $TFPN_obj->{genes_predicted_correct} = 1;
        } else {
            #print "$ev_type, diff structure.\n";
        }
        
        # score true positives, false positives, and negatives:
        
        for (my $i= 0; $i <= $to_pos-$from_pos; $i++) {
            
            ## TP, FP, TN, FN
            if ($gold_std_array[$i] != 0 && $prediction_array[$i] != 0) {
                $TFPN_obj->{TP}++;
            } elsif ($gold_std_array[$i] != 0 && $prediction_array[$i] == 0) {
                $TFPN_obj->{FN}++;
            } elsif ($gold_std_array[$i] == 0 && $prediction_array[$i] != 0) {
                $TFPN_obj->{FP}++;
            } elsif ($gold_std_array[$i] == 0 && $prediction_array[$i] == 0) {
                $TFPN_obj->{TN}++;
            }
            
            ## PP, PN (predicted positive, predicted negative)
            if ($prediction_array[$i]) {
                $TFPN_obj->{PP}++;
            } else {
                $TFPN_obj->{PN}++;
            }
            
        }
    }
}



1; #EOM
