package Feat_name_generator;
use strict;
use warnings;

my $asmbl_id;
my $next_TU_num;
my $next_model_num;
my $next_exon_num;
my $next_cds_num;


my $padding_size;

sub new {
    my $packagename = shift;
    my ($next_TU, $next_model, $next_exon, $next_cds) = @_;
    
    #print "next model: $next_model, next exon: $next_exon\n";
    
    $next_model =~ /^(\d+)\.m(\d+)$/;
    $asmbl_id = $1;
    my $model_string = $2;
    $padding_size = length("$model_string");
    print "padding_size: $padding_size\n";
    $next_model_num = int($model_string);
    
    $next_exon =~ /\.e(\d+)$/;
    $next_exon_num = int($1);
    
    $next_TU =~ /\.t(\d+)$/;
    $next_TU_num = int ($1);

    $next_cds =~ /\.c(\d+)$/;
    $next_cds_num = int($1);

    unless ($next_TU_num && $next_model_num && $next_exon_num && $next_cds_num) {
        die "Error extracting numeric fields: " . 
            "TU: $next_TU_num from $next_TU\n"
            . "model: $next_model_num from $next_model\n"
            . "exon: $next_exon_num from $next_exon\n"
            . "cds: $next_cds_num from $next_cds\n";
    }
    

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


sub next_TU_feat_name {
    my $self = shift;
    my $TU_feat_name = $self->_create_feat_name("t", $next_TU_num);
    $next_TU_num++;
    return ($TU_feat_name);
}


sub next_cds_feat_name {
    my $self = shift;
    my $cds_feat_name = $self->_create_feat_name("c", $next_cds_num);
    $next_cds_num++;
    return ($cds_feat_name);
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
    #print STDERR "created feat_name: $feat_name\n";
    return ($feat_name);
}




1; #EOM













