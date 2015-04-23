#!/usr/local/bin/perl

use strict;
use warnings;
use lib ($ENV{EUK_MODULES}, $ENV{EGC_SCRIPTS});
use Gene_obj_indexer;
use Egc_library;
use CdbTools;

my $usage = "usage: $0 gene_obj_inx genome_fasta_db [prot(default)|CDS|cDNA]\n\n";
my $gene_obj_inx_file = $ARGV[0] or die $usage;
my $genome_fasta_db = $ARGV[1] or die $usage;
my $type = $ARGV[2] || "prot";

unless ($type =~ /prot|CDS|cDNA/) {
    die "Error, don't understand type: $type\n";
}

my $gene_obj_indexer = new Gene_obj_indexer( { "use" => $gene_obj_inx_file } );

my @gene_ids = $gene_obj_indexer->get_keys();

my %asmbl_id_to_gene_list;

foreach my $gene_id (@gene_ids) {
    #print "$gene_id\n";
    my $gene = $gene_obj_indexer->get_gene($gene_id);
    my $asmbl_id = $gene->{asmbl_id};

    #print "\t$asmbl_id\n";

    my $gene_list_aref = $asmbl_id_to_gene_list{$asmbl_id};
    unless (ref $gene_list_aref) {
        $gene_list_aref = $asmbl_id_to_gene_list{$asmbl_id} = [];
    }
    
    push (@$gene_list_aref, $gene_id);
    
}

foreach my $asmbl_id (keys %asmbl_id_to_gene_list) {
    
    my $genome_seq = cdbyank_linear($asmbl_id, $genome_fasta_db);

    my @gene_ids = @{$asmbl_id_to_gene_list{$asmbl_id}};

    foreach my $gene_id (@gene_ids) {
        
        my $gene_obj_ref = $gene_obj_indexer->get_gene($gene_id);
        $gene_obj_ref->create_all_sequence_types(\$genome_seq);
        
        foreach my $gene_obj ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {

            my $model_id = $gene_obj->{Model_feat_name};
            
            my $seq = "";
            if ($type eq "prot") {
                $seq = $gene_obj->get_protein_seq();
            }
            elsif ($type eq "CDS") {
                $seq = $gene_obj->get_CDS_sequence();
            }
            elsif ($type eq "cDNA") {
                $seq = $gene_obj->get_cDNA_sequence();
            }
            
            $seq = &make_FASTA_format($seq);
            chomp $seq;

            print ">$model_id\n$seq\n";
        }
    }
}



exit(0);


