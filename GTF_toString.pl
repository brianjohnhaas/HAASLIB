#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ($FindBin::Bin);
use Gene_obj;
use GTF_utils;

my $usage = "usage: $0 annots.gtf\n\n";

my $gtf_file = $ARGV[0] or die $usage;

main: {
    
    
    my $gene_obj_indexer_href = {};

    my $contig_to_gene_list_href;

    print STDERR "-parsing gene annotations file\n";
    # GTF mode
    $contig_to_gene_list_href = &GTF_utils::index_GTF_gene_objs_from_GTF($gtf_file, $gene_obj_indexer_href);
    
    foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
        
        my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
        
        foreach my $gene_id (@gene_ids) {
            my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};
            
            print $gene_obj_ref->toString() . "\n";
        }
    }
    

    exit(0);
    
}
