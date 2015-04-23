#!/usr/local/bin/perl


use lib ($ENV{EUK_MODULES}, $ENV{EGC_SCRIPTS});
use Egc_library;
use strict;
use DBI;
use Gene_obj;
use Data::Dumper;
use Gene_obj_indexer;


my $usage = "usage: $0 gff3_file gene_obj_output_inx_file\n\n";
my $gff3_file = $ARGV[0] or die $usage;
my $gene_obj_inx_file = $ARGV[1] or die $usage;


my $gene_obj_indexer = new Gene_obj_indexer();
$gene_obj_indexer->make_index_file($gene_obj_inx_file);

my %transcript_to_gene;
my %gene_coords;


my $counter = 0;
open (my $fh, $gff3_file) or die $!;
while (<$fh>) {
    chomp;
    my @x = split (/\t/);
    my ($asmbl_id, $feat_type, $lend, $rend, $orient, $gene_info) = ($x[0], $x[2], $x[3], $x[4], $x[6], $x[8]);    
        
    unless ($feat_type =~ /^(mRNA|CDS|exon)$/) { next;}

    $gene_info =~ /ID=(\S+);/;
    my $id = $1 or die "Error, couldn't get the id field $_";
    
    $gene_info =~ /Parent=(\S+);?/;
    my $parent = $1 or die "Error, couldn't get the parent info $_";
    
    if ($feat_type eq 'mRNA') {
        ## just get the identifier info
        $transcript_to_gene{$id} = $parent;
        next;
    }

    my $transcript_id = $parent;
    my $gene_id = $transcript_to_gene{$transcript_id} or die "Error, no gene_id for $parent";
    

    my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);

    $gene_coords{$asmbl_id}->{$gene_id}->{$transcript_id}->{$feat_type}->{$end5} = $end3;
    # print "$asmbl_id, $gene_id, $transcript_id, $feat_type, $end5, $end3\n";
    

    $counter++;
    #if ($counter > 300) { last;}

}
close $fh;


## 
foreach my $asmbl_id (sort {$a<=>$b} keys %gene_coords) {
    
    #my $genome_seq = &get_seq($dbproc, $asmbl_id);
    
    my $genes_href = $gene_coords{$asmbl_id};
    foreach my $gene_id (keys %$genes_href) {
        
        my $transcripts_href = $genes_href->{$gene_id};
        
        my @gene_objs;
        
        foreach my $transcript_id (keys %$transcripts_href) {
            
            my $cds_coords_href = $transcripts_href->{$transcript_id}->{CDS};
            my $exon_coords_href = $transcripts_href->{$transcript_id}->{exon};
            
            unless (ref $cds_coords_href && ref $exon_coords_href) {
                print STDERR Dumper ($transcripts_href);
                die "Error, missing cds or exon coords for $transcript_id, $gene_id\n";
            }
            

            my $gene_obj = new Gene_obj();

            $gene_obj->populate_gene_obj($cds_coords_href, $exon_coords_href);
            
            $gene_obj->{Model_feat_name} = $transcript_id;
            $gene_obj->{TU_feat_name} = $gene_id;
            $gene_obj->{asmbl_id} = $asmbl_id;
            
            push (@gene_objs, $gene_obj);
            
            
        }

        my $template_gene_obj = shift @gene_objs;
        foreach my $gene_obj (@gene_objs) {
            $template_gene_obj->add_isoform($gene_obj);
        }
        
        $gene_obj_indexer->store_gene($gene_id, $template_gene_obj);
        print "stored $gene_id\n";
        
    }
}

exit(0);


    
