#!/usr/local/bin/perl

use strict;
use warnings;
use lib ($ENV{EUK_MODULES}, $ENV{EGC_SCRIPTS});

use Mummer_delta_parser;
use Mummer_remap;


my $delta_file = "delta";

my $mummer_remap = new Mummer_remap($delta_file);
my $mummer_delta_parser = $mummer_remap->get_delta_parser(); 

## print out descriptions of all the alignments parsed from the delta file:
my @alignments = $mummer_delta_parser->get_alignments();
foreach my $alignment (@alignments) {
    print $alignment->toString() . "\n\n";
    
    ## some useful methods
    my $reference_accession = $alignment->get_reference_accession();
    my $query_accession = $alignment->get_query_accession();

    my $reference_alignment_length = $alignment->get_reference_coord_span_length();
    my $query_alignment_length = $alignment->get_query_coord_span_length();

    ## remember, the above alignments may contain indels. 
    ## if you want the ungapped portions of each alignments, get them like so:
    my @msps = $alignment->get_MSPs();
    foreach my $msp (@msps) {
        # the indel-containing alignment inherits from the msp class, 
        # so you can use the same methods here for getting lengths
        # and coordinates
        
        # ie. print $msp->toString(); (not done again here, since already done in above
        #    $alignment->toString() call.
    }


}

# based on the delta parsing, it appears:
# Ref: gec1-200, Query: ntec02-2




## do coordinate conversions:
open (my $fh, "last_remap") or die $!;
<$fh>; #remove header

while (<$fh>) {
    chomp;
    
    if (/^\#/) { 
        next; 
    }
    
    my @x = split (/\t/);
    
    my ($ref_db, $ref_asmbl_id, $ref_feat_type, $ref_feat_name, $ref_end5, $ref_end3, 
        $targ_db, $targ_asmbl_id,  $targ_feat_name, $targ_end5, $targ_end3,
        $problem) = split (/\t/, $_, 12);
    
    unless ($targ_db && $targ_asmbl_id) {
        ## nothing to match it to...
        next;
    }
    
    my ($ntec_end5, $ntec_end3) = ($targ_end5, $targ_end3);
    my ($gec_end5, $gec_end3) = ($ref_end5, $ref_end3);
    
    ## must specify reference and query accessions here.  In multi-fasta vs. multi-fasta, these are not constant
    
    my $accs_href = { reference_accession =>  "$ref_db-$ref_asmbl_id",  # accs used in the nucmer search
                      query_accession => "$targ_db-$targ_asmbl_id"
                      };
    

    ## going gec to ntec coordinate conversion

    my ($adj_gec_end5, $adj_gec_end3) = $mummer_remap->transform_coordinates($Mummer_remap::REFERENCE_TO_QUERY,
                                                                             $accs_href,
                                                                             ($gec_end5, $gec_end3));
    
    print "transformed: ($gec_end5, $gec_end3): ($adj_gec_end5, $adj_gec_end3)\n";
    
}




exit(0);

