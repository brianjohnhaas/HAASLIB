package main;

our ($IMAGE_X_SIZE, $DEBUG);

package Genes_to_image;
use strict;
use Sequence_coords_image;
use CDNA::CDNA_alignment;
use Gene_obj;


my $consensus = 10;
my $consensus_color = "blue";

sub create_image {
    my (@genes) = @_;
    my @data;
    foreach my $obj (@genes) {
	my $struct = &get_struct_from_gene_obj($obj);
    	push (@data, $struct);
    }
    
    my %info;
    
    my @c;
    my $highlight = 0;
    foreach my $gene_struct (@data) {
	my $coords_aref = $gene_struct->{coords};
	foreach my $href (@$coords_aref) {
	    my $c1 = $href->{end5};
	    my $c2 = $href->{end3};
	    if ($gene_struct->{type} eq "gene" && (!$highlight)) {
		$info{$c1} = 1;
		$info{$c2} = 1;
	    }
	    push (@c, $c1, $c2);
	}
	$highlight = 1; #only track coords for first gene.
    }
    @c = sort {$a<=>$b} @c;
    my $min_end5 = shift @c;
    $min_end5 -= 1;
    
    if ($DEBUG) {
	print "---\nMIN_COORD: $min_end5\n----\n";
	print "CREATING IMAGE\n";
    }
    
    
    my $obj = new Sequence_coords_image(
					IMAGE_X_SIZE => $IMAGE_X_SIZE || 750,  #default image length in pixels
					DRAW_PANEL_SCALER => 0.6, #percentage of image to draw matches, rest for text-annotation of match
					ELEMENT_VERTICAL_SPACING => 15, #amount of vertical space consumed by each element
					TICKER_TOGGLE => 1 #toggle for displaying a ticker for protein length, default is on.
					#SEQ_START => $min_end5
					);
    
    foreach my $alignment (@data) {
	my $accession = $alignment->{name};
	my @coordsets = @{$alignment->{coords}};
	
	my @list; #initialize 
	foreach my $coordset (@coordsets) {
	    my $end5 = $coordset->{end5};
	    my $end3 = $coordset->{end3};
	   	    
	    my $end5_agree = $info{$end5};
	    my $end3_agree = $info{$end3};
	    $end5 -= $min_end5;
	    $end3 -= $min_end5;
	    if ($DEBUG) {
		print "END5: $end5\tEND3: $end3\t$accession\n";
	    }
	    my $color = ($accession =~ /\sFL/) ? "119:119:119" : "black";
	    push (@list, $end5, "full", $end3, "full", $color);
	    if ((my $cds_end5 = $coordset->{cds_end5}) && (my $cds_end3 = $coordset->{cds_end3})) {
		$cds_end5 -= $min_end5;
		$cds_end3 -= $min_end5;
		push (@list, $cds_end5, "full", $cds_end3, "full", "red");
	    }
	    if ($end5_agree) {
		push (@list, $end5, "full", ($end5 + $consensus) , "full", $consensus_color);
	    }
	    if ($end3_agree) {
		push (@list, ($end3 - $consensus), "full", $end3, "full", $consensus_color);
	    }
	    
	    # element params (end5, end5_type, end3, end3_type, color, [... repeat for each match ....])
	    #  end_types: (arrow, full, partial)
	    # colors ( (white, blue, red, black, or green) or ("R:G:B") )
	}
	my $element = new Sequence_coords_image::Element(@list);
	$element->set_text($accession);
	$obj->add_element($element);
    }
    
    #print image
    my $img = $obj->create_image();
    return ($img);
}
 

####
sub get_struct_from_gene_obj {
    my $gene_obj = shift;
    my @coords;
    my $strand = $gene_obj->{strand};
    
    my @exons = $gene_obj->get_exons();
    my %cds;
    my $i = 0;
    foreach my $exon (@exons) {
	my ($end5, $end3) = sort {$a<=>$b} $exon->get_coords();
	$coords[$i]->{end5} = $end5;
	$coords[$i]->{end3} = $end3;
	if (my $cds_ref = $exon->get_CDS_obj()) {
	    my ($end5, $end3) = sort {$a<=>$b} $cds_ref->get_coords();
	    $coords[$i]->{cds_end5} = $end5;
	    $coords[$i]->{cds_end3} = $end3;
	}
	$i++;
    }
    my $com_name = $gene_obj->{com_name} || "";
    my $struct = { type => 'gene',
		   strand => $strand,
		   coords => \@coords,
		   name => "($strand)" . $gene_obj->{Model_feat_name} . " $com_name"};
    return ($struct);
    
}

1; #EOM
