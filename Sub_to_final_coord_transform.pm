package main;

our $SEE; # main global.

=head1 NAME

Sub_to_final_coord_transform

=cut

=head1 DESCRIPTION

Given a $dbproc to an Eukaryotic annotation database and given population of the sub_to_final table for a given pseudomolecule construct, this module provides facilities to transform coordinates from the pseudochromosome to the assembly component it was constructed from and vice-versa.

=cut


package Sub_to_final_coord_transform;

use strict;
use DBI;
use Egc_library;

my %bac_to_chromo;
my %chromo_to_bac;





=over 4

=item new()

B<Description:> Instantiates a Sub_to_final_coord_transform object. 

B<Parameters:> $dbproc

$dbproc is a DBI connection handle to an Eukaryotic annotation database

B<Returns:> thisObj


The sub_to_final table is queried and scaffolds are stored in memory.

=back

=cut

sub new {
    my $classname = shift;
    my $dbproc = shift;
    unless (ref $dbproc) {
	die "Sub_to_final_coord_transform: ERROR: I need a db-connection ref.\n";
    }
    
    ## query sub_to_final;
    my $query = "select s.asmbl_id, s.asm_lend, s.asm_rend, s.sub_asmbl_id, s.sub_asm_lend, s.sub_asm_rend from sub_to_final s, clone_info c where s.asmbl_id = c.asmbl_id and c.final_asmbl = c.asmbl_id";
    my @results = &do_sql_2D($dbproc, $query);
    foreach my $result (@results) {
	my ($asmbl_id, $asm_lend, $asm_rend, $sub_asmbl_id, $sub_asm_lend, $sub_asm_rend) = @$result;
	
	$chromo_to_bac{$asmbl_id}->{$sub_asmbl_id} = {
	    asmbl_id=>$asmbl_id,
	    asm_lend=>$asm_lend,
	    asm_rend=>$asm_rend,
	    sub_asmbl_id=>$sub_asmbl_id,
	    sub_asm_lend=>$sub_asm_lend,
	    sub_asm_rend=>$sub_asm_rend
	    };
	
	$bac_to_chromo{$sub_asmbl_id}->{$asmbl_id} = { 
	    asmbl_id=>$asmbl_id,
	    asm_lend=>$asm_lend,
	    asm_rend=>$asm_rend,
	    sub_asmbl_id=>$sub_asmbl_id,
	    sub_asm_lend=>$sub_asm_lend,
	    sub_asm_rend=>$sub_asm_rend
	    };
    }

    my $self = {};
    bless $self, $classname;
    return ($self);
}





=over 4

=item convert_bac_to_chromo()

B<Description:> Object method converts bac-based coordinates to a chromosome coordinates.

B<Parameters:> bac_asmbl_id, chromo_asmbl_id, coord1, coord2, ...

The first parameter is the asmbl_id corresponding to the tiling path component
The following parameters represent a list of coordinates to be converted.  At least one coordinate should
be provided.


B<Returns:> @ret_coords

@ret_coords is a list of the converted coordinates in the same order as the inputted list of bac coordinates.

=back

=cut


sub convert_bac_to_chromo {
    my $self = shift;
    my $bac_asmbl = shift;
    my $chromo_asmbl = shift;
    my @coords = @_;
    
    my $info_ref = $bac_to_chromo{$bac_asmbl}->{$chromo_asmbl};
    unless (ref $info_ref) {
	die "Can't find information linking bac:$bac_asmbl to chromo:$chromo_asmbl\nFatal.\n";
    }
    my ($bac_lend, $bac_rend, $chromo_lend) = ($info_ref->{sub_asm_lend}, 
					       $info_ref->{sub_asm_rend},
					       $info_ref->{asm_lend});
    

    my $converter_routine = &get_bac_to_chromo_converter_routine($bac_lend, $bac_rend, $chromo_lend);
    
    my @ret_coords;
    foreach my $coord (@coords) {
	my $ret_coord = &$converter_routine($coord);
	push (@ret_coords, $ret_coord);
    }
    
    return (@ret_coords);

}




=over 4

=item convert_chromo_to_bac()

B<Description:> object method converts chromosome coordinates to bac coordinates.

B<Parameters:> chromo_asmbl_id, coord1, coord2, ...

Chromosome asmbl_id is required as the first parameter, followed by the list of coordinates to be converted.

B<Returns:> @converted_coord_structs

@converted_coord_structs contains a list of struct elements (hashrefs), each struct referring to the converted coordinate data in the same order as the inputted list of coordinates.
The struct has the following structure:

    $struct = {  
	bac_asmbl_id => $bac_asmbl_id,
	coord => $converted_coordinate 
	}

The coordinates are converted only to the BAC and sequence range built into the pseudochromosome sequence and specified by the sub_to_final table.


=back

=cut


sub convert_chromo_to_bac {
    my $self = shift;
    my $chromo_asmbl = shift;
    my @coords = @_;

    my @ret_coords;

    foreach my $coord (@coords) {
	
	my $data_ref = &map_to_bac($chromo_asmbl, $coord);
	
	unless (ref $data_ref) {
	    die "Error: Sub_to_final_coord_transform: can't find coord ($coord) on chromosome.\n";
	}
	my ($bac_lend, $bac_rend, $chromo_lend, $chromo_rend, $bac_asmbl) = ($data_ref->{sub_asm_lend},
									     $data_ref->{sub_asm_rend},
									     $data_ref->{asm_lend},
									     $data_ref->{asm_rend},
									     $data_ref->{sub_asmbl_id});
	
	my $converter_routine = &get_chromo_to_bac_converter_routine($bac_lend, $bac_rend, $chromo_lend);
	my $ret_coord = &$converter_routine($coord);
	push (@ret_coords, {bac_asmbl_id => $bac_asmbl, coord => $ret_coord});
	print "chromo($chromo_asmbl) coord($coord) --> bac($bac_asmbl) coord($ret_coord)\n" if $SEE;
	
    }
    return (@ret_coords);
}



###################
## Helper subs
###################

sub map_to_bac {
    my $chromo_asmbl = shift;
    my $coord = shift;
    foreach my $asmbl (keys %{$chromo_to_bac{$chromo_asmbl}}) {
	my $data_ref = $chromo_to_bac{$chromo_asmbl}->{$asmbl};
	my ($bac_lend, $bac_rend, $chromo_lend, $chromo_rend) = ($data_ref->{sub_asm_lend},
								 $data_ref->{sub_asm_rend},
								 $data_ref->{asm_lend},
								 $data_ref->{asm_rend});
	
	print "map_to_bac: examining bac($asmbl) chr($chromo_lend, $chromo_rend) coordinate $coord.\n" if $SEE;
	if ($coord >= $chromo_lend && $coord <= $chromo_rend) {
	    return ($data_ref);
	}
    }
}




sub get_bac_to_chromo_converter_routine {
    my ($bac_lend, $bac_rend, $chromo_lend) = @_;
    
    my $adj_sub;
    if ($bac_lend <= $bac_rend) {
	
	$adj_sub = sub {
	    my ($coord) = @_;
	    my $diff =  $coord - $bac_lend;
	    my $new_coord = $diff + $chromo_lend;
	    print "$coord --> $new_coord, method: 1\n" if $SEE;
	    return ($new_coord);
	};
    } elsif ($bac_lend > $bac_rend) {
	#reverse strand.  Must revcomp the coordinates.
	$adj_sub = sub {
	    my ($coord) = @_;
	    my $diff = $coord - $bac_lend;
	    my $new_coord = $chromo_lend - $diff; 
	    print "$coord --> $new_coord, method: 2\n" if $SEE;
	    return ($new_coord);
	}
    }
    
    return ($adj_sub);
}



sub get_chromo_to_bac_converter_routine {
    my ($bac_lend, $bac_rend, $chromo_lend) = @_;
    my $adj_sub;
    if ($bac_lend <= $bac_rend) {
	
	$adj_sub = sub {
	    	    
	    ## reverse: given new_coord, calculate $coord
	    my $new_coord = shift;
	    
	    my $coord = $new_coord - $chromo_lend + $bac_lend;
	    
	    print "$new_coord --> $coord, method 3.\n" if $SEE;

	    return ($coord);

	};
    } elsif ($bac_lend > $bac_rend) {
	#reverse strand.  Must revcomp the coordinates.
	$adj_sub = sub {
	    my $new_coord = shift;

	    my $coord = $chromo_lend - $new_coord + $bac_lend;
	   
	    print "$new_coord --> $coord, method 4.\n" if $SEE;
	    
	    return ($coord);
	    
	}
    }
    
    return ($adj_sub);
}





1;
