

package Loci_assigner;
use strict;

my $DEBUG;

## Globals
my @genes;

my %loci;
my $core_locus;
my $AVG_DIST;
my $PAD_TOTAL = 0;

####
sub assign_loci {
    my ($gene_struct_list_aref, $preexisting_loci_href) = @_;
    ## gene struct should have the following format:
    #   gene = { locus => (value || undef),
    #            feat_name => $feat_name,
    #            end5 => $end5,                         
    #            end3 => $end3 }
    #
    ## Attributes added herein:
    #    loc_pos    #integer value of trailing locus ID
    #    core_locus #prefix to locus assignment minus ".loc_pos" suffix
    #    newLocus   #value of newly assigned locus
    #    index      #int specifying index in the sorted gene list
    

    ## init globals
    $core_locus = undef;
    my @gene_structs = @$gene_struct_list_aref;
    %loci = %$preexisting_loci_href;
    $AVG_DIST = 10; #default;
    
    ## sort structs by midPt
    @gene_structs = sort {  ( ($a->{end5} + $a->{end3})/2)
				<=>
				( ($b->{end5} + $b->{end3}) / 2) } @gene_structs;
    
    
    ## assign indices:
    for (my $i=0; $i <= $#gene_structs; $i++) {
	$gene_structs[$i]->{index} = $i;
    }
    
    foreach my $gene_struct (@gene_structs) {
	my $locus = $gene_struct->{locus};
	if ($locus && $locus =~ /^(.*)\.(\d+)$/) {
	    $loci{$locus} = 1; #store that locus already exists.
	    my $core_loc_val = $1;
	    my $loc_pos = $2;
	    if ($loc_pos =~ /^0/) { # padded number
		$PAD_TOTAL = length ($loc_pos);
	    }
	    $gene_struct->{loc_pos} = $loc_pos; ## add attribute here.
	    $gene_struct->{core_locus} = $core_loc_val;
	} elsif ($locus) {
	    die "Error, locus ($locus) lacks the expected format.\n";
	}
    }
    
    ## determine avg distance between assigned locus positions.:
    my $count;
    my $sum;
    for (my $i=0; $i < $#gene_structs; $i++) {
	my $curr_pos = $gene_structs[$i]->{loc_pos};
	my $next_pos = $gene_structs[$i+1]->{loc_pos};
	
	my $curr_core = $gene_structs[$i]->{core_locus};
	my $next_core = $gene_structs[$i]->{core_locus};
	
	unless ($curr_core eq $next_core) {
	    next;
	}
	
	if (defined ($curr_pos) && defined($next_pos)) {
	    $sum += abs ($next_pos - $curr_pos);
	    $count++;
	}
    }
    if ($count) {
	$AVG_DIST = int ($sum/$count + 0.5);
    }
        
    ## Perform Loci suggestions.
    my $startVal = 0;
    my @loci_to_assign;
    for (my $i=0; $i <= $#gene_structs; $i++) {
	my $struct = $gene_structs[$i];
	my $locus = $struct->{locus};
	if (!defined($locus)) {
	    push (@loci_to_assign, $struct);
	} else {
	    my $endVal = $struct->{loc_pos};
	    if (@loci_to_assign) {
		&process_new_loci(\@gene_structs, \@loci_to_assign, $startVal, $endVal);
	    }
	    @loci_to_assign = (); #reinitialize
	    $startVal = $endVal;
	}
    }
    
    if (@loci_to_assign) {
	&process_new_loci (\@gene_structs, \@loci_to_assign, $startVal, undef);
    }
    
    
    return (@gene_structs);
    
    exit(0);
}


####
sub process_new_loci {
    my ($gene_structs_aref, $list_to_assign_aref, $start,$end) = @_;
    my @list = @$list_to_assign_aref;
    my $num_elements = $#list + 1;
    print "Number of elements: $num_elements\n" if $DEBUG;
    my $increment = 0;
    print "start: $start, end: $end\n" if $DEBUG;
    
    ## check for consistent core_locus (adjacent genes should be list boundaries or genes with assigned loci)
    my $core_locus;
    my $consistent_core_locus = 1;
    my $left_assigned_locus;
    my $right_assigned_locus;
    my $first_entry_index = $list[0]->{index};
    my $last_entry_index = $list[$#list]->{index};
    if ($first_entry_index > 0) {
	$left_assigned_locus = $gene_structs_aref->[$first_entry_index-1]->{core_locus};
    }
    if (my $last_struct = $gene_structs_aref->[$last_entry_index+1]) {
	$right_assigned_locus = $last_struct->{core_locus};
    }
    unless ($left_assigned_locus || $right_assigned_locus) {
	die "Error, no locus adjacent to list of genes requiring assignments\n";
    }
    if ($left_assigned_locus && $right_assigned_locus && (lc($left_assigned_locus) ne lc($right_assigned_locus))) {
	$consistent_core_locus = 0;
    }
    if ($left_assigned_locus) {
	$core_locus = $left_assigned_locus;
    } else {
	$core_locus = $right_assigned_locus;
    }
        
    if ($consistent_core_locus && $end) {
	## determine incrementation:
	my $diff = $end - $start;
	print "diff\t$diff\n" if $DEBUG;
        $increment = int($diff/($num_elements+1) + 0.5);
	print "increment: $increment\n" if $DEBUG;
	if ($increment < 1) {
	    $increment = 1;
	}
    } else {
	$increment = $AVG_DIST; #default
    }
    print "increment: $increment\n" if $DEBUG;
    my $begin = $start;
    foreach my $struct (@list) {
	$begin += $increment;
	
	if ($PAD_TOTAL && (length("$begin") < $PAD_TOTAL)) {
	    my $num_pad = $PAD_TOTAL - length("$begin");
	    $begin = ('0' x $num_pad) . $begin;
	}
	
	my $tentative_locus = $core_locus . "." . $begin;
	print "tentative_locus: $tentative_locus\n" if $DEBUG;
	while ($loci{$tentative_locus}) {
	    $begin += 1;
	    if ($PAD_TOTAL && (length("$begin") < $PAD_TOTAL)) {
		my $num_pad = $PAD_TOTAL - length("$begin");
		$begin = ('0' x $num_pad) . $begin;
	    }
	    $tentative_locus = $core_locus . "." . $begin;
	    print "\ttentative locus in use. Trying new one: $tentative_locus\n" if $DEBUG;
	}
	$struct->{newLocus} = $tentative_locus;
	$loci{$tentative_locus} = 1;
    }
}

1;

				
    
    

		   

