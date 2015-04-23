## subroutines in common to TIGR_GENBANK_comparator.dbi and TIGR_MIPS_com







####
sub get_overlap_length {
    my ($c1, $c2, $x1, $x2) = @_;
    my $xlength = $x2 - $x1 + 1;
    my $delta1 = ( ($x2 - $c2) < 0) ? 0 : ($x2 - $c2);
    my $delta2 = ( ($c1 - $x1) < 0) ? 0 : ($c1 - $x1);
    return ($xlength - $delta1 - $delta2);
}

####
sub decent_overlap {
    my ($length, $overlap) = @_;
    ## avoid using corrupt data
    if ($length > 50000) {
	return (0);
    }
    if ( ( ($overlap/$length) * 100) >= $OVERLAP_REQUIRED) {
	print "\tacceptable\n";
	return (1);
    } else {
	print "\tnot acceptable\n";
	return (0);
    }
}


1; #end of library file.
