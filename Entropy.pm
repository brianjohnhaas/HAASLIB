package Entropy;


####
sub compute_entropy {
	my ($string) = @_;
        
	my @chars = split(//, $string);
        
	my %char_counter;
	foreach my $char (@chars) {
		$char_counter{$char}++;
	}
        
	my $entropy = 0;
        
	my $num_chars = length($string);
	foreach my $char (keys %char_counter) {
		my $count = $char_counter{$char};
		my $prob = $count/$num_chars;
                
		my $val = $prob * log(1/$prob)/log(2);
		$entropy += $val;
	}
        
	return($entropy);
}


1; #EOM
