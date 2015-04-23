##########################################################################
# Written by Paolo Amedeo, slightly modified by bhaas (1/7/2005) #########
##########################################################################


=head1

Rainbow.pm

=cut


package Rainbow;

use strict;


sub new {
    my $packagename = shift;
    my $self = {};
    bless ($self, $packagename);
    return ($self);
}



sub GetRainbow_K_R { # $no_rainbow_cols = GetRainbow_B_R(\@rainbow_colors, $no_rainbow_cols)
	shift();
	my ($r_rainbow, $col_no) = @_;
	die "\n\nGetRainbow: It is possible to handle only a positive number of colors\n\n" unless $col_no > 0;
	my ($red, $green, $blue) = (0) x 3;
	my @step = (0) x 5;	# 5 are the steps for generating the different colors:
				#   0,   0,   0  ->   0,   0, 255 (black  -> blue)
				#   0,   0, 255  ->   0, 255, 255 (blue   -> cyan)
				#   0, 255, 255  ->   0, 255,   0 (cyan   -> green)
				#   0, 255,   0  -> 255, 255,   0 (green  -> yellow)
				# 255, 255,   0  -> 255,   0,   0 (yellow -> red)

	# "Rude way of distributing the colors...."
	my $n = @step; # 1 + last index in the array @step_no.  We want to have more gradients on the top than on the bottom => stepping faster out of the black
	
	while (--$col_no){ # until we don't run out of colors
		++$step[--$n];
		$n = @step unless $n; # new round of distribution	
	}

	# Converting the number of steps in steps
	
	foreach my $incr (@step){
		$incr = 255 / $incr;
	}
	push(@{$r_rainbow}, [$red, $green, $blue]); # black
	
	# black -> blue
	
	while($blue < 255){
		$blue += $step[0];
		$blue = 255 if 255 - $blue < $step[0];
		push(@{$r_rainbow}, [$red, $green, int($blue)]);
	}
	# blue -> cyan
	
	while ($green < 255){
		$green += $step[1];
		$green = 255 if 255 - $green < $step[1];
		push(@{$r_rainbow}, [$red, int($green), $blue]);
	}
	# cyan -> green
	
	while ($blue > 0){
		$blue -= $step[2];
		$blue = 0 if $blue < $step[2];
		push(@{$r_rainbow}, [$red, $green, int($blue)]);
	}
	# green -> yellow
	
	while ($red < 255){
		$red += $step[3];
		$red = 255 if 255 - $red < $step[3];
		push(@{$r_rainbow}, [int($red), $green, $blue]);
	}
	# yellow -> red
	
	while ($green > 0){
		$green -= $step[4];
		$green = 0 if $green < $step[4];
		push(@{$r_rainbow}, [$red, int($green), $blue]);
	}
	return(scalar(@{$r_rainbow}));
}

sub GetRainbow_B_R { # $no_rainbow_cols = GetRainbow_B_R(\@rainbow_colors, $no_rainbow_cols)
	shift();
	my ($r_rainbow, $col_no) = @_;
	die "\n\nGetRainbow: It is possible to handle only a positive number of colors\n\n" unless $col_no > 0;
	my ($red, $green, $blue) = (0, 0, 255);
	my @step = (0) x 4;	# 4 are the steps for generating the different colors:
				#   0,   0, 255  ->   0, 255, 255 (blue   -> cyan)
				#   0, 255, 255  ->   0, 255,   0 (cyan   -> green)
				#   0, 255,   0  -> 255, 255,   0 (green  -> yellow)
				# 255, 255,   0  -> 255,   0,   0 (yellow -> red)

	# "Rude way of distributing the colors...."
	my $n = @step; # 1 + last index in the array @step_no.  We want to have more gradients on the top than on the bottom => stepping faster out of the black
	
	while (--$col_no){ # until we don't run out of colors
		++$step[--$n];
		$n = @step unless $n; # new round of distribution
	}

	# Converting the number of steps in steps
	
	foreach my $incr (@step){
		$incr = 255 / $incr;
	}
	push(@{$r_rainbow}, [$red, $green, $blue]); # blue
	
	# blue -> cyan
	
	while ($green < 255){
		$green += $step[0];
		$green = 255 if 255 - $green < $step[0];
		push(@{$r_rainbow}, [$red, int($green), $blue]);
	}
	# cyan -> green
	
	while ($blue > 0){
		$blue -= $step[1];
		$blue = 0 if $blue < $step[1];
		push(@{$r_rainbow}, [$red, $green, int($blue)]);
	}
	# green -> yellow
	
	while ($red < 255){
		$red += $step[2];
		$red = 255 if 255 - $red < $step[2];
		push(@{$r_rainbow}, [int($red), $green, $blue]);
	}
	# yellow -> red
	
	while ($green > 0){
		$green -= $step[3];
		$green = 0 if $green < $step[3];
		push(@{$r_rainbow}, [$red, int($green), $blue]);
	}
	return(scalar(@{$r_rainbow}));
}


1; #EOM

__END__


=over 4

=item SampleProgram

#!/usr/local/bin/perl

use strict;
use Rainbow;


my $ptt = Rainbow->new();

my $col_no = 128;

my @col_list = ();

# Number_of_allocated colors = $ptt->GetRainbow_B_R(\@col_list, Number_of_requested_Colors)

$ptt->GetRainbow_B_R(\@col_list, $col_no) || die "\n\nproblems assigning the colors\n\n";

print "<html>\r\n<head>\r\n<title>Color Test</title>\r\n</head>\r\n",
	"<body>\r\n<table cols=128 valign=\"center\" align=\"center\">\r\n\t<tr>\r\n";

foreach my $color (@col_list){
	my $hex_col = sprintf("%02x%02x%02x", @{$color}); # Converting in a string of hexadecimal pairs

	print "\t\t<td bgcolor=\"#$hex_col\">$hex_col</td>\r\n";
}
print "\t</tr>\r\n</table>\r\n</body>\r\n</html>\r\n";

=back

=cut

    ;


