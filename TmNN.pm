#!/usr/local/bin/perl

package main; 
our $SEE = 0;


package TmNN;

use strict;
use warnings;
use Carp;

my $SALT_CONCENTRATION = 50e-3; # 50 mM Na
my $PROBE_CONCENTRATION = 1e-6; # 1 uM 


# public static methods:
####
sub set_SALT_concentration {
    my ($salt) = @_;
    unless ($salt =~ /\d/) {
        confess "Error, salt concentration must be numeric.\n";
    }
    
    $SALT_CONCENTRATION = $salt;
    return;
}

####
sub set_PROBE_concentration {
    my ($probe) = @_;
    unless ($probe =~ /\d/) { 
        confess "Error, probe concentration must be numeric.\n";
    }
    
    $PROBE_CONCENTRATION = $probe;
    return;
}


## below are object methods.

#############################################################################################
# Energies are the unified values provided in:
# SantaLucia PNAS 1998;95;1460-1465
# A unified view of polymer, dumbbell, and oligonuclitide DNA nearest-neighbor thermodynamics
##############################################################################################


my %deltaH = ( AA => -7.9, TT => -7.9,
               AT => -7.2,
               TA => -7.2,
               CA => -8.5, TG => -8.5,
               GT => -8.4, AC => -8.4,
               CT => -7.8, AG => -7.8,
               GA => -8.2, TC => -8.2,
               CG => -10.6, 
               GC => -9.8,
               GG => -8.0, CC => -8.0 );

my %deltaS = ( AA => -22.2, TT => -22.2,
               AT => -20.4,
               TA => -21.3,
               CA => -22.7, TG => -22.7,
               GT => -22.4, AC => -22.4,
               CT => -21.0, AG => -21.0,
               GA => -22.2, TC => -22.2,
               CG => -27.2,
               GC => -24.4,
               GG => -19.9, CC => -19.9,
    );


my %deltaG = ( AA => -1.00, TT => -1.00,
               AT => -0.88,
               TA => -0.58,
               CA => -1.45, TG => -1.45,
               GT => -1.44, AC => -1.44,
               CT => -1.28, AG => -1.28,
               GA => -1.30, TC => -1.30,
               CG => -2.17,
               GC => -2.24,
               GG => -1.84, CC => -1.84,
    );


my $TEMPERATURE_K = 310.15; # room temp in Kelvin

my $R = 1.987; # cal / K mol;  natural gas constant.

my $init_w_term_GC_deltaH = 0.1;
my $init_w_term_GC_deltaS = -2.8;
my $init_w_term_GC_deltaG = 0.98;

my $init_w_term_AT_deltaH = 2.3;
my $init_w_term_AT_deltaS = 4.1;
my $init_w_term_AT_deltaG = 1.03;



sub calcTm {
    my ($seq_record) = @_;
    
    my $sequence = $seq_record->{sequence} or die "Error, no sequence provided";
    my $probe_concentration = $seq_record->{probe_concentration} || $PROBE_CONCENTRATION; # 1 uM default
    my $salt_concentration = $seq_record->{salt_concentration} || $SALT_CONCENTRATION; # 50 mM default. 
    

    my @seq_chars = split (//, uc $sequence);
    
    my $sum_deltaH = 0;
    my $sum_deltaS = 0;
    my $sum_deltaG = 0;
    
    for (my $i = 0; $i < $#seq_chars; $i++) {
        my $dinuc = $seq_chars[$i] . $seq_chars[$i+1];
        my $dH = $deltaH{$dinuc} or die "Error, no deltaH for $dinuc";
        my $dS = $deltaS{$dinuc} or die "Error, no deltaS for $dinuc";
        my $dG = $deltaG{$dinuc} or die "Error, no deltaG for $dinuc";
        
        $sum_deltaH += $dH;
        $sum_deltaS += $dS;
        $sum_deltaG += $dG;
    }
    
    ## check first base pair
    if ($seq_chars[0] =~ /GC/) {
        $sum_deltaH += $init_w_term_GC_deltaH;
        $sum_deltaS += $init_w_term_GC_deltaS;
        $sum_deltaG += $init_w_term_GC_deltaG;
    }
    else {
        $sum_deltaH += $init_w_term_AT_deltaH;
        $sum_deltaS += $init_w_term_AT_deltaS;
        $sum_deltaG += $init_w_term_AT_deltaG;
    }

    ## check last base pair
    if ($seq_chars[$#seq_chars] =~ /GC/) {
        $sum_deltaH += $init_w_term_GC_deltaH;
        $sum_deltaS += $init_w_term_GC_deltaS;
        $sum_deltaG += $init_w_term_GC_deltaG;
    }
    else {
        $sum_deltaH += $init_w_term_AT_deltaH;
        $sum_deltaS += $init_w_term_AT_deltaS;
        $sum_deltaG += $init_w_term_AT_deltaG;
    }
            
    my $Tm = ( 1000*$sum_deltaH / ($sum_deltaS + $R * log($probe_concentration/4.0))) - 273.15 
        + 12.0 * (log($salt_concentration) / log(10)); 
    
    my $re_dG = sprintf ("%.2f", &recalc_deltaG($sum_deltaH, $sum_deltaS));
    my $re_dH = sprintf ("%.2f", &recalc_deltaH($sum_deltaS, $sum_deltaG));
    my $re_dS = sprintf ("%.2f", &recalc_deltaS($sum_deltaH, $sum_deltaG));
    
    if ($SEE) {
        print "dG = $sum_deltaG, recalc: $re_dG\n";
        print "dH = $sum_deltaH, recalc: $re_dH\n";
        print "dS = $sum_deltaS, recalc: $re_dS\n";
        
        print "dG: $sum_deltaG\tdH: $sum_deltaH\tdS: $sum_deltaS\tTm: $Tm\n";
    }

    
    return (sprintf("%.2f", $Tm));
    
}


####
sub compare_deltaGHS {
    
    foreach my $dinuc (keys %deltaG) {
        my $dG = $deltaG{$dinuc};
        my $dH = $deltaH{$dinuc};
        my $dS = $deltaS{$dinuc};

        my $recalc_dG = sprintf ("%.2f", &recalc_deltaG($dH, $dS));
        
        my $recalc_dH = sprintf ("%.2f", &recalc_deltaH($dS, $dG)); 
        
        my $recalc_dS = sprintf ("%.2f", &recalc_deltaS($dH, $dG));
        
        print "$dinuc\tdG: $dG, $recalc_dG\tdH: $dH, $recalc_dH\tdS: $dS, $recalc_dS\n";
    }
}


####
sub recalc_deltaG {
    my ($dH, $dS) = @_;
    return ($dH - $TEMPERATURE_K * $dS/1000);
}

####
sub recalc_deltaH {
    my ($dS, $dG) = @_;
    return ($dG + $TEMPERATURE_K * $dS/1000);
}

####
sub recalc_deltaS {
    my ($dH, $dG) = @_;
    return ( ($dH - $dG) / $TEMPERATURE_K * 1000);
}


## Testing:
package main;

use strict;
use warnings;
use File::Basename;

if (basename($0) eq "TmNN.pm") {
    
    # &TmNN::compare_deltaGHS(); exit(0);
    
    my $usage = "usage: $0 sequence probeConc saltConc\n\n";
    my $sequence = $ARGV[0] or die $usage;
    my $probe_conc = $ARGV[1] or die $usage;
    my $salt_conc = $ARGV[2] or die $usage;
    
    my $Tm = &TmNN::calcTm( { sequence => $sequence,
                              probe_concentration => $probe_conc,
                              salt_concentration => $salt_conc,
                          }
                            );
    
    print "Tm($sequence) = $Tm\n";


    


    exit(0);
}


1;
