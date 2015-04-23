#!/usr/local/bin/perl

## Requirements
## bac-based loci must be unique for tigr genes.
## if tigr loci match genbank loci, the genes must also match each other.
## already assigned tigr genes will preferentially be used if it meets the above constraints.

use strict;
use Data::Dumper;
use vars qw ($SEE);

use DBI;
use lib ($ENV{EGC_SCRIPTS});
use Egc_library;


my ($dbproc) = &ConnectToDb("SYBTIGR","Sybase","access","access","ath1");

$/ = "\nEND";

my $SEE = $ARGV[0];
my $DEBUG = 0;

my ($absolute_static, $absolute_total, $absolute_unequal, $absolute_equal) = (0,0,0);

main: {
    while (<STDIN>) {
	## Process data one bac at a time.
	my (%TIGR_struct, #stores all tigr feat_names as keys.  'assigned' == 1 if assigned, == 0 if not assigned.
	    %Genbank_loci); #stores all bac-based loci present within this bac.
	
	my ($clone_name, $gb_acc, $asmbl_id);
	foreach my $line (split (/\n/)) {
	    if ($line =~ /^BEGIN/) {
		#print "$line\n\n";
		
		($gb_acc, $clone_name, $asmbl_id) = &process_BEGIN_statement($line);
	    } elsif ($line =~ /^GENBANK_INFO/) {
		#print "$line\n\n";
		&process_genbank_info($line, \%Genbank_loci);
	    } elsif ($line =~ /^TIGR_feat/) {
		#print "$line\n\n";
		
		&process_tigr_info ($line, $clone_name, \%TIGR_struct);
		#print Dumper(\%TIGR_struct);
	    }
	}
	
	if ( (!$clone_name) && (%TIGR_struct)) { die "No clone_name here!!\n";}
	&analyze_data (\%TIGR_struct, \%Genbank_loci, $clone_name);
	#die;
    }
    print "***************\nComplete Analsis:\n";
    print "Total: $absolute_total\nEqual: $absolute_equal\nNum static: $absolute_static\nNum to update: $absolute_unequal\n"; 
}



$dbproc->disconnect;


####
sub analyze_data {
    my ($T_struct_ref, $G_loci_ref, $clone_name) = @_;
    ## local copies of data structures.
    my %TIGR_struct = %$T_struct_ref;
    my %Genbank_loci = %$G_loci_ref;
    my (%T_LOCI_ASSIGNED); #holds tigr loci that have been assigned.
    my %Preexisting_TIGR_loci; #holds initial tigr loci assignments.

    print "GENBANK LOCI for $clone_name\n" . Dumper($G_loci_ref) if $SEE;

    ## First Step of Analysis
    ##     -if tigr loci matches genbank loci, keep it if that loci has not 
    ##      already been assigned to another tigr gene.
    foreach my $tigr_gene_feat (keys %TIGR_struct) {
	my $tigr_locus = $TIGR_struct{$tigr_gene_feat}->{t_locus};
	$Preexisting_TIGR_loci{$tigr_locus} = 1;
	my $genbank_locus = $TIGR_struct{$tigr_gene_feat}->{g_locus};
	if ( ($tigr_locus) && ($genbank_locus) && 
	     ($tigr_locus eq $genbank_locus) &&
	     (!$T_LOCI_ASSIGNED{$tigr_locus}) ){
	    $TIGR_struct{$tigr_gene_feat}->{assigned} = 1;
	    $T_LOCI_ASSIGNED{$tigr_locus} = 1;
	    $TIGR_struct{$tigr_gene_feat}->{case} = 0;
	    print "CASE 0: IDENTICAL ASSIGNMENTS: t: $tigr_locus g: $genbank_locus\n" if $SEE;
	}
    }

    ## Now for the ones that were not assigned correctly.  
    ## If genbank locus exists and it hasn't already been assigned to another tigr gene, assign it.
       
    foreach my $tigr_gene_feat (keys %TIGR_struct) {
	if ($TIGR_struct{$tigr_gene_feat}->{assigned}) {next;}
	my $tigr_locus = $TIGR_struct{$tigr_gene_feat}->{t_locus};
	my $genbank_locus = $TIGR_struct{$tigr_gene_feat}->{g_locus};
	if ( ($genbank_locus) && 
	     (!$T_LOCI_ASSIGNED{$genbank_locus}) && !$Preexisting_TIGR_loci{$genbank_locus}){
	    print "\tCASE 1: TIGR Gene being assigned to genbank locus: t: $tigr_locus g: $genbank_locus\n" if $SEE;
	    $TIGR_struct{$tigr_gene_feat}->{assigned} = 1;
	    $T_LOCI_ASSIGNED{$genbank_locus} = 1;
	    $TIGR_struct{$tigr_gene_feat}->{new_locus} = $genbank_locus;
	    $TIGR_struct{$tigr_gene_feat}->{case} = 1;
	    print Dumper($TIGR_struct{$tigr_gene_feat}) if $SEE;
	}
    }

    ## For remaining tigr genes, if tigr locus exists and is not already assigned to another tigr gene
    ## and is not a current genbank locus, keep it.  Otherwise, give it a new locus
    my $x = 1;
    foreach my $tigr_gene_feat (keys %TIGR_struct) {
	if ($TIGR_struct{$tigr_gene_feat}->{assigned}) {next;}
	my $tigr_locus = $TIGR_struct{$tigr_gene_feat}->{t_locus};
	if ( ($tigr_locus) && 
	     (!$T_LOCI_ASSIGNED{$tigr_locus}) &&
	     (!$Genbank_loci{$tigr_locus}) &&
	     ($tigr_locus =~ /$clone_name/i) ) {
	    print "\tCASE 2: tigr loci exists, not already assigned, and not a preexisting genbank locus\n\tt: $tigr_locus\n" if $SEE;
	    $TIGR_struct{$tigr_gene_feat}->{assigned} = 1;
	    $TIGR_struct{$tigr_gene_feat}->{case} = 2;
	    $T_LOCI_ASSIGNED{$tigr_locus} = 1;
	    print Dumper($TIGR_struct{$tigr_gene_feat}) if $SEE;
	} else {
	    print "\tCASE 3: Current tigr locus requires replacement\n" if $SEE;
	    ## give tigr gene a new assignment that isn't already extant.
	    my $new_locus = $clone_name . ".$x";
	    while ( $T_LOCI_ASSIGNED{$new_locus} || $Genbank_loci{$new_locus} || $Preexisting_TIGR_loci{$new_locus} ) {
		$x++;
		$new_locus = $clone_name . ".$x";
	    }
	    $TIGR_struct{$tigr_gene_feat}->{assigned} = 1;
	    $T_LOCI_ASSIGNED{$new_locus} = 1;
	    $TIGR_struct{$tigr_gene_feat}->{new_locus} = $new_locus;
	    $TIGR_struct{$tigr_gene_feat}->{case} = 3;
	    print "\tt: $tigr_locus\n" if $SEE;
	    print Dumper($TIGR_struct{$tigr_gene_feat}) if $SEE;
	}
    }
    #print Dumper(\%TIGR_struct);
    
    ## Now, assignments have been made and updates to db are required.
    my $num_unequal = 0;
    my $total_num = 0;
    my $num_equal = 0;
    my $num_static = 0; 

    foreach my $tigr_gene_feat (keys %TIGR_struct) {
	$total_num++;
	my $gene_data_ref = $TIGR_struct{$tigr_gene_feat};
	print "REPORT:\t$tigr_gene_feat\t$gene_data_ref->{t_locus}\t$gene_data_ref->{g_locus}\n";
	my $new_locus;
	if ($new_locus = $gene_data_ref->{new_locus}) {
	    $num_unequal++;
	    print $tigr_gene_feat . "\n" . Dumper ($gene_data_ref) if $DEBUG;
	    &update_gene_info($tigr_gene_feat, $new_locus);
	} elsif (!$TIGR_struct{$tigr_gene_feat}->{case}) {
	    $num_equal++; ## perfect preexisting loci match
	} else {
	    $num_static++; ## apparently, already assigned tigr locus is acceptable.
	}
    }
    
    print "Total: $total_num\nEqual: $num_equal\nDiff: $num_unequal\n" if $SEE;
    $absolute_total += $total_num;
    $absolute_unequal += $num_unequal;
    $absolute_equal += $num_equal;
    $absolute_static += $num_static;


    ## print out comparative data
    print "//CLONE: $clone_name\n";
    foreach my $tigr_feat_name (keys  %TIGR_struct) {
	my $tigr_locus = $TIGR_struct{$tigr_feat_name}->{new_locus};
	unless ($tigr_locus) {
	    $tigr_locus = $TIGR_struct{$tigr_feat_name}->{t_locus};
	}
	my $genbank_text = $TIGR_struct{$tigr_feat_name}->{g_text};
	my $gen_locus =  $TIGR_struct{$tigr_feat_name}->{g_locus};
	unless ($gen_locus) {
	    $gen_locus = "NONE";
	}
	print "$tigr_feat_name\t$tigr_locus\t$gen_locus\t$genbank_text\n";
    }

}


####
sub update_gene_info {
    my ($tigr_gene_feat, $new_locus) = @_;
    ## Only updating bac_locus if doesn't already exist.
    my $query = "select locus from ident where feat_name = \"$tigr_gene_feat\"";
    my $result = &first_result_sql($dbproc, $query);
    unless ($result =~ /\w/) {
	    my $query = "update ident set locus = \"$new_locus\" where feat_name = \"$tigr_gene_feat\"\n\n";
	    print $query;
	}
}


####
sub process_BEGIN_statement {
    my ($line) = @_;
    my @x = split (/\t/, $line);
    shift @x;
    print "\n\n*** Processed BEGIN info: [@x]\n" if $SEE;
    return (@x);
}


#### 
sub process_genbank_info {
    my ($line, $Gen_loc_ref) = @_;
    my @data = split (/\t/, $line);
    shift @data; #rid genbank_info identifier
    my $clone_name = shift @data;
    shift @data;
    my @texts;
    push (@texts,  shift @data);
    shift @data;
    push (@texts, @data);
    my $locus = &findme_loci($clone_name, \@texts);
    if ($locus) {
	&add_locus($locus, $Gen_loc_ref);
    }
}


####
sub add_locus {
    my ($locus, $hash_ref) = @_;
    $locus =~ tr/_/\./;
    print "Adding locus: $locus\n" if $DEBUG;
    $hash_ref->{$locus} = 1;
}


####
sub process_tigr_info {
    my ($line, $clone_name, $tstruct_ref) = @_;
    my @x = split (/\t/, $line);
    my ($tigr_feat, $t_loc, $t_ploc, @genbank_text) = ($x[1], $x[3], $x[5], $x[7], $x[9], $x[10]);
    my @tigr_text = ($t_loc, $t_ploc);
    my $t_locus = &get_locus ($tigr_feat, $clone_name);
    $t_locus =~ tr/_/./;
    my $genbank_text = join (" ", @genbank_text);
    
    #my $t_locus = &findme_loci($clone_name, \@tigr_text);
    my $g_locus = &findme_loci($clone_name, \@genbank_text);
    $tstruct_ref->{$tigr_feat} = { t_locus => $t_locus,
				   g_locus => $g_locus};  
    $tstruct_ref->{$tigr_feat}->{g_text} .= $genbank_text;
}


####
sub findme_loci {
    my ($clone_name, $text_ref) = @_;
    foreach my $text (@$text_ref) {
	print "\tLooking for loci matching clone: $clone_name in [$text]\n" if $DEBUG;
	while ($text =~ /(($clone_name)[\._]\d+)/g) {
	    my $locus = $1;
	    $locus =~ tr/_/\./;
	    print "\t\tFound [$locus]\n" if $DEBUG;
	    return ($locus);
	}
    }
} 





####
sub get_locus {
    my ($tu_feat_name, $clone_name) = @_;
    my $query = "select locus from ident where feat_name = \"$tu_feat_name\"\n";
    my $db_locus =&first_result_sql ($dbproc, $query);
    if ($db_locus =~ /$clone_name/i) {
	return ($db_locus);
    } else {
	return ("");
    }
}











