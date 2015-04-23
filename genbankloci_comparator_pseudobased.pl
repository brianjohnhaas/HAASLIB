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

$/ = "\n";

my $SEE = $ARGV[0];
my $DEBUG = 0;


my %clone_name_via_asmbl_id;

my $query = "select asmbl_id, clone_name from clone_info";
my @results = &do_sql ($dbproc, $query);
foreach my $result (@results) {
    my ($asmbl_id, $clone_name) = split (/\t/, $result);
    $clone_name_via_asmbl_id{$asmbl_id} = $clone_name;
}



my ($absolute_total, $absolute_unequal, $absolute_equal, $absolute_static) = (0,0,0,0);

my (%TIGR_struct, #stores all tigr feat_names as keys.  'assigned' == 1 if assigned, == 0 if not assigned.
    %Genbank_loci, #stores all bac-based loci present within this bac. Basically, a quick lookup table.
    %Genbank_struct); #stores loci keyed by clone_name for comparison with tigr struct.

 main: {
     my $x = 0; 
     while (<STDIN>) {
	 ## Process data one bac at a time.
	 my ($clone_name, $gb_acc, $asmbl_id);
	 foreach my $line (split (/\n/)) {
	     if ($line =~ /^GENBANK_INFO/) {
		 #print "$line\n\n";
		 &process_genbank_info($line, \%Genbank_struct);
	     } elsif ($line =~ /^TIGR_feat/) {
		 #print "$line\n\n";
		 
		 &process_tigr_info ($line, \%TIGR_struct);
		 $x++;
		 #print Dumper(\%TIGR_struct);
	     }
	 }
	 #if ($x > 10) {last;}
     }
     &analyze_data (\%TIGR_struct, \%Genbank_struct);
     print "***************\nComplete Analsis:\n";
     print "Total: $absolute_total\nEqual: $absolute_equal\nNum static: $absolute_static\nNum to update: $absolute_unequal\n";    
 }


$dbproc->disconnect;


####
sub analyze_data {
    my ($T_struct_ref, $G_struct_ref) = @_;
    ## local copies of data structures.
    

    foreach my $clone_name (keys %$T_struct_ref) {
	my %TIGR_struct = %{$T_struct_ref->{$clone_name}};
	my %Genbank_loci = ();
	if (exists ($G_struct_ref->{$clone_name})) {
	    %Genbank_loci = %{$G_struct_ref->{$clone_name}};
	}
	my (%T_LOCI_ASSIGNED); #holds tigr loci that have been assigned.
	my %Preexisting_TIGR_loci; #holds initial tigr loci assignments.

	print "GENBANK LOCI for $clone_name\n" . Dumper($G_struct_ref->{$clone_name}) if $SEE;
	
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
		print "\tCASE 0: correctly assigned genes. t: $tigr_locus g: $genbank_locus\n" if $SEE;
		$TIGR_struct{$tigr_gene_feat}->{assigned} = 1;
		$T_LOCI_ASSIGNED{$tigr_locus} = 1;
		$TIGR_struct{$tigr_gene_feat}->{case} = 0;
		#print Dumper($TIGR_struct{$tigr_gene_feat}) if $SEE;
	    }
	}
	
	## Now for the ones that were not assigned correctly.  
	## If genbank locus exists and it hasn't already been assigned to another tigr gene, assign it.
	
	foreach my $tigr_gene_feat (keys %TIGR_struct) {
	    if ($TIGR_struct{$tigr_gene_feat}->{assigned}) {next;}
	    my $tigr_locus = $TIGR_struct{$tigr_gene_feat}->{t_locus};
	    my $genbank_locus = $TIGR_struct{$tigr_gene_feat}->{g_locus};
	    if ( ($genbank_locus) && 
		 (!$T_LOCI_ASSIGNED{$genbank_locus}) ){
		print "\tCASE 1: tigr obtains genbank locus. t: $tigr_locus g: $genbank_locus\n" if $SEE;
		$TIGR_struct{$tigr_gene_feat}->{assigned} = 1;
		$T_LOCI_ASSIGNED{$genbank_locus} = 1;
		$TIGR_struct{$tigr_gene_feat}->{new_locus} = $genbank_locus;
		$TIGR_struct{$tigr_gene_feat}->{case} = 1;
		#print Dumper($TIGR_struct{$tigr_gene_feat}) if $SEE;
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
		 (!$Genbank_loci{$tigr_locus}) ){
		print "\tCASE 2: tigr loci exists, not already assigned, and not a preexisting genbank locus\n" if $SEE;
		$TIGR_struct{$tigr_gene_feat}->{assigned} = 1;
		$TIGR_struct{$tigr_gene_feat}->{case} = 2;
		$T_LOCI_ASSIGNED{$tigr_locus} = 1;
		print "\tt: $tigr_locus\n" if $SEE; 
		#print Dumper($TIGR_struct{$tigr_gene_feat}) if $SEE;
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
		print "\tt: $tigr_locus new: $new_locus\n" if $SEE;
		#print Dumper($TIGR_struct{$tigr_gene_feat}) if $SEE;
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
	    my $new_locus;
	    if ($new_locus = $gene_data_ref->{new_locus}) {
		$num_unequal++; ## some db update is required
		print $tigr_gene_feat . "\n" . Dumper ($gene_data_ref) if $DEBUG;
		&update_gene_info($tigr_gene_feat, $new_locus);
	    } elsif (!$TIGR_struct{$tigr_gene_feat}->{case}) {
		$num_equal++; ## perfect preexisting loci match
	    } else {
		$num_static++; ## apparently, already assigned tigr locus is acceptable.
	    }
	}
	
	print "Total: $total_num\nEqual: $num_equal\nNum Stati: $num_static\nNum to update: $num_unequal\n" if $SEE;
	$absolute_total += $total_num;
	$absolute_unequal += $num_unequal;
	$absolute_static += $num_static;
	$absolute_equal += $num_equal;
	
    }

    ## print out comparative data
    foreach my $clone_name (keys %$T_struct_ref) {
	print "//CLONE: $clone_name\n";
	foreach my $tigr_feat_name (keys %{$T_struct_ref->{$clone_name}}) {
	    my $tigr_locus = $T_struct_ref->{$clone_name}->{$tigr_feat_name}->{new_locus};
	    unless ($tigr_locus) {
		 $tigr_locus = $T_struct_ref->{$clone_name}->{$tigr_feat_name}->{t_locus};
	     }
	    my $genbank_text = $T_struct_ref->{$clone_name}->{$tigr_feat_name}->{g_text};
	    my $gen_locus =  $T_struct_ref->{$clone_name}->{$tigr_feat_name}->{g_locus};
	    unless ($gen_locus) {
		$gen_locus = "NONE";
	    }
	    print "$tigr_feat_name\t$tigr_locus\t$gen_locus\t$genbank_text\n";
	}
    }
    
}

####
sub update_gene_info {
    my ($tigr_gene_feat, $new_locus) = @_;
    my $query = "update ident set locus = \"$new_locus\" where feat_name = \"$tigr_gene_feat\"\n\n";
    print $query;
}


#### 
sub process_genbank_info {
    my ($line, $Gen_struct_ref) = @_;
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
	$Gen_struct_ref->{$clone_name}->{$locus} = 1;
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
    my ($line, $tstruct_ref) = @_;
    my @x = split (/\t/, $line);
    my ($tigr_feat, $t_loc, $t_ploc, @genbank_text) = ($x[1], $x[3], $x[5], $x[7], $x[9], $x[10]);
    my @chromo_syns = &get_chromo_syns($tigr_feat);
    foreach my $chromo_syn (@chromo_syns) {
	my $true_tigr_locus = &get_locus($chromo_syn);
	$true_tigr_locus =~ tr/_/\./;
	my $clone_name = &get_clone_name($chromo_syn);
	my $g_locus = &findme_loci($clone_name, \@genbank_text);
	$tstruct_ref->{$clone_name}->{$chromo_syn}->{t_locus} = $true_tigr_locus;
	if ($g_locus =~ /$clone_name/) {
	    $tstruct_ref->{$clone_name}->{$chromo_syn}->{g_locus} = $g_locus;
	}
	my $genbank_text = join (" ", @genbank_text);
	$tstruct_ref->{$clone_name}->{$chromo_syn}->{g_text} .= $genbank_text;
    }
}


####
sub get_locus {
    my ($tu_feat_name) = @_;
    my $query = "select locus from ident where feat_name = \"$tu_feat_name\"\n";
    return (&first_result_sql ($dbproc, $query));
}


#### 
sub get_clone_name {
    my ($feat_name) = @_;
    $feat_name =~ /^(\d+)\./;
    my $asmbl_id = $1;
    if (!$asmbl_id) { die "Can't determine asmbl_id based on feat_name $feat_name\n";}
    return ($clone_name_via_asmbl_id{$asmbl_id});
}




####
sub get_chromo_syns {
    my ($tigr_feat) = @_;
    my $query = "select syn_feat_name from chromo_syn_link where feat_name = \"$tigr_feat\"\n";
    return (&do_sql ($dbproc, $query));
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











