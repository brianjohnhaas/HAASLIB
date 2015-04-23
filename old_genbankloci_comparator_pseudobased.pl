#!/usr/local/bin/perl

use strict;
use lib ($ENV{EGC_SCRIPTS});
use Egc_library;
use DBI;

$| = 1;

my $QUERY = @ARGV;

my $dbproc = &ConnectToDb("SYBTIGR","Sybase","bhaas","tim1bird","ath1");

open (QUERIES, ">queries") if $QUERY;

my %loci_tracker; #ensures unique loci assignments
my %query_tracker; #track the 1-1 assignments for generating queries.
my %feat_name_tracker; #ensures one gene gets one loci assignment
my %problem_loci_tracker; #holds multiple genes with same bac locus ( gene->locus )
my %problem_gene_tracker; #holds multiple loci assigned to same gene (locus->gene)
my $number_equal = 0; #number of loci identical
my $number_unequal = 0; #number of loci not identical



while (<STDIN>) {
    unless (/\w/) {next;}
    my @x = split (/\s+/);
    my ($tigr_feat, $t_loc, $t_ploc, $g_loc, $g_ploc, $g_comname) = ($x[1], $x[3], $x[5], $x[7], $x[9], $x[10]);
    $g_loc .= ":${g_ploc}:${g_comname}"; #look at everything just in case.
    my @feat_names = (split (/\:/, $tigr_feat));
    my @tigr_loci =  (split (/\:/, $t_loc));
    if ($#feat_names > 1) { #more than two tigr genes matching a single genbank entry
	print "Too many TIGR genes matched by genbank entry : $g_loc, $g_ploc, $g_comname\n";
	next;
    }
    foreach my $feat_name (@feat_names) {
	my $tigr_bac_locus = shift (@tigr_loci);
	my $query = "select syn_feat_name from chromo_syn_link where feat_name = \"$feat_name\"\n";
	my $bac_feat = &first_result_sql ($dbproc, $query);
	my $asmbl_id;
	if ($bac_feat =~ /^(\d+)\./) {
	    $asmbl_id = $1;
	} else {
	    print  "Can't determine asmbl_id based on feat_name $bac_feat, chromofeat: $feat_name\n";
	    next;
	}
	my $clone_name = &get_clone_name($asmbl_id);
	if ($tigr_bac_locus !~ /$clone_name/) {
	    print  "$tigr_bac_locus of $bac_feat doesn't include clone_name $clone_name\n";
	}
	
	
	
	## Here is where every genbank locus field is compared to the TIGR gene.
	
	my %seen; #look at each gb_bac_locus in $g_loc only once.
	foreach my $gb_bac_locus (split (/\:/, $g_loc)) { #not sure where it is (locus:publocus:com_name).
	    $gb_bac_locus =~ s/>//;
	    $gb_bac_locus =~ s/\~.*//;
	     if ($seen{$gb_bac_locus}) {
		 next;
	     } 
	     $seen{$gb_bac_locus} = 1;
	     if ($gb_bac_locus =~ /$clone_name/) {
		 print "$bac_feat\t$tigr_bac_locus\t$gb_bac_locus\n";
		 ## data checking. Acknowledge multiple genes with same bac locus.
		 if ( ($loci_tracker{$gb_bac_locus}) && ($bac_feat ne $loci_tracker{$gb_bac_locus})) {
		     print "Problem, multiple feats assigned same bac_locus: gb:$gb_bac_locus\tTIGRfeat:$bac_feat\tAlready assigned:$loci_tracker{$gb_bac_locus}\n";
		     delete ($query_tracker{$gb_bac_locus});
		     $problem_loci_tracker{$bac_feat} = $gb_bac_locus;
		     $problem_loci_tracker{$loci_tracker{$gb_bac_locus}} = $gb_bac_locus;
		     next;
		 } else {
		     $loci_tracker{$gb_bac_locus} = $bac_feat;
		 }
		 ## data checking. Make sure the gene hasn't already been assigned a locus
		 if ( ($feat_name_tracker{$bac_feat}) && ($feat_name_tracker{$bac_feat} ne $gb_bac_locus)) {
		     print "Problem with locus $gb_bac_locus: gene $bac_feat has already been assigned a locus $feat_name_tracker{$bac_feat}\n";
		     delete  ($query_tracker{$gb_bac_locus});
		     $problem_gene_tracker{$gb_bac_locus} = $bac_feat;
		     $problem_gene_tracker{$feat_name_tracker{$bac_feat}} = $bac_feat;
		     next;
		 } else {
		     $feat_name_tracker{$bac_feat} = $gb_bac_locus;
		 }
		 my $temp_tigr_bac_locus = $tigr_bac_locus;
		 $temp_tigr_bac_locus =~ s/_/\./;
		 my $temp_gb_bac_locus = $gb_bac_locus;
		 $temp_gb_bac_locus =~ s/_/\./;
		 if ($temp_gb_bac_locus eq $temp_tigr_bac_locus) {
		     $number_equal++;
		 } else {
		     $number_unequal++;
		 }
		 if ($tigr_bac_locus ne $gb_bac_locus) {		     
		     $query_tracker{$gb_bac_locus} = $bac_feat;
		 }
	     }
	 }
     }
}

$dbproc->disconnect;
   
## print out query data 

print "////// Number equal: $number_equal, Number conflict: $number_unequal\n";


if ($QUERY) {
    foreach my $locus (keys %query_tracker) {
	my $bac_feat = $query_tracker{$locus};
	print QUERIES "update ident set locus = \"$locus\" , date = getdate() where feat_name = \"$bac_feat\"\n\n" if $QUERY;
    }
}

## print out problem data sets:

print "\n\n#####\nMultiple genes assigned to same loci:\n";
foreach my $gene (sort {$problem_loci_tracker{$a} cmp $problem_loci_tracker{$b}} keys %problem_loci_tracker) {
    print "$problem_loci_tracker{$gene}\t$gene\n";
}

print "\n\n#####\nMultiple loci assigned to same gene\n";
foreach my $locus (sort {$problem_gene_tracker{$a} cmp $problem_gene_tracker{$b}} keys %problem_gene_tracker) {
    print "$locus\t$problem_gene_tracker{$locus}\n";
}

exit;




####
sub get_clone_name {
    my ($asmbl_id) = @_;
    my $query = "select clone_name from clone_info where asmbl_id = $asmbl_id\n";
    my $clone_name = &first_result_sql ($dbproc, $query);
    return ($clone_name);
}












