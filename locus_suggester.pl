#!/usr/local/bin/perl


use lib ($ENV{EUK_MODULES}, $ENV{EGC_SCRIPTS});
use Egc_library;
use strict;
use DBI;
use Loci_assigner;
use Data::Dumper;


my $usage = "usage: $0 db asmbl_id_list_file [locus|pub_locus] preexistingLociFile\n\n";

my $db = $ARGV[0] or die $usage;
my $asmbl_id_file_list = $ARGV[1] or die $usage;
my $locus_type = $ARGV[2] or die $usage;
my $preexisting_loci_file = $ARGV[3];

my ($dbproc) = &ConnectToDb("SYBTIGR","Sybase","access","access","$db");

#############################################################################
# begin program

my %preexistingLoci;

if ($preexisting_loci_file) {
    open (FILE, $preexisting_loci_file) or die "Error, cannot read $preexisting_loci_file\n\n";
    while (<FILE>) {
	my $locus = $_;
	$locus =~ s/\s//g;
	if ($locus) {
	    $preexistingLoci{$locus} = 1;
	}
    }
    close FILE;
}



open (FILE, $asmbl_id_file_list) or die "Error, cannot open file $asmbl_id_file_list\n";
while (<FILE>) {
    my $asmbl_id = $_;
    $asmbl_id =~ s/\s//g;
    
    my $query = "select i.feat_name, i.$locus_type, a.end5, a.end3 from ident i, asm_feature a where a.feat_type = 'TU' and a.feat_name = i.feat_name and a.asmbl_id = $asmbl_id";
    
    my @genes;

    my @results = &do_sql_2D($dbproc, $query);
    foreach my $result (@results) {
	my ($feat_name, $locus, $end5, $end3) = @$result;
	unless ($locus =~ /\w/) {
	    $locus = undef;
	}
	push (@genes, { feat_name => $feat_name,
			locus => $locus,
			end5 => $end5,
			end3 => $end3 } );
    }
    
    eval {
    
	@genes = &Loci_assigner::assign_loci(\@genes, \%preexistingLoci);
	
	print "\n\n#############################################\nLOCI for asmbl_id: $asmbl_id:\n";
	foreach my $gene (@genes) {
	    my ($feat_name, $locus, $end5, $end3) = ($gene->{feat_name},
						     $gene->{locus},
						     $gene->{end5},
						     $gene->{end3});
	    
	    my ($new_locus) = $gene->{newLocus};
	    print "$feat_name\t$end5\t$end3\t$locus\t";
	    if ($new_locus) {
		print "****\t$new_locus\n";
	    } else {
		print "\n";
	    }
	}
    };

    if ($@) {
	open (ERROR, ">$asmbl_id.$db.loci_suggestions.error") or die "Cannot open error file for $asmbl_id, $db.\n";
	print ERROR Dumper (\@genes);
    }
}


