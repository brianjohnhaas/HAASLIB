package Run_Align;

use strict;
use Cwd;
require Exporter;

our @ISA = qw (Exporter);
our @EXPORT = qw (align align0);

sub align {
    my @fastaSeqs = @_; ## two raw peptide sequences
    my $prog = "align";
    
    return (&runAlignment($prog, @fastaSeqs));
}


sub align0 {
    my @fastaSeqs = @_;
    my $prog = "align0";
    
    return (&runAlignment($prog, @fastaSeqs));
}



## private
sub runAlignment {
    my ($prog, $seq1, $seq2) = @_;
    my $host = $ENV{HOST};
    my $tmpDir = "/tmp/align.$host.$$";
    mkdir ($tmpDir) or die "Cannot make temp dir. $tmpDir\n";
    my $cwd = cwd;
    
    chdir $tmpDir or die "Cannot cd to $tmpDir\n";
    
    my $tmpFile1 = "prot_1";
    my $tmpFile2 = "prot_2";
    my $outFile = "$prog.output";
    
    my $len1 = length($seq1);
    my $len2 = length($seq2);
    my $maxLen = ($len1 < $len2) ? $len1 : $len2;
    
    open (FILE, ">$tmpFile1") or die "Error, cannot open file $tmpFile1\n";
    print FILE ">prot_1\n$seq1\n";
    close FILE;

    open (FILE, ">$tmpFile2") or die "Error, cannot open file $tmpFile2\n";
    print FILE ">prot_2\n$seq2\n";
    close FILE;
    
    my $cmd = "$prog $tmpFile1 $tmpFile2 > $outFile";
    my $ret = system ($cmd);
    if ($ret) {
	# die "Error, cmd: $cmd\n returned exit status ($ret)\n";
	## apparently, cannot trust the return value.
    }
    my @seqArray1;
    my @seqArray2;
    open (FILE, $outFile);
    while (<FILE>) {
	chomp;
	if (/^prot_(\d)\s+([A-Z\-]+)/i) {
	    my $seqNum = $1;
	    my $sequence = $2;
	    if ($seqNum == 1) {
		push (@seqArray1, split (//, $sequence));
	    } elsif ($seqNum == 2) {
		push (@seqArray2, split (//, $sequence));
	    }
	} 
    }
    close FILE;
    
    unless (@seqArray1 && @seqArray2) {
	die "Error, no alignments read from $tmpDir/$outFile\n";
    }
    
    unlink ($tmpFile1);
    unlink ($tmpFile2);
    unlink ($outFile);
    chdir ($cwd) or die "Cannot get back to $cwd\n";
    rmdir ($tmpDir) or die "Cannot rmdir $tmpDir\n";
    
    my $first_char_doublet_index = -1;
    my $last_char_doublet_index = -1;
    my $num_identities = 0;

    for (my $i = 0; $i < $#seqArray1; $i++) {
	my $char_1 = $seqArray1[$i];
	my $char_2 = $seqArray2[$i];
	if ($char_1 =~ /[a-z]/i && $char_2 =~ /[a-z]/i) {
	    ## chars at both alignment positions:
	    if ($char_1 eq $char_2) {
		$num_identities++;
	    }
	    
	    if ($first_char_doublet_index == -1) {
		$first_char_doublet_index = $i;
	    }

	    $last_char_doublet_index = $i;
	}
    }

    if ($first_char_doublet_index == -1) {
	die "Error, didn't read any alignment.  See file $tmpDir/$outFile\n";
    }

    my $alignLength = $last_char_doublet_index - $first_char_doublet_index + 1;
    
    my $percent_identity = $num_identities / $alignLength * 100;
    
    return ($alignLength, $percent_identity, [\@seqArray1, \@seqArray2]);
}


	    
