

=head1 NAME

package WuXdget

=cut



=head1 DESCRIPTION

    routines for extracting entries from Fasta file using the Wash-U's xdformat and xdget tools.

=cut

    ;

package main;
our $SEE;


package WuXdget;

use strict;
use warnings;
require Exporter;
use Carp;

our @ISA = qw (Exporter);
our @EXPORT = qw (xdget linearize xdget_linear);

## xdget and xdformat must be in path, otherwise the system will die.

=over 4

=item xdget()

B<Description:> Retrieves a fasta sequence entry from a fasta database

B<Parameters:> accession, fastaFilename, dbType

B<Returns:> fastaEntry

dbType = p || n

p for protein, n for nucleotide database type

use the linearize method to extract the fasta entry components

=back

=cut


    ;

sub xdget {
    my ($accession, $fastaFile, $type) = @_;

    if ($type eq 'p') {
        
        unless (-s "$fastaFile.xpi") {
            ## regenerate index file:
            my $cmd = "xdformat -p $fastaFile -I > /dev/null ";
            my $ret = system $cmd;
            if ($ret) {
                die "Error, couldn't create index file: $cmd, ret($ret)\n";
            }
        }
    }
    elsif ($type eq 'n') {
        unless (-s "$fastaFile.xni") {
            ## regenerate index file:
            my $cmd = "xdformat -n $fastaFile -I > /dev/null ";
            my $ret = system $cmd;
            if ($ret) {
                die "Error, couldn't create index file: $cmd, ret($ret)\n";
            }
        }
    }
    else {
        confess "Error, don't understand dbtype($type)\n";
    }
    
    my $cmd = "xdget -$type $fastaFile \'$accession\'";
    
    if ($SEE) {
        print "CMD: $cmd\n";
    }
    
    my $fastaEntry = `$cmd`;
    if ($?) {
        die "Error, couldn't run xdget: $cmd, ret($?)\n";
    }
    
    unless ($fastaEntry) {
        die "Error, no fasta entry retrieved by accession: $accession\n";
    }
    
    return ($fastaEntry);
}


=over 4

=item linearize()

B<Description:> breaks down a fasta sequence into its components 

B<Parameters:> fastaEntry

B<Returns:> (accession, header, linearSequence)

=back

=cut

    ;

sub linearize {
    my ($fastaEntry) = @_;
    
    unless ($fastaEntry =~ /^>/) {
        die "Error, fasta entry lacks expected format starting with header '>' character.\nHere's the entry\n$fastaEntry\n\n";
    }
    
    my @lines = split (/\n/, $fastaEntry);
    my $header = shift @lines;
    my $sequence = join ("", @lines);
    $sequence =~ s/\s+//g;
    
    $header =~ />(\S+)/;
    my $accession = $1;
    
    return ($accession, $header, $sequence);
}



=over 4

=item xdget_linear()

B<Description:> same as calling xdget (), and chasing it with linearize(), but only the sequence is returned

B<Parameters:> accession, fasta_db, dbtype

B<Returns:> linearSequence

see xdget() description

=back

=cut

    ;


sub xdget_linear {
    my ($acc, $fasta_db, $dbtype) = @_;
    
    my $fasta_entry = xdget($acc, $fasta_db, $dbtype);
    
    my ($acc2, $header, $genome_seq) = linearize($fasta_entry);

    return ($genome_seq);
}


1; #EOM
    
