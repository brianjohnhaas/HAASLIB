#!/usr/local/bin/perl

package main;
our $SEE;

=head1 NAME

CDNA::Btab_manager

=head1 DESCRIPTION

Module used to parse btab files generated resulting from a blast search.

=cut


package BLAST_btab_manager;
use strict;


=item new()

=over 4

B<Description:>Instantiate a Btab_manager object. 

Must instantiate a Btab_manager object and then use the parse_btab_file() method before using access methods.
There is a requirement that the btab file represents the search data matching a single query sequence.
B<Parameters:> none

B<Returns:> Btab_manager_obj



=back

=cut


sub new {
    my $packagename = shift;
    my $self = { query_acc_to_matches => {}, #hash keyed on accession points to list of btab chains.
		 db_acc_to_matches => {}
	     };
    bless ($self, $packagename);
}


=item get_accessions()

=over 4

B<Description:>Method retrieves and returns the database accessions found in the btab file. 

B<Parameters:> none

B<Returns:> @accessions

=back

=cut


sub get_query_accessions {
    my $self = shift;
    my @accessions = keys %{$self->{query_acc_to_matches}};
    return (@accessions);
}


sub get_db_accessions {
    my $self = shift;
    my @accessions =  keys %{$self->{db_acc_to_matches}};
    return (@accessions);
}




=item get_matches_via_acc()

=over 4

B<Description:> Returns all the matches that exist for a given database accession.

B<Parameters:> $accession

$accession is a string

B<Returns:> @Btab_lines

@Btab_chains is a list of Btab_line objects. 

=back

=cut

sub get_matches_via_query_acc {
    my ($self, $acc) = @_;
    my $list_ref = $self->{query_acc_to_matches}->{$acc};
    if (ref ($list_ref)) {
	return (@$list_ref);
    } else {
	return ();
    }
}



sub get_matches_via_db_acc {
    my ($self, $acc) = @_;
    my $list_ref = $self->{db_acc_to_matches}->{$acc};
    if (ref ($list_ref)) {
	return (@$list_ref);
    } else {
	return ();
    }
}










=item get_top_scoring_match_via_acc()

=over 4

B<Description:>Retrieves the highest scoring alignment chain for a given accession.  Score is the sum of the (per_id*length) for all alignment segments in a single alignment chain.

B<Parameters:>none

B<Returns:> Btab_chain_object

Btab_chain_object is of type CDNA::Btab_chain (see below).

=back

=cut




sub get_top_scoring_match_via_query_acc {
    my ($self, $acc) = @_;
    my @matches = $self->get_matches_via_query_acc($acc);
    @matches = sort {$a->{bit_score}<=>$b->{bit_score}} @matches;
    my $top_match = pop @matches;
    return ($top_match);
}


sub get_top_scoring_match_via_db_acc {
    my ($self, $acc) = @_;
    my @matches = $self->get_matches_via_db_acc($acc);
    @matches = sort {$a->{bit_score}<=>$b->{bit_score}} @matches;
    my $top_match = pop @matches;
    return ($top_match);
}


#private
sub add_btab_line {
    my $self = shift;
    my $query_acc = shift;
    my $db_acc = shift;
    my $btab_line = shift;
    ## Add to query list
    if (my $listref = $self->{query_acc_to_matches}->{$query_acc}) {
	push (@$listref, $btab_line);
    } else {
	$self->{query_acc_to_matches}->{$query_acc} = [$btab_line];
    }

    ## Add to db list
    if (my $listref = $self->{db_acc_to_matches}->{$db_acc}) {
	push (@$listref, $btab_line);
    } else {
	$self->{db_acc_to_matches}->{$db_acc} = [$btab_line];
    }
}

=item parse_btab_file()

=over 4

B<Description:> Method used to parse a btab file. This method should be used following the instantiation of a Btab_manager object.

B<Parameters:> btab_filename

btab_filename is the file name of the btab file.

B<Returns:> none.

=back

=cut

sub parse_btab_file {
    my ($self, $btabfilename) = @_;

    my @current_chain;
    my $current_key = ""; #join query_acc db_acc and chain_num

    open (BTAB, "$btabfilename") or die "Can't open $btabfilename";
    while (<BTAB>) {
	chomp;
	unless (/\w/) {next;}
	my $entry = CDNA::Btab_line->new($_);
	my $db_acc = $entry->{db_acc};
	my $query_acc = $entry->{query_acc};
	$self->add_btab_line($query_acc, $db_acc, $entry);
    }
    close BTAB;
}

=item dataDump()

=over 4

B<Description:>Method prints the btab file by calling the toString() method on each alignment chain. 

B<Parameters:> none.

B<Returns:> none.

=back

=cut

sub dataDump {
    my $self = shift;
    my $acc_to_matches = $self->{acc_to_matches};
    foreach my $db_acc (keys %$acc_to_matches) {
	print "\n\n// Entries for acc: $db_acc\n";
	my $match_list = $acc_to_matches->{$db_acc};
	foreach my $match (@$match_list) {
	    print $match->toString() . "\n";
	}
    }
}



####################################

package CDNA::Btab_line;
use strict;

=head1 NAME

CDNA::Btab_line

=head1 DESCRIPTION

Provides an object representing a single line within a btab file (ie. single alignment segment).


=item new()

=over 4

B<Description:> Instantiates a Btab_line object via a single line of a btab file (single alignment segment). 

B<Parameters:> btab_text_line

btab_text_line is a string containing a single line of a btab file.

B<Returns:> Btab_line_obj

Btab_line_obj is an object of type CDNA::Btab_line

The following fields are available:

B<query_acc>  (accession of the query sequence (genomic sequence))

B<date>  (date the search was run)

B<query_length>  (length of the query sequence (genomic sequence))

B<program>  (the name of the alignment program)

B<database>  (the Fasta-file database that was searched)

B<db_acc>  (accession of the database match (cDNA))

B<query_end5>  (the 5-prime coordinate of the query sequence match)

B<query_end3>  (the 3-prime coordinate of the query sequence match)

B<db_end5>  (the 5-prime coordinate of the database match)

B<db_end3>  (the 3-prime coordinate of the database match)

B<per_id>  (the percent identity of the match)

B<per_sim> (the percent similarity of the match) [optional]

B<bit_score> (alignment score, sometime not in bits...) [optional]

B<chain_num> (the alignment chain number in the btab file, uniquely identifying an alignment chain)

B<segment_num> (the alignment segment number within the chain)

B<db_header>  (the database header from the Fasta-file searched)

B<match_length> (the length of the alignment for a given alignment segment)

B<e_value> (the E-value (expect value) [optional])

B<p_value> (the P-value (probability value) [optional])

=back

=cut

sub new {
    my $packagename = shift;
    my $line = shift;
    my @x = split (/\t/, $line); #tab delimeted btab line
    my $self = {
	query_acc => $x[0],
	date => $x[1],
	query_length => $x[2],
	program => $x[3],
	database => $x[4],
	db_acc => $x[5],
	query_end5 => $x[6],
	query_end3 => $x[7],
	db_end5 => $x[8],
	db_end3 => $x[9],
	per_id => $x[10],
	per_sim => $x[11],
	bit_score => $x[12],
	db_header => $x[15],
	match_length => $x[18],
	e_value => $x[19],
	p_value => $x[20] };
    bless ($self, $packagename);
    return ($self);
}


=item toString()

=over 4

B<Description:> Returns a blurb of text describing the btab line. 

B<Parameters:> none

B<Returns:> text_line

=back

=cut


sub toString {
    my $self = shift;
    my $text = "query_acc: " . $self->{query_acc} . ", " 
	. "date: " . $self->{date} . ", "
	    . "query_length: " . $self->{query_length} . ", "
		. "program: " . $self->{program} . ", "
		    . "database: " . $self->{database} . ", "
			. "db_acc: " . $self->{db_acc} . ", "
			    . "query_end5: " . $self->{query_end5} . ", "
				. "query_end3: " . $self->{query_end3} . ", "
				    . "db_end5: " . $self->{db_end5} . ", "
					. "db_end3: " . $self->{db_end3} . ", "
					    . "per_id: " . $self->{per_id};
    if (my $x = $self->{per_sim}) {
	$text .= ", per_sim: " . $x;
    }
    if (my $x = $self->{bit_score}) {
	$text .= ", bit_score: " . $x;
    }
    
    if (my $x = $self->{db_header}) {
	$text .= ", header: " . $x;
    }
    if (my $x = $self->{match_length}) {
	$text .= ", match_length: " . $x;
    }
    if (my $x = $self->{e_value}) {
	$text .= ", e-value: " . $x;
    } 
    if (my $x = $self->{p_value}) {
	$text .= ", p-value: " . $x;
    }
    return ($text);
}

1; #EOM


