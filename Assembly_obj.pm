#!/usr/local/bin/perl

package Gene_obj;

use strict;

sub new {
    shift;
    my $self = { accession => 0,     #genbank accession
		 seq_group => 0,     #seq_group
		 tigr_asmbl_id => 0, #asmbl_id annotation db
		 clone_id => 0,      #clone_id (tracking info)
		 gb_desc => 0,       #genbank description text
		 gb_comment => 0,    #genbank comment text
		 sequence => 0       #nucleotide sequence
	     };
    bless($self);
    return ($self);
}

1;


