package Asmbl_map_node;
use strict;
use lib ($ENV{EGC_SCRIPTS});
use Egc_library;

sub new {
    my $class = shift;
    my $dbproc = shift;
    my $asmbl_id = shift;
    my $map_type = shift;

    my $self = { asmbl_id=>$asmbl_id,
		 from_asmbl=>0,
		 to_asmbl=>0,
		 from_overlap=>0,
		 to_overlap=>0,
		 from_overlap_size=>0,
		 orientation=>0,
		 from_overhang_size=>0,
		 to_overhang_size=>0,
		 map_type=>$map_type};
    bless ($self, $class);

    $self->populate_data($dbproc, $asmbl_id, $map_type);

    return ($self);

}


sub populate_data {
    my $self = shift;
    my $dbproc = shift;
    my $asmbl_id = shift;
    my $map_type = shift;
    
    my $query = "select from_asmbl, to_asmbl, from_overlap, to_overlap, orientation, from_overhang_size, to_overhang_size, from_overlap_size from asmbl_map where asmbl_id = $asmbl_id and map_type = \"$map_type\"\n";
    my $result = &first_result_sql ($dbproc, $query);
    my ($from_asmbl, $to_asmbl, $from_overlap, $to_overlap, $orientation, $from_overhang_size, $to_overhang_size, $from_overlap_size) = split (/\t/, $result);
    $self->{from_asmbl} = $from_asmbl;
    $self->{to_asmbl} = $to_asmbl;
    $self->{from_overlap} = $from_overlap;
    $self->{to_overlap} = $to_overlap;
    $self->{orientation} = $orientation;
    $self->{from_overhang_size} = $from_overhang_size;
    $self->{to_overhang_size} = $to_overhang_size;
    $self->{from_overlap_size} = $from_overlap_size;
}

1; #EOM
    
