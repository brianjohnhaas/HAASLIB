#!/usr/local/bin/perl

use Mummer_coord_converter;

my $converter = new Mummer_coord_converter();
$converter->parse_Mummer_outputfiles('mum_forward', 'mum_reverse');
#$converter->data_dumper();
my $num;
unless ($num = $ARGV[0]) {
    $num = 1;
}

my $converted_num = $converter->transform_coordinate_1_to_2($num);
print "$num converted to $converted_num\n";
