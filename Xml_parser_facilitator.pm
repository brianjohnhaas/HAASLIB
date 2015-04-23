#!/usr/local/bin/perl

package Xml_parser_facilitator;

my $DEBUG = 0;

####
sub tree_parser {
    my ($node_name, $node_ref, $handler_ref) = @_;
    print "TreeParser: $node_name\n" if $DEBUG;
    my ($node_attributes, @tag_content_pairs) = @{$node_ref};    
    if (exists($handler_ref->{$node_name})) {
	my $cur_index = $#{$handler_ref->{$node_name}};
	$cur_index++;
	$handler_ref->{$node_name}->[$cur_index] = $node_ref;
    } #else {
	while (@tag_content_pairs) {
	    my $next_node_name = shift (@tag_content_pairs);
	    my $next_node_ref = shift (@tag_content_pairs);
	    if ($next_node_name) {
		&tree_parser ($next_node_name, $next_node_ref, $handler_ref);
	    }
	}
   # }
}
	

####
sub branch_parser {
    my ($node_name, $node_ref, $obj_ref) = @_;
    print "BranchParser: $node_name\n" if $DEBUG;
    my $want_text = 0;
    my ($node_attributes, @tag_content_pairs) = @{$node_ref};
    print "content: @tag_content_pairs\n" if $DEBUG;
    #if node_name is specified, text and attributes will be assigned to node struct.
    if (exists($obj_ref->{$node_name})) {
	#get text
	$want_text = 1;
	#assign attributes
	foreach my $key (keys %{$node_attributes}) {
	    my $tempkey = $node_name . "_atts";
	    $obj_ref->{$tempkey}->{$key} = $node_attributes->{$key};
	}
    }
    while (@tag_content_pairs) {
	my $next_node_name = shift (@tag_content_pairs);
	my $next_node_ref = shift (@tag_content_pairs);
	if ($next_node_name) {
	    &branch_parser ($next_node_name, $next_node_ref, $obj_ref);
	} elsif ($want_text) { #grab text if avail.
	    if (ref ($obj_ref->{$node_name}) ne "ARRAY") {
		if ($next_node_ref =~ /\w+/) {
		    if ($obj_ref->{$node_name}) {
			#append to preexisting text
			$obj_ref->{$node_name} .= " / $next_node_ref";
		    } else {
			#assign text
			$obj_ref->{$node_name} = $next_node_ref; #assign text to element name
		    }
		}
	    } else {
		my $current_index = $#{$obj_ref->{$node_name}};
		$current_index++;
		$obj_ref->{$node_name}->[$current_index] = $next_node_ref;
	    }
	}
    }
    # could be looking for presence of an element only
    if ( (ref ($obj_ref->{$node_name}) ne "ARRAY") && (exists ($obj_ref->{$node_name})) && (!($obj_ref->{$node_name})) ) {
	$obj_ref->{$node_name} = 1;
    }
}


1;
