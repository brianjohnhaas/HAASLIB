package GenericGraph;

use strict;
use warnings;
use Carp;
use GenericNode;

sub new {
	my $packagename = shift;
	
	my $self = { 
		_nodes => {}, # node_name => node_reference
		_edge_counter => {},  # prev_node -> after_node = count
	};

	bless ($self, $packagename);

	return($self);
}


sub get_or_create_node {
	my $self = shift;
	my $node_name = shift;
	
	if ($self->node_exists($node_name)) {
		return($self->get_node($node_name));
	}
	else {
		# instantiate it, add it to the graph
		return($self->create_node($node_name));
	}
}
	

sub node_exists {
	my $self = shift;
	my $node_name = shift;
	
	unless ($node_name =~ /\w/) {
		confess "Error, node_name required";
	}
	
	if (exists $self->{_nodes}->{$node_name}) {
		return(1);
	}
	else {
		return(0);
	}
}

sub get_node {
	my $self = shift;
	my $node_name = shift;

	if (! $self->node_exists($node_name)) {
		confess "Error, $node_name doesn't exist in graph";
	}

	my $node = $self->{_nodes}->{$node_name};
	
	return($node);
}

sub get_all_nodes {
	my $self = shift;
	
	return(values %{$self->{_nodes}});
}


sub create_node {
	my $self = shift;

	my $node_name = shift;

	if ($self->node_exists($node_name)) {
		confess "Error, node $node_name already exists in the graph";
	}
    
	my $node = new GenericNode($node_name);
	
	$self->{_nodes}->{$node_name} = $node;
	
	return($node);
}


sub link_adjacent_nodes {
	my $self = shift;
	my ($before_node, $after_node, $edge_increment) = @_;

	unless (ref $before_node && ref $after_node) {
		confess "Error, need both before and after nodes for linking";
	}
	
	$before_node->add_next_node($after_node);
	$after_node->add_prev_node($before_node);

	if (defined($edge_increment)) {
		if ($edge_increment =~ /^\d+/ && $edge_increment > 0) {
			$self->{_edge_counter}->{$before_node}->{$after_node} += $edge_increment;
		}
	}
	else {
		$self->{_edge_counter}->{$before_node}->{$after_node}++;
	}
	
	
	return;
}


sub get_edge_count {
	my $self = shift;
	my ($prev_node, $next_node) = @_;

	my $edge_count = $self->{_edge_counter}->{$prev_node}->{$next_node} || 0;
	
	return($edge_count);
}


####
sub prune_nodes_from_graph {
	my $self = shift;
	my @nodes = @_;

	my $graph_nodes_href = $self->{_nodes};

	foreach my $node (@nodes) {

		delete ($self->{_edge_counter}->{$node}); # remove edge counts starting at current node.
		
		my $node_name = $node->get_value();
		delete $graph_nodes_href->{$node_name};
		
		my @next_nodes = $node->get_all_next_nodes();
		my @prev_nodes = $node->get_all_prev_nodes();
		
		foreach my $prev_node (@prev_nodes) {
			$prev_node->delete_next_node($node);
			$node->delete_prev_node($prev_node);
			
			delete ($self->{_edge_counter}->{$prev_node}->{$node}); # remove edge counts for prev nodes linking to current node.

		}
		
		foreach my $next_node (@next_nodes) {
			
			$next_node->delete_prev_node($node);
			$node->delete_next_node($next_node);
		}


	}
	
	
	return;
}

####
sub prune_edge {
	my $self = shift;
	my ($prev_node, $node) = @_;

	# Sever connection between nodes.

	$prev_node->delete_next_node($node);
	$node->delete_prev_node($prev_node);
	
	delete ($self->{_edge_counter}->{$prev_node}->{$node});

	return;
}


####
sub print_path {
	my @nodes = @_;

	my $counter = 0;
	foreach my $node (@nodes) {
		$counter++;
		printf("%4s", $counter);
		print " " . $node->get_value() . "\n";
	}
	print "\n";
	
	return;
}

####
sub toString {
	my $self = shift;

	my @nodes = $self->get_all_nodes();

	my $text = "";
	foreach my $node (@nodes) {
		$text .= $node->toString() . "\n";
	}

	return($text);
}

####
sub get_root_nodes {
	my $self = shift;

	my @roots;
	foreach my $node ($self->get_all_nodes()) {
		unless ($node->get_all_prev_nodes()) {
			push (@roots, $node);
		}
	}

	return(@roots);
}


####
sub toGraphViz {
    my $self = shift;
    my ($node_label_function_fref) = @_;
        
    my @nodes = $self->get_all_nodes();
    
    my $text = "digraph {\n";
    
    $text .= 
        #"node [width=0.1,height=0.1,fontsize=10,shape=point];\n"
        "node [width=0.1,height=0.1,fontsize=10];\n"
        . "edge [fontsize=12];\n"
        . "margin=1.0;\n"
        . "rankdir=LR;\n"
        . "labeljust=l;\n";
    
    
    my $node_counter = 0;
    
    foreach my $node (@nodes) {
    
        $node_counter++;
        #if ($node_counter > 1000) {
        #    last;
        #}
        
        my $id = hex($node->get_ID());
        
        my $node_label = "$id";
        if (ref $node_label_function_fref) {
            $node_label = &$node_label_function_fref($node);
        }
        
        $text .= "\t$id \[label=\"$node_label\"];\n";
        
        
        my @next_nodes = $node->get_all_next_nodes();
        
        foreach my $next_node (@next_nodes) {
                        
            my $next_id = hex($next_node->get_ID());
            
            my $edge_weight = $self->get_edge_count($node, $next_node);
            
            $text .= "\t$id->$next_id [label=$edge_weight];\n";
            
            
        }
    }
    
    $text .= "}\n";
        
    return($text);
}


####
sub to_sif_format {
    my $self = shift;
    
    my $ret_text = "";
    
    my  @nodes = $self->get_all_nodes();
    
    foreach my $node (@nodes) {
        
        my $id = hex($node->get_ID());
        
        my @next_nodes = $node->get_all_next_nodes();

        my @next_node_ids;
        foreach my $next_node (@next_nodes) {
            my $next_node_id = hex($next_node->get_ID());
            push (@next_node_ids, $next_node_id);
        }

        $ret_text .= join("\t", $id, "xx", @next_node_ids) . "\n";
    }


    return($ret_text);
}

1; #EOM
