use Egc_library;

####
sub delete_current_GO_assignments {
    my ($dbproc, $TU_feat_name) = @_;
    
    ## first delete the evidence:
    my $query = "delete go_evidence from go_role_link grl, go_evidence ge where grl.feat_name = \"$TU_feat_name\" and grl.id = ge.role_link_id";
    print "$query\n";
    
    &RunMod($dbproc, $query);


    my $query = "delete go_role_link where feat_name = \"$TU_feat_name\"";
    print "$query\n";
    &RunMod($dbproc, $query);

}


####
sub copy_GO_assignments {
    my ($dbproc, $from_TU, $destination_TU) = @_;
    
    my $query = "select gl.go_id, gl.id from go_role_link gl where gl.feat_name = \"$from_TU\" ";
    
    my @results = &do_sql ($dbproc, $query);
    foreach my $result (@results) {
        print "$result\n";
        my ($go_id, $gorolelink_id) = split (/\t/, $result);
        my $query = "insert go_role_link (feat_name, go_id, assigned_by, date) select \"$destination_TU\", \"$go_id\", assigned_by, date from go_role_link where id = $gorolelink_id\n";
        
        print "$query\n";
        &RunMod($dbproc, $query);
        #get the new row identifier
        my $query = "select id from go_role_link where feat_name = \"$destination_TU\" and go_id = \"$go_id\"\n";
        my $new_gorolelink_id = &first_result_sql ($dbproc, $query);
        unless ($new_gorolelink_id) {
            die "ERROR, couldn't retrieve data based on insert derived from $result.
\n";
        }
        my $query = "insert go_evidence (role_link_id, ev_code, evidence, with_ev) select $new_gorolelink_id, ev_code, evidence, with_ev from go_evidence where role_link_id = $gorolelink_id\n";
        print "$query\n";
        &RunMod($dbproc, $query);
        
    }
}

1; #require success.

