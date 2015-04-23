#!/usr/local/bin/perl


package main;
our $DEBUG;


## simple interface to the HTCRequest module
package Run_Condor;
# to use condor
#use lib ("/home/condor/lib");

# To use SGE
use lib qw(/home/sgeworker/lib);
use TIGR::HTCRequest;
use strict;
use Data::Dumper;

sub launch_condor_job {
    ## simplest thing to do is to write a perl script which does all the hard work, and present this perl script name as the cmd parameter here, in cases where parameter parsing is a problem.
    
    my ($workdir, $cmd, $parameter_list_aref, $inputFile_param_opt, $inputFiles_list_aref) = @_;
    
	
    my $request = TIGR::HTCRequest->new(group => "BHaas-EukAnnot",
                                        initialdir => $workdir,
                                        opsys=>"Linux");
    
    $request->set_command($cmd);
    
    foreach my $param (@$parameter_list_aref) {
        $request->add_param($param);
    }
    
    $request->add_param({key => "$inputFile_param_opt \$(Name)", value => $inputFiles_list_aref, type => "ARRAY"});
    
    $request->set_output("\$(Name).stdout");
    $request->set_error("\$(Name).stderr");
    $request->set_getenv(1);
    
    $request->length('long');
    
    if($DEBUG) {
        my $xx = $request->to_xml();
        print $xx;
    }
    
    my $id = $request->submit();
    print "Request id was $id \nDirectory: $workdir\n\n";
    
    
    $request->wait_for_request();
    open (my $log_fh, ">search.$$.log") or die $!;
    print $log_fh Dumper($request->get_tasks());
    close $log_fh;
    my $message = $request->get_message();
    if($request->get_state() eq "FAILURE") {
        print " Condor failed: $message \n";
        return (0);
    } else {
        print " Request finished with state " . $request->get_state() . " and $message \n";
        return(1);
    }
    
}



1; #EOM
