#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../");
use UGER_task_arrayer;
use Cwd;


main: {

    my $logdir = cwd() . "/uger_log_$$";
    
    my $uger_runner = new UGER_task_arrayer($logdir);
    
    my @cmds;
    foreach my $i (1..10) {
        
        push (@cmds, "echo hello world, task: $i");
        
    }
    

    my $fail_cmd = "no_such_command this_should_fail";

    print STDERR "** Initial submission of 11 jobs, 1 should fail\n";
    $uger_runner->run(@cmds, $fail_cmd);

    
    print STDERR "** And submitting a second round... should just try the failed command this time.\n";
    $uger_runner->run(@cmds, $fail_cmd);


    print STDERR "** Now running with just the commands that succeded earlier, plus one more for fun.\n";
    $uger_runner->run(@cmds, "echo one more for fun");
    
    print STDERR "\n\tDone testing.\n\n";
    
    exit(0);
}


