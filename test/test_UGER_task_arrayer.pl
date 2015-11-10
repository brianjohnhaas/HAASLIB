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
    
    $uger_runner->run(@cmds);



    exit(0);
}


