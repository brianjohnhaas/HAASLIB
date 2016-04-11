#!/usr/bin/env perl

package UGER_task_arrayer;

use strict;
use warnings;
use Carp;
use File::Path;


sub new {
    my ($packagename, $logdir, $resume_mode) = @_;

    if ($resume_mode) {
        if (! -d $logdir) {
            confess "Error, in resume mode but logdir $logdir doesnt exist yet!";
        }
    }
    elsif (! -d $logdir) {
        mkpath($logdir);
    }
    
    
    my $self = { logdir => $logdir,
                 resume => $resume_mode || 0,
                 bash_header => "#!/bin/bash -l\n\nset -e\n\n",
                 queue => "short", # long|short
                 memory => "4g", # -l m_mem_free=${memory}g
                 threads => 1, # -pe smp 1
                 name => "noname",
                 project_name => "",  # set to 'regevlab' if in regevlab
    };
    
    bless ($self, $packagename);
    
    return($self);
}


####
sub set_bash_header {
    my ($self, $bash_header) = @_;

    $self->{bash_header} = $bash_header;
    
    return;
}

####
sub set_queue {
    my ($self, $queue_name) = @_;

    $self->{queue} = $queue_name;
    
    return;
}

####
sub set_memory {
    my ($self, $memory) = @_;
    
    unless ($memory =~ /^\d+g$/) {
        confess "Error, memory format must be NUMg (ie. 4g) ";
    }

    $self->{memory} = $memory;

    return;
}

####
sub set_threads {
    my ($self, $thread_count) = @_;

    $self->{threads} = $thread_count;

    return;
}

####
sub set_name {
    my ($self, $name) = @_;

    $self->{name} = $name;
    
    return;
}

####
sub set_project_name {
    my ($self, $project_name) = @_;

    $self->{project_name} = $project_name;

    return;
}


    
####
sub run {
    my $self = shift;
    my @cmds = @_;


    my $logdir = $self->{logdir};

    my $cmds_dir = "$logdir/cmds";
    if (! -d $cmds_dir) {
        mkpath($cmds_dir);
    }
    my $retvals_dir = "$logdir/ret";
    if (! -d $retvals_dir) {
        mkpath($retvals_dir);
    }

    my $runner_script = $self->_init_runner($logdir, $cmds_dir, $retvals_dir);


    ## Examine cache from previous run
    my $cache_success_file = "$logdir/cache_SUCCESS";
    my %cache_success;
    if (-s $cache_success_file) {
        open (my $fh, $cache_success_file) or die "Error, cannot open file $cache_success_file";
        while (<$fh>) {
            chomp;
            $cache_success{$_} = 1;
        }
        close $fh;
    }
    

    
    ##################################
    ## write individual command files:
    
    my $counter = 0;
    my @indices_to_track;
    foreach my $cmd (@cmds) {
        $counter++;

        if ($cache_success{$cmd}) {
            next; # already completed
        }
        
        push (@indices_to_track, $counter);
        
        my $cmd_file = "$cmds_dir/$counter.cmd";
        if (! -s $cmd_file) {
            # writing bash script for command
            open (my $ofh, ">$cmd_file") or die "Error, cannot write to $cmd_file";
            print $ofh $self->{bash_header};
            print $ofh "\n$cmd\n";
            print $ofh "exit \$?\n";
            close $ofh;

            chmod(0775, $cmd_file);
        }
    }
    
    my $num_cmds = scalar(@cmds);
    # launch the process and wait for results.
    
    my $project_info = "";
    if (my $project_name = $self->{project_name}) {
        $project_info = "-P $project_name";
    }

    my $qsub_cmd = "qsub -V -cwd -b y -sync y "
        . " -N " . $self->{name} . " $project_info" 
        . " -e $logdir/qsub.err -o $logdir/qsub.out "
        . " -q " . $self->{queue} . " " 
        . " -l m_mem_free=" . $self->{memory} . " "
        . " -pe smp " . $self->{threads} . " " 
        . " -t 1-$num_cmds $runner_script";
    

    print STDERR "CMD: $qsub_cmd\n";
    my $ret = system($qsub_cmd);

    ## audit successes/failures/unknowns.
    my $num_success = 0;
    my $num_error = 0;
    my $num_unknown = 0;

    open (my $cache_success_ofh, ">>$cache_success_file") or die "Error, cannot append to file $cache_success_file";
    my @failed_cmds;
    foreach my $index (@indices_to_track) {
        my $cmd = $cmds[$index - 1];
        my $retval_file = "$logdir/ret/$index.ret";
        if (-s $retval_file) {
            my $ret = `cat $retval_file`;
            chomp $ret;
            if ($ret == 0) {
                # success
                $num_success++;
                print $cache_success_ofh "$cmd\n";
            }
            else {
                push (@failed_cmds, $cmd);
                $num_error++;
            }
        }
        else {
            $num_unknown++;
            push (@failed_cmds, $cmd);
        }
    }
    if (@failed_cmds) {
        my $failed_cmds_file = "$logdir/FAILED_cmds.txt";
        open (my $ofh, ">$failed_cmds_file") or die "Error, cannot write to $failed_cmds_file";
        print $ofh join("\n", @failed_cmds) . "\n";
        close $ofh;
        
        print STDERR "WARNING: some commands failed.  Status = success($num_success), error($num_error), unknown($num_unknown)\n";
        
        return ($num_error + $num_unknown);
    
    }
    else {
        print STDERR "All $num_success commands succeeed. :) \n\n";
        return(0);
    }
    
}


####
sub _init_runner {
    my ($self, $logdir, $cmds_dir, $retvals_dir) = @_;


    my $runner_script = "$logdir/runner.pl";
    
    if (! -s $runner_script) {
        $self->_write_runner_script($runner_script, $logdir, $cmds_dir, $retvals_dir);
        
    }

    return($runner_script);

}


####
sub _write_runner_script {
    my ($self, $runner_script, $logdir, $cmds_dir, $retvals_dir) = @_;
    
    open (my $ofh, ">$runner_script") or die "Error, cannot write to $runner_script";
    print $ofh <<__EOSCRIPT;
#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
    
my \$SGE_TASK_ID = \$ENV{SGE_TASK_ID};
if (! defined \$SGE_TASK_ID) {
    confess "Error, env var SGE_TASK_ID not set";
}

    my \$script_name = "$logdir/cmds/\$SGE_TASK_ID.cmd";
    my \$retval_file = "$logdir/ret/\$SGE_TASK_ID.ret";

    my \$ret;
    if (-e \$retval_file) {
        \$ret = `cat \$retval_file`;
        chomp \$ret;
        if (\$ret == 0) {
            # no reason to rerun the command, already finished succesfully in earlier run.
            exit(0);
        }
    }
    
    \$ret = system(\$script_name);
    system("echo \$ret > \$retval_file");
    
    exit(\$ret);



__EOSCRIPT

    ;

    close $ofh;
    chmod (0775, $runner_script);

    return;
}



1; #EOM
