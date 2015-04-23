package main;
our ($SEE, $DEBUG);

package Sqlite_connect;

require 5.6.0;
require Exporter;
use Carp;
use strict;
use Data::Dumper;
use DBI;

our @ISA = qw(Exporter);

## export the following for general use
our @EXPORT = qw ($QUERYFAIL do_sql_2D connect_to_db RunMod first_result_sql very_first_result_sql do_sql);

our $QUERYFAIL = 0; #intialize.  Status flag, indicating the success of a query.

############### DATABASE CONNECTIVITY ################################
####
sub connect_to_db {
    my ($db_name) = @_;
    
    my $dbproc = DBI->connect("DBI:SQLite:$db_name");
    unless (ref $dbproc) {
        croak "Cannot connect to $db_name: $DBI::errstr";
    }
    $dbproc->{RaiseError} = 1; #turn on raise error.  Must use exception handling now.
    return($dbproc);
}


sub do_sql {
    my ($dbproc, $query) = @_;

    my @results_2D = &do_sql_2D($dbproc, $query);
   

    my @results;
    foreach my $result (@results_2D) {

        my $line = join("\t", @$result);
        push (@results, $line);
    }


    return(@results);
}



## return results in 2-Dimensional array.
sub do_sql_2D {
    my ($dbproc,$query, @values) = @_;
    my ($statementHandle,@x,@results);
    my ($i,$result,@row);
    
    ## Use $QUERYFAIL Global variable to detect query failures.
    $QUERYFAIL = 0; #initialize
    print "QUERY: $query\tVALUES: @values\n" if($::DEBUG||$::SEE);
    $statementHandle = $dbproc->prepare($query);
    if ( !defined $statementHandle) {
        print "Cannot prepare statement: $DBI::errstr\n";
        $QUERYFAIL = 1;
    } else {
        
        # Keep trying to query thru deadlocks:
        do {
            $QUERYFAIL = 0; #initialize
            eval {
                $statementHandle->execute(@values);
                while ( @row = $statementHandle->fetchrow_array() ) {
                    push(@results,[@row]);
                }
            };
            ## exception handling code:
            if ($@) {
                print STDERR "failed query: <$query>\tvalues: @values\nErrors: $DBI::errstr\n";
                $QUERYFAIL = 1;
            }
            
        } while ($statementHandle->errstr() =~ /deadlock/);
        #release the statement handle resources
        $statementHandle->finish;
    }
    if ($QUERYFAIL) {
        die "Failed query: <$query>\tvalues: @values\nErrors: $DBI::errstr\n";
    }
    return(@results);
}

sub RunMod {
    my ($dbproc,$query, @values) = @_;
    my ($result);
    if($::DEBUG||$::SEE) {print "QUERY: $query\tVALUES: @values\n";}
    if($::DEBUG) {
        $result = "NOT READY";
    } else {
        eval {
            $dbproc->do($query, undef, @values);
        };
        if ($@) { #error occurred
            die "failed query: <$query>\tvalues: @values\nErrors: $DBI::errstr\n";
            
        }
    }
}


sub first_result_sql {
    my ($dbproc, $query, @values) = @_;
    my @results = &do_sql_2D ($dbproc, $query, @values);
    return ($results[0]);
}

sub very_first_result_sql {
    my ($dbproc, $query, @values) = @_;
    my @results = &do_sql_2D ($dbproc, $query, @values);
    if ($results[0]) {
        return ($results[0]->[0]);
    } else {
        return (undef());
    }
}



1; #EOM
