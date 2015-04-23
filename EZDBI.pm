package EZDBI;
use strict;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
use Exporter;
use DBI;
$VERSION = .01;
@ISA = qw(Exporter DBI CGI);
@EXPORT= qw(&connectToDb &doSimpleSQL &doSQL &printDebug);
@EXPORT_OK= qw(&connectToDb &doSimpleSQL &doSQL &printDebug);
%EXPORT_TAGS = (
		DB => [qw(&connectToDb &doSimpleSQL &doSQL &printDebug) ]
		);

##################
#Database access methods
##################
sub connectToDb {
    my($server,$dbtype,$user,$passwd,$db) = @_;
    my($dbproc);
    $dbproc = (DBI->connect("dbi:$dbtype:server=$server",$user, $passwd));
    if ( !defined $dbproc ) {
        die "Cannot connect to Sybase server: $server $DBI::errstr\n";
    }
    $dbproc->do("use $db");
    return($dbproc);
}
sub doSimpleSQL{
    my($dbproc,$query) = @_;
    my($statementHandle);
    $statementHandle = $dbproc->prepare($query);
    if ( !defined $statementHandle) {
	die "Cannot prepare statement: $DBI::errstr\n";
    }
    $statementHandle->execute() || die "failed query: $query\n";
    return $statementHandle->fetchrow();
}
sub doSQL{
    my($dbproc,$query) = @_;
    my($statementHandle,@results,$row);
    $statementHandle = $dbproc->prepare($query);
    if ( !defined $statementHandle) {
	die "Cannot prepare statement: $DBI::errstr\n";
    }
    $statementHandle->execute() || die "failed query: $query\n";
    if($statementHandle->{syb_more_results} NE "") {
        while ( $row = $statementHandle->fetchrow_hashref ) {
            push(@results,$row);
        }
    }
    return \@results;
}
##########################################
# DEBUGGING 
##########################################
sub printDebug{
    my($q) = new CGI;
    my($envkey);
    print $q->header();
    print "<body>";
    print $q->dump,"<pre>";
    print "BROWSER:user_agent=",$q->user_agent,"\n";
    print "PATH:path_info=",$q->path_info,"\n";
    print "PYS_PATH:path_translated=",$q->path_translated,"\n";
    print "HOST:remote_host=",$q->remote_host,"\n";
    print "SCRIPT:script_name",$q->script_name,"\n";
    print "SERVER:server_name=",$q->server_name,"\n";
    foreach $envkey (keys %ENV){
	print "*$envkey:",$ENV{$envkey},"\n";
    }
    print "</pre></body></html>\n";
}
1;






