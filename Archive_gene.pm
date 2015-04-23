package Archive_gene;

use strict;
use lib ($ENV{EUK_MODULES}, $ENV{EGC_SCRIPTS});
use Egc_library;
use DBmodel_to_geneobj;

sub archive {
    my ($dbproc, $TU, $comment, $operation) = @_;
    my $gene_xml = DBmodel_to_geneobj::TU_to_XML($dbproc, $TU);
    $gene_xml =~ s/\'//g; #remove hard quotes.
    $gene_xml =~ s/\n|\t+/ /g; #remove tabs and newlines.
    ## get rest info 
    my $query = "select locus, pub_locus, alt_locus, asmbl_id, feat_type from ident i, asm_feature a where a.feat_name = i.feat_name and i.feat_name = \"$TU\"\n";
    my $identdata = &first_result_sql ($dbproc, $query);
    my ($locus, $pub_locus, $alt_locus, $asmbl_id, $feat_type) = split (/\t/, $identdata);
    ## Can only reference one locus id in feature_history.  Choose appropriately.
    if ($pub_locus) {
        $locus = $pub_locus;
    } elsif ($alt_locus && !$locus) {
        $locus = $alt_locus;
    }

    unless ($feat_type) {
        $feat_type = "na";
    }
    
    my $s = { TU_feat_name => $TU,
              locus =>$locus,
              feat_type => $feat_type,
              asmbl_id => $asmbl_id,
              operation => $operation,
              gene_xml => $gene_xml,
              assignby => "script",
              comment => $comment,
          };

    eval {
        &archive_gene_xml ($dbproc, $s);
    };

    if ($@) {
        print "Error trying to archive gene $TU\n";
    }
}


1; #EOM
