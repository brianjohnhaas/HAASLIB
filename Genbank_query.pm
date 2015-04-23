package main;

our $SEE;

package Genbank_query;
use strict;
require Exporter;
our @ISA = qw (Exporter);
our @EXPORT = qw (ESearch EPost ESummary EFetch ELink);



=head1 NAME

package Genbank_query

=cut

=head1 DESCRIPTION

This package includes the E-utilities provided by Genbank.

http://eutils.ncbi.nlm.nih.gov/entrez/query/static/eutils_help.html

 
=over 4

=item Utility Listing
   
*B<ESearch:>  Searches and retrieves primary IDs (for use in EFetch, ELink and ESummary) and term translations, and optionally retains results for future use in the user\'s environment.
   
*B<EPost>:  Posts a file containing a list of primary IDs for future use in the user\'s environment to use with subsequent search strategies.
       
*B<ESummary:> Retrieves document summaries from a list of primary IDs or from the user\'s environment.
   
*B<EFetch:>  Retrieves records in the requested format from a list of one or more primary IDs or from the user\'s environment.
   
*B<ELink:>  Checks for the existence of an external or Related Articles link from a list of one or more primary IDs.  Retrieves primary IDs and relevancy scores for links to Entrez databases or Related Articles;  creates a hyperlink to the primary LinkOut provider for a specific ID and database, or lists LinkOut URLs and Attributes for multiple IDs.
   
=back


User System Requirements
Do not overload NCBI\'s systems. Users intending to send numerous queries and/or retrieve large numbers of records from Entrez should comply with the following:

    * Run retrieval scripts on weekends or between 9 PM and 5 AM ET weekdays for any series of more than 100 requests.
    * Make no more than one request every 3 seconds.
    * NCBI\'s Disclaimer and Copyright notice must be evident to users of your service.  NLM does not claim the copyright on the abstracts in PubMed; however, journal publishers or authors may. NLM provides no legal advice concerning distribution of copyrighted materials, consult your legal counsel.

Database primary IDs:
	
Genome          Genome ID
Nucleotide      GI number
OMIM 	        MIM number
PopSet 	        GI number
Protein 	GI number
PubMed 	        PMID
Structure 	MMDB ID
Taxonomy 	TAXID




=cut


    ;

######################################################################################################################    


=head1 DESCRIPTION

ESearch
Last updated: July 21, 2003

ESearch:  Searches and retrieves primary IDs (for use in EFetch, ELink and ESummary) and term translations, and optionally retains results for future use in the user\'s environment.

    * URL Parameters
          o Database
          o History  Web Environment  Query_key Tool  E-mailAddress
                + PubMed
                      # Search Terms  Search Field  Relative Dates  Date Ranges  Date Type  Display Numbers  Retrieval Mode  Retrieval Type   Examples  Sample XML Retrieval
                +  Journals
    * User System Requirements
    * Entrez DTDs
    * Demonstration Program
    * Announcement Mailing List
    * Help Desk

URL parameters:

Utility parameters may be case sensitive, therefore, use lower case characters in all parameters except for WebEnv.

Base URL: http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?

Database:

    db=database name

Current database values by category:

    Literature databases:

        omim - The OMIM database including the collective data is the property of the Johns Hopkins University, which holds the copyright.
        pubmed - Journal publishers hold the copyright on the abstracts in PubMed. NLM provides no legal advice concerning distribution of copyrighted materials.
        journals

    Sequence databases:

        genome
        nucleotide
        protein
        popset 

    3D database:

        structure

    Taxonomy database:

        taxonomy

History: Requests utility to maintain results in user\'s environment. Used in conjunction with WebEnv.

    usehistory=y

Web Environment: Value previously returned in XML results from ESearch or EPost. This value may change with each utility call. If WebEnv is used, History search numbers can be included in an ESummary URL, e.g., term=cancer+AND+%23X (where %23 replaces # and X is the History search number).

    WebEnv=WgHmIcDG]B\`&gt;&gt; etc.

Query_key:  The value used for a history search number or previously returned in XML results from ESearch or EPost.

    query_key=6

Note: WebEnv is similar to the cookie that is set on a user\'s computers when accessing PubMed on the web.  If the parameter usehistory=y is included in an ESearch URL both a WebEnv (cookie string) and query_key (history number) values will be returned in the results. Rather then using the retrieved PMIDs in an ESummary or EFetch URL you may simply use the WebEnv and query_key values to retrieve the records. WebEnv will change for each ESearch query, but a sample URL would be as follows:

http://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed
&WebEnv=%3D%5DzU%5D%3FIJIj%3CC%5E%5DA%3CT%5DEACgdn%3DF%5E%3Eh
GFA%5D%3CIFKGCbQkA%5E_hDFiFd%5C%3D
&query_key=6&retmode=html&rettype=medline&retmax=15

Tool: A string with no internal spaces that identifies the resource which is using Entrez links (e.g., tool=flybase). This argument is used to help NCBI provide better service to third parties generating Entrez queries from programs. As with any query system, it is sometimes possible to ask the same question different ways, with different effects on performance. NCBI requests that developers sending batch requests include a constant \'tool\' argument for all requests using the utilities.

    tool=

E-mail Address: If you choose to provide an email address, we will use it to contact you if there are problems with your queries or if we are changing software interfaces that might specifically affect your requests. If you choose not to include an email address we cannot provide specific help to you,  but you can still sign up for utilities-announce to receive general announcements.

    email=

PubMed:

Search terms: This command uses search terms or phrases with or without Boolean operators.  See the PubMed Help for information about search term qualification.

    term=search strategy

For example:

            term=asthma[mh]+OR+hay+fever[mh]

You may also qualify search terms using field=qualifier.

Search Field: Use this command to specify a specific search field.

    field=

PubMed fields: affl, auth, ecno, jour, iss, mesh, majr, mhda, page, pdat, ptyp, si, subs, subh, tiab, word, titl, lang, uid, fltr, vol

Relative Dates: Limit items a number of days immediately preceding today\'s date.

    reldate=

For example:

      reldate=90
      reldate=365

Date Ranges:  Limit results bounded by two specific dates. Both mindate and maxdate are required if date range limits are applied using these variables.

    mindate=
    maxdate=

For example:

      mindate=2001
      maxdate=2002/01/01

Date Type:  Limit dates to a specific date field based on database.

    datetype=

For example:

    datetype=edat 

Display Numbers:

    retstart=x  (x= sequential number of the first record retrieved - default=0 which will retrieve the first record)
    retmax=y  (y= number of items retrieved)

Retrieval Mode:

    retmode=xml

Use your web browser\'s View Page Source function to display results.

Retrieval Type:

    rettype=

PubMed values:

    count
    uilist (default)

Examples:
Search in PubMed for the term cancer for the entrez date from the last 60 days and retrieve the first 100 IDs and translations using the history parameter:
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=cancer&reldate=60&datetype=edat&retmax=100&usehistory=y

Search in PubMed for the journal PNAS Volume 97, and retrieve 6 IDs starting at ID 7 using a tool parameter:
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=PNAS[ta]+AND+97[vi]&retstart=6&retmax=6&tool=biomed3

Search in Journals for the term obstetrics:
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=journals&term=obstetrics

=cut


    ;




=over 4

=item ESearch()

B<Description:> Search specified GenBank Database using search parameters.

B<Parameters:> ($UserAgent, $URLvals_href)

$UserAgent is a LWP::UserAgent object

$URLvals_href is a reference to a hash containing values for URL keys.  Available keys include:

    db
    usehistory
    WebEnv
    query_key
    tool
    email
    term
    field
    reldate
    mindate
    maxdate
    datetype
    retstart
    retmax
    retmode
    rettype
    
see details above for expected contents.


B<Returns:> $text

$text includes the unparsed data returned from GenBank based on the search request. 

=back

=cut


sub ESearch {
    my ($ua, $URLvals_href) = @_;
    
    my $base_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?";

    my @valid_keys = qw(db
			usehistory
			WebEnv
			query_key
			tool
			email
			term
			field
			reldate
			mindate
			maxdate
			datetype
			retstart
			retmax
			retmode
			rettype);

    return (&process_request($ua, $base_url, \@valid_keys, $URLvals_href));
    
}


#######################################################################################################

=head1 DESCRIPTION

    EPost
Updated: July 21, 2003

EPost:  Posts a file containing a list of UIs for future use in the user\'s environment to use with subsequent search strategies.

    * URL Parameters
          o Database
          o Record Identifier  Retrieval Mode  Web Environment Query_key  Tool E-mail Address
                + PubMed
                      # Example
                + Protein, Nucleotide, Structure, Genome, PopSet, OMIM, Taxonomy, Books, ProbeSet, 3D Domains, UniSTS, Domains, SNP, Journals, UniGene
    * User System Requirements
    * Entrez DTDs
    * Demonstration Program
    * Announcement Mailing List
    * Help Desk

URL parameters:

Utility parameters may be case sensitive, therefore, use lower case characters in all parameters except for WebEnv.

Base URL: http://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi?

Database:

    db=database name

Current database values by category:

    Literature databases:

        omim - The OMIM database including the collective data is the property of the Johns Hopkins University, which holds the copyright.
        pubmed -  Journal publishers hold the copyright on the abstracts in PubMed. NLM provides no legal advice concerning distribution of copyrighted materials.

    Sequence databases:

        genome
        nucleotide
        protein
        popset
        *sequences- Composite name including nucleotide, protein, popset and genome.

    3D database:

        structure

    Taxonomy database:

        taxonomy

*Not yet available

Record Identifier: UIs required if web environment (i.e., WebEnv=) is not used.

    id=11877539,11822933,11871444

Current values:

    PubMed ID
    MEDLINE UI
    GI number
    MMDB ID
    TaxID
    MIM number

Retrieval Mode:

    retmode=xml (default)

Note: Use your web browser\'s View Page Source function to display results.

Web Environment: Value previously returned in XML results from ESearch. Web environment is required in place of a primary ID result list.

    WebEnv=WgHmIcDG]B\`&gt;&gt; etc.

Query_key: The value used for a history search number or previously returned in XML results from ESearch or EPost.

    query_key=6 

Note: WebEnv is similar to the cookie that is set on a user\'s computers when accessing PubMed on the web.  If the parameter usehistory=y is included in an ESearch URL both a WebEnv (cookie string) and query_key (history number) values will be returned in the results. Rather then using the retrieved PMIDs in an EPost URL you may simply use the WebEnv and query_key values to retrieve the records. WebEnv will change for each ESearch query.

Tool: A string with no internal spaces that identifies the resource which is using Entrez links (e.g., tool=flybase). This argument is used to help NCBI provide better service to third parties generating Entrez queries from programs. As with any query system, it is sometimes possible to ask the same question different ways, with different effects on performance. NCBI requests that developers sending batch requests include a constant \'tool\' argument for all requests using the utilities.

    tool=

E-mail Address:  If you choose to provide an email address, we will use it to contact you if there are problems with your queries or if we are changing software interfaces that might specifically affect your requests. If you choose not to include an email address we cannot provide specific help to you,  but you can still sign up for utilities-announce to receive general announcements.

    email=

PubMed Example:
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi?db=pubmed&id=11237011

=cut




=over 4

=item EPost()

B<Description:> Posts a file containing a list of UIs for future use in the user\'s environment to use with subsequent search strategies.

B<Parameters:> ($UserAgent, $URLvals_href)


$UserAgent is a LWP::UserAgent object

$URLvals_href is a reference to a hash containing values for URL keys.  Available keys include:

db
id
retmode
WebEnv
query_key
tool
email


    
see details above for expected contents.


B<Returns:> $text

$text includes the unparsed data returned from GenBank based on the post request. 

=back

=cut


####
sub EPost {
    my ($ua, $URLvals_href) = @_;
    
    my $base_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi?";

    my @valid_keys = qw(db
			id
			retmode
			WebEnv
			query_key
			tool
			email);
			

    return (&process_request($ua, $base_url, \@valid_keys, $URLvals_href));
}



    ;

#################################################################################################################

=head1 DESCRIPTION

    ESummary
    Last update: July 21, 2003

ESummary:  Retreives document Summaries from a list of primary IDs or from the user's environment.

    * URL Parameters
          o Database
          o History  Web Environment Query_key  Tool E-mail Address
                + PubMed
                      # Record Identifier  Display Numbers  Retrieval Mode  Sample XML Retrieval Examples
                + Journals
                + Sequence Databases Examples
    * User System Requirements
    * Entrez Database ESummary Fields
    * Entrez DTDs
    * Demonstration Program
    * Announcement Mailing List
    * Help Desk

URL parameters:

Utility parameters may be case sensitive, therefore, use lower case characters in all parameters except for WebEnv.

Base URL: http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?

Database:

    db=database name

Current database values: pubmed, protein, nucleotide, structure, genome, pmc, omim, taxonomy, books, probeset,  domains, unists, cdd, snp, journals, unigene, popset

omim - The OMIM database including the collective data is the property of the Johns Hopkins University, which holds the copyright.
pubmed - Journal publishers hold the copyright on the abstracts in PubMed. NLM provides no legal advice.
sequence
History: Requests utility to maintain results in History server, used in conjunction with WebEnv.

    usehistory=y

Web Environment: Value previously returned in XML results from ESearch and EPost and used with ESummary in place of primary ID result list.

    WebEnv=WgHmIcDG]B`&gt;&gt; etc.

Query_key: The value used for a history search number or previously returned in XML results from ESearch or EPost.

    query_key=6 

Note: WebEnv is similar to the cookie that is set on a user's computers when accessing PubMed on the web.  If the parameter usehistory=y is included in an ESearch URL both a WebEnv (cookie string) and query_key (history number) values will be returned in the results. Rather then using the retrieved PMIDs in an ESummary URL you may simply use the WebEnv and query_key values to retrieve the records. WebEnv will change for each ESearch query.

Tool: A string with no internal spaces that identifies the resource which is using Entrez links (e.g. tool=igm or tool=flybase). This argument is used to help NCBI provide better service to third parties generating Entrez queries from programs. As with any query system, it is sometimes possible to ask the same question different ways, with different effects on performance. NCBI requests that developers sending batch requests include a constant 'tool' argument for all requests using the utilities.

    tool=

E-mail Address: If you choose to provide an email address, we will use it to contact you if there are problems with your queries or if we are changing software interfaces that might specifically affect your requests. If you choose not to include an email address we cannot provide specific help to you,  but you can still sign up for utilities-announce to receive general announcements.

    email=

PubMed

Record Identifier: Required if WebEnv is not used.

    id=12345,92932

Current values:

    PubMed ID
    MEDLINE UI

Display Numbers: Used when the results from EPost or ESearch are maintained in the user's environment. The maximum number of retrieved records is 10,000.

    retstart=x  (x= sequential number of the first record retrieved - default=0 which will retrieve the first record)
    retmax=y  (y= number of records retrieved - default=20)

Retrieval Mode:

    retmode=xml

Note: Use your web browser's View Page Source function to display results.

Example:
In PubMed display records for PMIDs 11850928 and 11482001 in xml retrieval mode:
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=11850928,11482001&retmode=xml

In Journals display records for journal IDs 27731,439,735,905:
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=journals&id=27731,439,735,905


Sequence Databases

Record Identifier: Required if WebEnv is not used.

    id=28864546,28800981

Current values:

    GI number
    MMDB ID (Structure database)
    TAX ID (Taxonomy database)

Example:
In Protein display records for GIs 28800982 and 28628843 in xml retrieval mode:
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=protein&id=28800982,28628843&retmode=xml

In Nucleotide display records for GIs 28864546 and 28800981 in xml retrieval mode:
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nucleotide&id=28864546,28800981&retmode=xml

In Structure display records for MMDB IDs 19923 and 12120 in xml retrieval mode:
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=structure&id=19923,12120&retmode=xml

In Taxonomy display records for TAXIDs 9913 and 30521 in xml retrieval mode:
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy&id=9913,30521&retmode=xml

In UniSTS display records for IDs 254085 and 254086 in xml retrieval mode:
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=unists&id=254085,254086&retmode=xml

=cut



=over 4

=item ESummary()

B<Description:> Posts a file containing a list of UIs for future use in the user\'s environment to use with subsequent search strategies.

B<Parameters:> ($UserAgent, $URLvals_href)


$UserAgent is a LWP::UserAgent object

$URLvals_href is a reference to a hash containing values for URL keys.  Available keys include:

db
usehistory
WebEnv
query_key
tool
email
id
retstart
retmax
retmode
    
see details above for expected contents.


B<Returns:> $text

$text includes the unparsed data returned from GenBank based on the summary request. 

=back

=cut



sub ESummary {

    my ($ua, $URLvals_href) = @_;
    
    my $base_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?";
    
    my @valid_keys = qw(db
			usehistory
			WebEnv
			query_key
			tool
			email
			id
			retstart
			retmax
			retmode
			);
    
    return (&process_request($ua, $base_url, \@valid_keys, $URLvals_href));
}



#############################################################################################################

=head1 DESCRIPTION

    EFetch Overview
Last updated: July 21, 2003

EFetch:  Retrieves records in the requested format from a list of one or more UIs or from user\'s environment. Click on a database below to display database specific documentation.

    * URL Parameters
          o Database
          o Web Environment Query_key  Tool E-mail Address
                + PubMed  Journals
                + Protein, Nucleotide, Taxonomy, Genome, PopSet
    * User System Requirements
    * Entrez DTDs
    * Demonstration Program
    * Announcement Mailing List
    * Help Desk

URL parameters:

Utility parameters may be case sensitive, therefore, use lower case characters in all parameters except for WebEnv.

Base URL: http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?

EFetch for the Sequence Databases
Last Updated: July 21, 2003

EFetch documenation is also available for the Literature, and Taxonomy databases.

EFetch:  Retrieves records in the requested format from a list of one or more unique identifiers.

    * URL Parameters
    * User System Requirements
    * Entrez DTDs
    * Demonstration Program
    * Announcement Mailing List

URL parameters:

Utility parameters may be case sensitive, therefore, use lower case characters in all parameters except for WebEnv.

Base URL: http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?

Database

    db=nucleotide

Current database values by category:

    Sequence databases:

        genome
        nucleotide
        protein
        popset
        sequences - Composite name including nucleotide, protein, popset and genome.

Web Environment: History link value previously returned in XML results from ESearch and used with EFetch in place of primary ID result list.

    WebEnv=WgHmIcDG], etc.

Query_key: The value used for a history search number or previously returned in XML results from Esearch or EPost.

    query_key=6 

Note: WebEnv is similar to the cookie that is set on a user's computers when accessing PubMed on the web.  If the parameter usehistory=y is included in an ESearch URL both a WebEnv (cookie string) and query_key (history number) values will be returned in the results. Rather then using the retrieved PMIDs in an ESummary URL you may simply use the WebEnv and query_key values to retrieve the records. WebEnv will change for each ESearch query.

Tool: A string with no internal spaces that identifies the resource which is using Entrez links (e.g., tool=flybase). This argument is used to help NCBI provide better service to third parties generating Entrez queries from programs. As with any query system, it is sometimes possible to ask the same question different ways, with different effects on performance. NCBI requests that developers sending batch requests include a constant 'tool' argument for all requests using the utilities.

    tool=

E-mail Address: If you choose to provide an email address, we will use it to contact you if there are problems with your queries or if we are changing software interfaces that might specifically affect your requests. If you choose not to include an email address we cannot provide specific help to you,  but you can still sign up for utilities-announce to receive general announcements.

    email=



B<Sequence Databases>

Record Identifier: IDs required if WebEnv is not used.

    id=123,U12345,U12345.1,gb|U12345|

Current values:

    NCBI sequence number (GI)
    genome ID
    accession
    accession.version
    fasta
    seqid

Display Numbers:

    retstart=x  (x= sequential number of the first id retrieved - default=0 which will retrieve the first record)
    retmax=y  (y= number of items retrieved)

Sequence Strand, Start, Stop and Complexity Parameters
strand= 	what strand of DNA to show (1=plus or 2=minus)
seq_start= 	show sequence starting from this base number
seq_stop= 	show sequence ending on this base number
complexity= 	gi is often a part of a biological blob, containing other gis

complexity regulates the display:

    0 - get the whole blob
    1 - get the bioseq for gi of interest (default in Entrez)
    2 - get the minimal bioseq-set containing the gi of interest
    3 - get the minimal nuc-prot containing the gi of interest
    4 - get the minimal pub-set containing the gi of interest

Retrieval Mode:

    retmode=output format

Current values:

    xml
    html
    text
    asn.1

Retrieval Type:

    rettype=output types based on database

Current values:

    native (full record)
    fasta
    gb
    gbwithparts
    est
    gss
    gp
    uilist 

Type descriptions:
native 	Default format for viewing sequences
fasta 	FASTA view of a sequence
gb 	GenBank view for sequences, constructed sequences will be shown as contigs (by pointing to its parts).  Valid for nucleotides.
gbwithparts 	GenBank view for sequences, the sequence will always be shown. Valid for nucleotides
est 	EST Report. Valid for sequences from dbEST database.
gss 	GSS Report. Valid for sequences from dbGSS database.
gp 	GenPept view. Valid for proteins.
seqid 	To convert list of gis into list of seqids
acc 	To convert list of gis into list of accessions 

Not all Retrieval Modes are possible with all Retrieval Types.

Sequence Options:
 
	native 	fasta 	gb 	gbwithparts 	est 	gss 	gp 	seqid 	acc
xml 	x 	x* 	x* 	TBI 	        TBI 	TBI 	x* 	TBI 	TBI
text 	x 	x 	x* 	x* 	        x* 	x* 	x* 	x 	x
html 	x 	x 	x* 	x* 	        x* 	x* 	x* 	x 	x
asn.1 	x 	n/a 	n/a 	n/a 	        n/a 	n/a 	n/a 	x 	n/a

x = retrieval mode available
*  - existence of the mode depends on gi type
TBI - to be implemented (not yet available)
n/a - not available
 

Examples:

http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=5

http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=5&complexity=0&rettype=fasta

http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=5&rettype=gb&seq_start=1&seq_stop=9

http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.cgi?db=nucleotide&id=5&rettype=fasta&seq_start=1&seq_stop=9&strand=2

http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=5&rettype=gb

http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=popset&id=12829836&rettype=gp

http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=8&rettype=gp

Entrez display format GBSeqXML:
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=5&rettype=gb&retmode=xml
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=8&rettype=gp&retmode=xml

Entrez display format TinySeqXML:
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=5&rettype=fasta&retmode=xml 



B<EFetch for Literature Databases>
Last updated: July 21, 2003

EFetch:  Retrieves records in the requested format from a list of one or more UIs or the user's environment.

EFetch documentation is also available for the Sequence, and Taxonomy databases.

    * URL Parameters
          o Literature Databases
          o Web Environment Query_key  Tool E-mail Address
                + PubMed
                      # Record Identifier  Display Numbers  Retrieval Mode Retrieval Type PubMed Retrieval Options Examples
                + Journals 
    * User System Requirements
    * Entrez DTDs
    * Demonstration Program
    * Announcement Mailing List
    * Help Desk

URL parameters:

Utility parameters may be case sensitive, therefore, use lower case characters in all parameters except for WebEnv.

Base URL: http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?

Literature Database

    db=pubmed

        pubmed - Journal publishers hold the copyright on the abstracts in PubMed. NLM provides no legal advice concerning distribution of copyrighted materials.
        journals

Web Environment: Value previously returned in XML results from ESearch and EPost and used with EFetch in place of a primary UI result list.

    WebEnv=WgHmIcDG], etc.

Query_key:  The value used for a history search number or previously returned in XML results from ESearch or EPost.

    query_key=6 

Note: WebEnv is similar to the cookie that is set on a user's computers when accessing PubMed on the web.  If the parameter usehistory=y is included in an ESearch URL both a WebEnv (cookie string) and query_key (history number) values will be returned in the results. Rather then using the retrieved PMIDs in an ESummary URL you may simply use the WebEnv and query_key values to retrieve the records. WebEnv will change for each ESearch query.

Tool: A string with no internal spaces that identifies the resource which is using Entrez links (e.g., tool=flybase). This argument is used to help NCBI provide better service to third parties generating Entrez queries from programs. As with any query system, it is sometimes possible to ask the same question different ways, with different effects on performance. NCBI requests that developers sending batch requests include a constant 'tool' argument for all requests using the utilities.

    tool=

E-mail Address: If you choose to provide an email address, we will use it to contact you if there are problems with your queries or if we are changing software interfaces that might specifically affect your requests. If you choose not to include an email address we cannot provide specific help to you,  but you can still sign up for utilities-announce to receive general announcements.

    email=

PubMed

Record Identifier: UIs required if WebEnv is not used.

    id=11877539, 11822933,11871444

Current values:

    PubMed ID
    MEDLINE UI

Display Numbers: Used when the results from EPost or ESearch are maintained in the user's environment. The maximum number of retrieved records is 10,000.

    retstart=x (x= sequential number of the first id retrieved - default=0 which will retrieve the first record)
    retmax=y  (y= number of items retrieved - default=20)

Retrieval Mode:

    retmode=output format

Current values:

    xml
    html
    text
    asn.1

Use your web browser's View Page Source function to display results in xml retrieval mode.

Retrieval Type:

    rettype=output types based on database

Current values:

    uilist
    abstract
    citation
    medline
    full  (journals only) 

Not all Retrieval Modes are possible with all Retrieval Types.

PubMed Options:

	uilist 	abstract 	citation 	medline
xml 	x 	x* 	        x* 	         x*
text 	x 	x 	        x 	         x
html 	x 	x 	        x 	         x
asn.1 	n/a 	x* 	        x* 	         x

x = retrieval mode available
*returned retrieval type is the complete record in the retrieval mode
n/a - not available

Examples:
In PubMed display PMIDs 12091962 and 9997 in html retrieval mode and abstract retrieval type:
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=12345,9997&retmode=html&rettype=abstract

In PubMed display PMIDs from history statement in html retrieval mode and medline retrieval type (where x is replaced by WebEnv and query_key values):
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&WebEnv=xxxx&query_key=x&retmode=html&rettype=medline

In PubMed display PMIDs in xml retrieval mode:
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=11748933,11700088&retmode=xml

In Journals display records for journal IDs 22682,21698,1490:
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=journals&id=22682,21698,1490&rettype=full 

=cut



=over 4

=item EFetch()

B<Description:> Retrieves records in the requested format from a list of one or more UIs or from user\'s environment. 


B<Parameters:> ($UserAgent, $URLvals_href)

$UserAgent is a LWP::UserAgent object

$URLvals_href is a reference to a hash containing values for URL keys.  Available keys include:

db
WebEnv
query_key
tool
email
id
retstart
retmax
strand
seq_start
seq_stop
complexity
retmode
rettype





see details above for expected contents.


B<Returns:> $text

$text includes the unparsed data returned from GenBank based on the fetch request. 

=back

=cut


sub EFetch {
    my ($ua, $URLvals_href) = @_;
    
    my $base_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?";
    
    my @valid_keys = qw(db
			WebEnv
			query_key
			tool
			email
			id
			retstart
			retmax
			strand
			seq_start
			seq_stop
			complexity
			retmode
			rettype
			);
    
    return (&process_request($ua, $base_url, \@valid_keys, $URLvals_href));
}


##############################################################################################################

    ;

=head1 DESCRIPTION

ELink
Last updated: July 21, 2003

ELink:  Checks for the existence of an external or Related Articles link from a list of one or more primary IDs;  retrieves IDs and relevancy scores for links to Entrez databases or Related Articles; creates a hyperlink to the primary LinkOut provider for a specific ID and database, or lists LinkOut URLs and attributes for multiple IDs.
 

    * URL Parameters
          o Databases
                + PubMed, Protein, Nucleotide, Structure, Genome, PopSet, OMIM, Taxonomy, ProbeSet, 3D Domains, UniSTS, Domains, SNP, UniGene.
          o Record Identifier  Relative Dates  Date Ranges  Date Type Search Limits  Retrieval Mode  To Database  From Database   Web Environment Query_key  Link Command Values  Tool E-mail Address
            Examples
    * System Requirements for Utility Users
    * Entrez DTDs
    * Demonstration Program
    * Announcement Mailing List
    * Questions?

URL parameters:

Utility parameters may be case sensitive, therefore, use lower case characters in all parameters except for WebEnv.

Base URL: http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?

Current database values by category:

    Literature databases:

        omim - The OMIM database including the collective data is the property of the Johns Hopkins University, which holds the copyright.
        pubmed - Journal publishers hold the copyright on the abstracts in PubMed. NLM provides no legal advice concerning distribution of copyrighted materials.

    Sequence databases:

        nucleotide
        protein
        genome
        popset
        snp
        unists
        unigene 

    Structure databases:

        structure
        domains - equivalent to Entrez 3D Domains
        cdd - equivalent to Entrez Domains 

    Taxonomy database:

        taxonomy

    Gene expression database:

        geo

Record Identifier: Required if WebEnv is not used.

    id=12345,92932,828282

Current ID values:

    PMID (pubmed)
    MIM number (omim)
    GI number (nucleotide, protein)
    Genome ID (genome)
    Popset ID (popset)
    SNP cluster ID (snp)
    UniSTS ID (unists)
    UniGene cluster ID (unigene)
    MMDB-ID (structure)
    PSSM-ID (cdd)
    3D SDI (domains)
    TAXID (taxonomy)
    GEO ID (geo)

Date Parameters (only valid for dbfrom=pubmed&cmd=neighbor)

Relative Dates: Limit items to the number of days immediately preceding today's date.

    reldate=

For example:

      reldate=90
      reldate=365

Date Ranges:  Limit results bounded by two specific dates. Both mindate and maxdate are required if date range limits are applied using these variables.

    mindate=
    maxdate=

For example:

      mindate=2001
      maxdate=2002/01/01

Date Type:  Limit dates to a specific date field based on database.

    datetype=

For example:

      datetype=edat
      datetype=mdat

Search Limits: Limit the linked results with an Entrez query (only valid for cmd=neighbor).

    term=query

For example:

    term=medline[sb]
    term=human[organism]
    term=biomol genomic[properties]

Retrieval Mode:

    retmode=mode

For example:

    retmode=xml (default)
    retmode=ref (only used with cmd=prlinks)

Use your web browser's View Page Source function to display XML (default) results. Retmode=ref is only available for one ID.

To Database: ELink 'linking to' or destination database, pubmed is the default db.

    db=database

From Database: ELink 'linking from' or origination database, pubmed is the default db.

    dbfrom=database

Web Environment: Value previously returned in XML results from ESearch and EPost and used with ELink in place of primary ID result list. Required if id is not used.

    WebEnv=WgHmIcDG]B`&gt;&gt; etc.

Query_key: The value used for a history search number or previously returned in XML results from ESearch or EPost. Required if id is not used.

    query_key=6 

Link Command Values: Default is neighbor.

    cmd=command

For example:

    cmd=prlinks - List the hyperlink to the primary LinkOut provider for multiple IDs and database.
    cmd=prlinks&retmode=ref - Create a hyperlink to the primary LinkOut provider for a single ID and database.
    cmd=llinks - List LinkOut URLs and Attributes, except PubMed libraries, for multiple IDs and database.
    cmd=llinkslib - List LinkOut URLs and Attributes for multiple IDs and database.
    cmd=lcheck - Check for the existence (Y or N) of an external link in for multiple IDs and database.
    cmd=ncheck - Check for the existence of a neighbor link for each ID within a database, e.g., Related Articles in PubMed.
    cmd=neighbor - Display neighbors and their scores by database and ID.

Note:

    * Each ID set for lcheck and llinks is processed separately whether the IDs are listed or WebEnv is used.
    * Date limits are not valid for Link Commands.
    * The default Command Value is cmd=neighbor.

Tool:  A string with no internal spaces that identifies the resource which is using Entrez links (e.g. tool=flybase). This argument is used to help NCBI provide better service to third parties generating Entrez queries from programs. As with any query system, it is sometimes possible to ask the same question different ways, with different effects on performance. NCBI requests that developers sending batch requests include a constant 'tool' argument for all requests using the utilities.

    tool=

E-mail Address: If you choose to provide an email address, we will use it to contact you if there are problems with your queries or if we are changing software interfaces that might specifically affect your requests. If you choose not to include an email address we cannot provide specific help to you,  but you can still sign up for utilities-announce to receive general announcements.

    email=

Examples:

To retrieve IDs and relevancy scores from PubMed for PMID  9298984 to the PubMed database:
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&id=9298984&cmd=neighbor

To retrieve IDs from Nucleotide for GI 18250303, 18250301, 18250299 to Protein:
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nucleotide&db=protein&id=18250303,18250307

To retrieve PubMed related articles for PMIDs 11812492 11774222 with a publication date from 1995 to the present:
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&id=11812492,11774222&db=pubmed&mindate=1995&datetype=pdat

To retrieve MEDLINE indexed only related articles for PMID 12242737:
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&id=12242737&db=pubmed&term=medline[sb]

To create a hyperlink to the first link available for PMID 10611131 in PubMed:
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&id=10611131&retmode=ref&cmd=prlinks

To list all available links in PubMed, except for libraries, for PMIDs 12085856 and 12085853:
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&id=12085856,12085853&cmd=llinks

To check for the existence of a Related Articles link for PMIDs 0611131, 111645 and 12068369:
http://eutils.ncbi.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&id=10611131+111645&id=12068369&cmd=ncheck

 
=cut
 

    ;


=over 4

=item ELink()

B<Description:> Checks for the existence of an external or Related Articles link from a list of one or more primary IDs;  retrieves IDs and relevancy scores for links to Entrez databases or Related Articles; creates a hyperlink to the primary LinkOut provider for a specific ID and database, or lists LinkOut URLs and attributes for multiple IDs. 


B<Parameters:> ($UserAgent, $URLvals_href)

$UserAgent is a LWP::UserAgent object

$URLvals_href is a reference to a hash containing values for URL keys.  Available keys include:


id
reldate
mindate
maxdate
datetype
term
retmode
db
dbfrom
WebEnv
query_key
cmd
tool
email
    
see details above for expected contents.


B<Returns:> $text

$text includes the unparsed data returned from GenBank based on the link request. 

=back

=cut


sub ELink {
    my ($ua, $URLvals_href) = @_;
    
    my $base_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?";
    
    my @valid_keys = qw(id
			reldate
			mindate
			maxdate
			datetype
			term
			retmode
			db
			dbfrom
			WebEnv
			query_key
			cmd
			tool
			email
			);

    
    return (&process_request($ua, $base_url, \@valid_keys, $URLvals_href));
    
}


####
#private

sub process_request {
    my ($ua, $base_url, $valid_keys_aref, $URLvals_href) = @_;
    
    my $url = $base_url;
    
    my $append_flag = 0;
    foreach my $key (@$valid_keys_aref) {
	if (my $value = $URLvals_href->{$key}) {
	    if ($append_flag) {
		$url .= '&';
	    }
	    $url .= "$key=$value";
	    $append_flag = 1;
	}
    }

    ## Process URL
    
    print "URL: $url\n" if $SEE;
    
    my $req = new HTTP::Request (GET => $url);
    my $res = $ua->request($req);
    my $text = $res->content;
    
    return ($text);
    
}









1; #EOM



