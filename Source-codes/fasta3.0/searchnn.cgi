#! /usr/local/bin/perl

use LWP::Simple;

require "cgi-lib.pl";
require "my-cgi.pl";

# programs and the associate command lines:

@libp = ("libpa", "libpb", "libpn", "libpu", "libpp", "libpd", "libpg","libps",
	 "libpr", "libph", "libpo", "libpy", "libpe", "libpi", "libpj", "libpl",
	 "libpm", "libpc", "libpw", "libpk" );

$lib_s = "%";

$host_cgi = $HOST_NAME . $CGI_DIR;

%pgm_line = ("FASTA", "fasta3_t",
	     "TFASTX", "tfastx3_t",
	     "FASTX", "fastx3_t",
	     "TFASTY", "tfasty3_t",
	     "FASTY", "fasty3_t",
	     "FASTF", "fastf3_t",
	     "TFASTF", "tfastf3_t",
	     "FASTS", "fasts3_t",
	     "TFASTS", "tfasts3_t",
	     "SSEARCH", "ssearch3_t");

%mat_line = ("Blosum50", "-s BL50",
	     "Blosum62", "-s BL62",
	     "Blosum80", "-s BL80",
	     "Pam250", "-s P250",
	     "Pam120", "-s P120",
	     "MD20", "-s M20",
	     "MD10", "-s M10" );

#if ( &CheckHost == 0) { &DenyHost; }

&ReadParse;

$ENV{'FASTLIBS'} = $FAST_LIBS;

if ($in{"searchx"}) {
#microbial genomes

    $ENV{'FASTLIBS'} = $FAST_GNMS;
    if ($in{"searchx"} eq "genomes") {
	@libp = ("libpe", "libpc", "libpi", "libpj", "libpl", 
		 "libpm", "libph", "libpy", "libpw",
		 "libpa", "libpb", "libpd", "libpf", "libpr", "libpt"
		 );
	$search_url = "<A HREF=\"".$host_cgi . "searchfg.cgi?query=%s&db=%c&lib=%s&pgm=%s&start=%ld&stop=%ld&n1=%d\">Re-search database</A>&nbsp;&nbsp;";
	$ENV{'SRCH_URL'} = $search_url;
    }
    elsif ($in{"searchx"} eq "gen_subs") {
	foreach $l ( "libpa","libpp","libpd","libpn","libpk",
		     "libpq","libpr","libps","libpw" ) {
	    if ($in{$l} ne "") {
		$ENV{'FASTLIBS'} = $FAST_LIBS;
		last;
	    }
	}
    }
    @libn = ("libne", "libnc", "libni", "libnj", "libnl", 
	     "libnm", "libnh", "libny",
	     "libna", "libnb", "libnd", "libnf", "libnr", "libnt"
	     );

    %libnN =  ("libne"=>"-N 5000", "libnc"=>"-N 5000",
	       "libni"=>"-N 5000", "libnj"=>"-N 5000",
	       "libnl"=>"-N 5000", "libnm"=>"-N 5000",
	       "libnh"=>"-N 5000", "libny"=>"-N 5000",
	       "libna"=>"-N 5000", "libnb"=>"-N 5000",
	       "libnd"=>"-N 5000", "libnf"=>"-N 5000",
	       "libnt"=>"-N 5000", "libnr"=>"-N 5000"
	       );
}
else {
# genbank
    $ENV{'FASTLIBS'} = $FAST_LIBS;
    @libn = ("libnp", "libnr", "libnm", "libnb", "libne", "libnh","libni",
	     "libnl", "libnt", "libns", "libnv", "libng", "libnz",
	     "libnf", "libnj", "libnc", "libno", "libny");

    %libnN =  ("libnp"=>"", "libnr"=>"", "libnm"=>"", "libnb"=>"",
	       "libne"=>"", "libnh"=>"", "libni"=>"", "libnl"=>"",
	       "libnt"=>"","libns"=>"", "libnv"=>"", "libng"=>"",
	       "libnz"=>"", "libny"=>"-N 5000");
}

$p_name = $in{"program"};

print &PrintHeader;

#Check to see what variables are being entered
#print &HtmlTop;
#print &PrintVariables(%in);
#print &HtmlBot;

$ssr = $in{"ssr"};
if ($ssr) {$ssr = ":".$ssr;}

$m_name = $in{"pmatrix"};
if ($m_name eq "Default") {$m_name="";}

$seqtype = $in{"seqtype"};

$r_name = $in{"dmatrix"};
if ($r_name eq "Default") {$r_name="";}
else {
    if ($m_name eq "" && $seqtype == 1) {$seqtype = 2;}
    if ($r_name eq "blastn2") {$r_name="+4/-12";}
}

if ($in{"ktup"}) {
    $ktup = "$in{\"ktup\"}";
    $ktup =~ s/\D//g;
}
else {$ktup = "";}

if ($in{"eval"}) {
    $eval = "$in{\"eval\"}";
    $eval =~ s/[^\d\.eE\-]//g;
    $eval = "-E ".$eval;
}
else {$eval = "";}

if ($in{"etop"}) {
    $etop = "$in{\"etop\"}";
    $etop =~ s/[^\d\.eE\-]//;
    $etop = "-F ".$etop;
}
else {$etop = "";}

if ($in{"best"}) {
    $best = "$in{\"best\"}";
    $best =~ s/[^\d]//g;
    $best = "-b ".$best;
}
else {$best = "";}

if ($in{"align"}) {
    $align = "$in{\"align\"}";
    $align =~ s/[^\d]//g;
    $align = "-d ".$align;
}
else {$align = "";}

if ($in{"out_opt"}) {
    $out_opt = $in{out_opt};
    $out_opt =~ s/;\|()//g;
}
else {$out_opt = ""};

if ($in{"gap"}) {
    $gap = "$in{\"gap\"}";
    $gap =~ s/[^\d\-]//g;
    $gap = "-f ".$gap;
}
else {$gap = "";}

if ($in{"ext"}) {
    $ext = "$in{\"ext\"}";
    $ext =~ s/[^\d\-]//g;
    $ext = "-g ".$ext;
}
else {$ext = "";}

$mat_type = "";
if ($seqtype == 1) {
    if ($m_name ne "") {$mat_type=$mat_line{$m_name};}
}
else {
    if ($r_name ne "") {$mat_type = "-r \"".$r_name."\"";}
}

$rv_type = "";
$db="p";

if ($p_name eq "FASTA" || $p_name eq "SSEARCH") {
#Check to see whether the sequence is DNA or Protein
    if ($seqtype == 2 || $seqtype == 3 || $seqtype == 4) {
	if ($seqtype == 3) {$rv_type .= " -3";}
	if ($seqtype == 4) {$rv_type .= " -i";}
	$db="n";
	$seqtype="-n ";
	foreach $l ( @libn ) {
	    if ("$in{$l}" ne "") { $lib_s = $lib_s."$in{\"$l\"}"; }
	}
    }
    else {
	$seqtype = "";
	foreach $l ( @libp ) {
	    if ("$in{$l}" ne "") { 
		$lib_s = $lib_s."$in{\"$l\"}";
	    }
	}
    }
}
elsif ($p_name eq "FASTF" || $p_name eq "FASTX" || $p_name eq "FASTY") {
    if ($p_name eq "FASTX" || $p_name eq "FASTY") {
	if ($seqtype == 3) {$rv_type .=" -3";}
	if ($seqtype == 4) {$rv_type .=" -i";}
	$db="n";
    }
    $seqtype = "";
    foreach $l ( @libp ) {
	if ("$in{$l}" ne "") { $lib_s = $lib_s."$in{\"$l\"}"; }
    }
}
else {   # TFASTF/TFASTX/TFASTY
    $seqtype = "";
    foreach $l ( @libn ) {
	if ("$in{$l}" ne "") { $lib_s = $lib_s."$in{\"$l\"}";
			       $Nopt = $libnN{$l};
 }
    }
}

if ($lib_s eq "%") {
    if ($seqtype eq "-n") {$lib_s .= "p"}
    else {$lib_s .= "q";}
}

select(STDOUT), $| = 1;

$seq_form = $in{"in_seq"};

$seq_in = $in{"sequence"};

# check if seq_in is really an accession number
if ($seq_form !~ /Accession/ && $seq_in !~ /^>/ && ($seq_in =~ /_/ || $seq_in =~/\|/)) {
    $seq_form = "Accession";
#    print "<p>Changing to Accession\n<p>\n";
}

#print "seq_form: $seq_form<br>\n";
#print "seq_in: $seq_in<br>\n";

if (($p_name =~ /FASTF/) && ($seq_form =~ /Accession/)) {
    print &HtmlTop ("$seq_in/FASTF ERROR");
    print "<h2>FASTF requires partial sequence data</h2><br>\n";
    print &HtmlBot;
    exit;
}

#print "p_name: $p_name\tseqtype: $seqtype\n";
if (($in{"searchx"} eq "") && (! &CheckHost)
    && ($p_name eq "SSEARCH" || $seqtype eq "-n")) {
    &DenyHost;
}

$pgm_cmd = $BIN_DIR . $pgm_line{$p_name};
$pgm_cmd = $pgm_cmd . " $Nopt $seqtype $rv_type $mat_type $eval $etop $best $align -w 80 -m 6 $gap $ext -q $out_opt \@$ssr $lib_s $ktup";

# old way
# $pgm_cmd =~ s/[;><&\*`\|]//g;
# better way
$pgm_cmd =~ s/[^$OK_CHARS]/_/go;

if ($seq_form =~ /^FASTA/) {

    print "<pre>$pgm_cmd<pre>\n";

    open(TEMP, "| $pgm_cmd");

    if ($seq_in !~ /^>/) {
	print TEMP ">QUERY sequence\n";
    }
    print TEMP $seq_in;
    print TEMP "\n";
    close(TEMP);
}
elsif ($seq_form =~ /Accession/) {
    $url="http://www.ncbi.nlm.nih.gov/htbin-post/Entrez/query?db=$db&form=6&Dopt=f&html=no&uid=$seq_in";

    $entry_line = get $url;
    if ($entry_line eq "") {goto Error1;}

    @entry = split(/\n/,$entry_line);

#    print &HtmlTop ("$seq_in ERROR");
#    print "<pre>\n";
#    print "<pre>%%$entry_line%%</pre>\n";
#    foreach $l ( @entry ) { print "***$l\n";}
#    print "</pre>\n";
#    print &HtmlBot;
#
    if ($entry_line =~/^\*\*\* No Documents Found \*\*\*/ ||
	$entry_line =~/ERROR/) {
Error1:
	print &HtmlTop ("$seq_in ERROR");

	print "\n<tt><pre>$url</tt></pre>\n";
	print "\n<h2>ERROR - $seq_in not found</h2>\n";
	print &HtmlBot;
    }
    else {
	print "<pre>$pgm_cmd</pre>\n";
#	print "<pre>FASTLIBS: $FAST_LIBS</pre><br>\n";

	open(TEMP, "| $pgm_cmd");

	$have_gt = 0;
	for $l (@entry) {
	    if ($l=~/^>/) {$have_gt = 1;}
	    if ($have_gt == 1) { print TEMP "$l\n";}
	}
	close (TEMP);
    }
}
else {
    print &HtmlTop;
    print "<h2>ERROR</h2>\n";
    print "invalid option: $seq_form\n";
    print &HtmlBot;
}

