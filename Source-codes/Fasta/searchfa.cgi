#! /usr/local/bin/perl
#
# copyright (c) 1997 William R. Pearson and the U. of Virginia
#
use LWP::Simple;

require "cgi-lib.pl";
require "my-cgi.pl";

%pgm_name = ("fx" => "fastx3_t","tx"=>"tfastx3_t",
	     "fa"=>"fasta3_t","tf"=>"tfasta3_t","sw"=>"ssearch3_t");

#if (! &CheckHost) {&DenyHost;}

# programs and the associate command lines:

&ReadParse;

# test conditions
#$in{"query"} = "GTM1_HUMAN";
#$in{"lib"} = "a";
#$in{"pgm"} = "fa";
#$in{"db"} = "p";

$ENV{'FASTLIBS'} = $FAST_LIBS;

#print &PrintHeader;
#print "<html>\n";
#print "<head>\n";
#print "<center>\n";
#print "<title> FASTX Results </title>\n";
#print "</title>\n";
#print "</center>\n";
#print "</head>\n";

#Check to see what variables are being entered

#print &PrintHeader;
#&HtmlTop ("Entrez: $query");
#print &PrintVariables(%in);

select(STDOUT), $| = 1;

#print "<tt><pre>\n";

$query = $in{"query"};
$db = $in{"db"};
$lib = $in{"lib"};
if ($opt_str = $in{"opt_str"}) {
    $opt_str =~ s/;\|()//g;
}
else {$opt_str = ""};

#$lib =~ s/%//;
$func=$in{"pgm"};
if (!func || $func eq "" ) {$pgm= $BIN_DIR . "fasta3_t";}
else {$pgm= $BIN_DIR . $pgm_name{$func};}

$sq_type = "";
if ($func eq "fa" && $db eq "n") {
    $sq_type="-n";
}

$rev_comp = "";
if ($query =~ /_r$/) {
    $rev_comp = "-i";
    $query =~ s/_r$//;
}

$query =~ s/_\d$//;

#
# check on the size, alignment start/stop
# older html files may not have this info
#

$q_str="@";
$n1 = $in{"n1"};
if ($n1 eq "" || $n1 >= 2000) {
    $start = $in{"start"};
    $stop = $in{"stop"};
    if ($start > $stop) {$tmp = $start; $start = $stop; $stop = $tmp;}
    if ($start < 200) {$start = 1;}
    else {$start -= 200};
    if ($n1 ne "" && $n1 - $stop < 200) {$stop = $n1;}
    else {$stop += 200;}
    $q_str .=":$start-$stop";
}

$acc = $in{"acc"};

#
# check that it is likely to be an ACC number
#
if ($acc && $acc =~ /^[A-Z][A-Z0-9]+/) {
    $url="http://www.ncbi.nlm.nih.gov/htbin-post/Entrez/query?db=$db&form=6&Dopt=f&html=no&uid=$acc&uid=$query";
}
else {
    $url="http://www.ncbi.nlm.nih.gov/htbin-post/Entrez/query?db=$db&form=6&Dopt=f&html=no&uid=$query";
}

#print "<pre>\n";
#print "$url\n";
#print "</pre><br>\n";


@entry = split(/\n/,(get $url));

#print "<pre>\n";
#foreach $n ( @entry ) {print "$n\n";}
#print "</pre>\n";

$pgm_cmd = "$pgm $sq_type -b 100 -d 50 -m 6 -w 80 $rev_comp $opt_str -q -w 80 $q_str $lib";
$pgm_cmd =~ s/[^$OK_CHARS]/_/go;

if ($entry[2] =~/^\*\*\* No Documents Found \*\*\*/ || $entry[0] =~/^ERROR/) {

    print &PrintHeader;

    &HtmlTop ("$query ERROR");

    print "\n<tt><pre>$url</tt></pre>\n";
    print "\n<h2>ERROR - $query not found</h2>\n";

    &HtmlBot;
}
else {
    print &PrintHeader;
    print "| $pgm_cmd\n";
    open(TEMP,"| $pgm_cmd");

    $have_gt = 0;
    for $l (@entry) {
	if ($l=~/^>/) {$have_gt = 1;}
	if ($have_gt == 1) { print TEMP "$l\n";}
    }

    close (TEMP);
}

#print "</pre></tt>\n";
#print &HtmlBot;

__END__
