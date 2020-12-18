#
# a simple routine to provide access checking
#

$OK_CHARS='\"\+\-a-zA-Z0-9_.@ \/%:';

$HOST_NAME="http://fasta.bioch.virginia.edu/";
$HOST_DIR="fasta/";
$CGI_DIR="fasta/cgi/";
$BIN_DIR="/seqprg/slib/bin/";
$FAST_LIBS="/slib0/lib/fastlibs";
$FAST_GNMS="/slib0/genomes/fastlibs";
$TMP_DIR="/www/doc/tmp/";
$GS_BIN="/usr/local/bin/gs";
$PPM_BIN="/usr/local/pbmplus/ppmtogif";

#@ALT_HOSTS = ("fasta.bioch.virginia.edu","wrpx1.bioch.virginia.edu");
@ALT_HOSTS = ("fasta.bioch.virginia.edu");

$msg1 = <<EOF1 ;
<html>
<head>
<title>FASTA Sequence Comparison</title></head>
<p>
&nbsp;
<h3> Sorry - some FASTA programs/databases are not available to: </h3>

EOF1
    ;
$msg2 = <<EOF2 ;
<p>
<font size=+1>
The FASTA package can be downloaded from <A HREF="ftp://ftp.virginia.edu/pub/fasta">
ftp://ftp.virginia.edu/pub/fasta</A></font>
<p>
<A HREF="mailto:wrp\@virginia.edu">wrp\@virginia.edu</A>
<p>

EOF2

sub CheckHost {
    my $host;

    $host = $ENV{'REMOTE_HOST'};
    if ($host eq "") {$host = $ENV{'REMOTE_ADDR'};}

    if ( $host =~ /virginia.edu$/ || $host =~/cstone.net$/ ) {return 1;} 
    else {return 0;}
}
#  
sub DenyHost {
    my $host;

    $host = $ENV{'REMOTE_HOST'};

    print $msg1;
    print "<h3> $host</h3>\n";
    print $msg2;
    print &HtmlBot; 
    exit();
}

sub Seq_Len 
{
    my $seq = @_;

    return length $seq;
    }

