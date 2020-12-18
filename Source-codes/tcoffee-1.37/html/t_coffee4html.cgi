#!/igs/local/bin/perl
#Perl Wrapper for the html version of T-Coffee
#Environement variables
 


$USER_BIN="$HOME/bin/";


#MAIN VARIABLES:
      use Env qw(HOME);
      if ( $HOME eq ""){$HOME=`pwd`;chop ($HOME);$HOME=~s/public_html.*//;}

      $USER_BIN="$HOME/bin/";
      $WWW_HOME="http://igs-server.cnrs-mrs.fr/~cnotred/";
      $WWW_FILE_HOME="$HOME/public_html/";
      $USER_BIN="$HOME/bin/";
      $html_page="/t_coffee.html";
      $programme="t_coffee";
      $SCRATCH_AREA="$WWW_FILE_HOME/tmp/";
      $CALLING_PAGE    ="$WWW_HOME/$html_page";



#Set the programme so that it uses the right bin
#[CUSTOMIZE HERE FOR YOUR ARCHITECTURE]
     $machine=`/bin/uname -s`;
   
     chop ($machine);
	$ENV{MAC}=$machine;
	use Env qw(MAC);
	if ($MAC eq "IRIX64" || $MAC eq "SGI64" || $MAC eq "SGI" || $MAC eq "IRIX") 
            {$ENV{OS}="Irix";}
	elsif ($MAC eq "SunOS" || $MAC eq "SUN")
            {$ENV{OS}="Sun";}
	elsif ($MAC eq "Linux")
            {$ENV{OS}="Linux";}    
        else
            {
	    &output_error_page ( "The machine $MAC on which this program should run is unknown[FATAL]\n");
	    }
	use Env qw(OS);

#Define the Programs and various addresses [CUSTOMIZE HERE]
    
#Your binaries

    $ENV{USER_BIN}="$HOME/bin/$OS/";#Your binaries
    $ENV{ALN_TAB}="$USER_BIN";
    use Env qw(USER_BIN);
    use Env qw(PATH);
    $ENV{PATH}=".:/usr/sbin:/usr/bsd:/sbin:/usr/bin:/bin:/usr/bin/X11:$USER_BIN"; #The Binary Paths
    $programme="$programme";


#check the required programs are installed
$t_coffee_exec=`which t_coffee`;
chop($t_coffee_exec);

if ( !(-e $t_coffee_exec))
    {
  
    &output_error_page ( "$programme  is not Installed home=[$HOME] [$t_coffee_exec] Path=[$PATH][FATAL]\n", $CALLING_PAGE);
    die;
    }


if ($ENV{'REQUEST_METHOD'} eq "GET")
    {
    $request=$ENV {'QUERRY_STRING'};
    @parameter_list=split(/&/,$request);
    foreach (@parameter_list)
        {
	($name, $value)=split (/=/);
	if (!$arg{$name}{number}) {$arg{$name}{number}=0;}

	$arg{$name}{$arg{$name}{number}}=$value;
	$arg{$name}{number}++;
        }
    $request=$ENV{'QUERRY_STRING'};
    }
elsif ($ENV{'REQUEST_METHOD'} eq "POST")
    {
     read(STDIN, $request, $ENV{'CONTENT_LENGTH'}) || output_error_page( "Could not read querry\n", $CALLING_PAGE);
     @parameter_list=split (/-{2,}\d+.*\n/, $request);
     foreach $val (@parameter_list) 
	 {	
	 $val=~s/(Content.*)\n.*\n//g;	 
	 $param_line=$1;
	 
	 if ($param_line=~/filename/) 
	     {	 
	     
	     $param_line=~/.*name="(.*)"; filename="(.*)".*/;
	     $param_name="$1";
	     if ($2 ne "" && $2 ne "no_value")
		 {
	         if (($file=&make_random_file($val)) ne "no_file")
	            {
		    $arg{$param_name}{$arg{$param_name}{number}++}=$file;
		    if ( $param_name eq "-in"){$sequence_in=1;}
		    }
		 }
	     }
	 else 
	     {
	     
	     $param_line=~/.*name="(.*)".*/;	
	     $param_name="$1";

	     if ($param_name eq "data") 
		 {
		 if (($file=&make_random_file($val)) ne "no_file") 
		     {
		     $arg{"-in"}{$arg{"-in"}{number}++}=$file;
		     $sequence_in=1;
		     }
	         }
	     elsif ($param_name=~/\B-[^-]*/)
		 {	         
	         $val=~s/\s//g;
		 if ( !$arg{flag}{$param_name}){$arg{flag}{$param_name}=1;$flag_list[$nflags++]=$param_name;}
		 if ( $val ne "" && $val ne "no_value") 
		     {
		     $arg{$param_name}{$arg{$param_name}{number}++}=$val;
		     }
		 }
	     else  
	         {
	         $val=~s/\s//g;
		 $arg{$param_name}{0}=$val;
		 }
	     
	     }
	 }
 }
 

#1 CHECK FOR COMPULSORY PARAMETERS
if ( $sequence_in!=1)
   {
   &output_error_page ( "You did not provide any Data", $CALLING_PAGE);
   die;
   }

if ( !$arg{'e-mail'}{0})
   {
   &output_error_page ( "You did not provide your Email", $CALLING_PAGE);
   die;
   }

#2 Set The Run Name
$flag_list[$nflags++]="-run_name";
$arg{-run_name}{number}=1;
$arg{-run_name}{0}=&make_random_file_name;
$run_name=$arg{-run_name}{0};

#2 PROCESS THE PARAMETER LIST

$command="#!\/bin\/sh\/n";
$command=$command."source $HOME/.cshrc\n";
$command=$command.$programme;
foreach $flag (@flag_list)
	{
        $command="$command $flag=";
        for ( $a=0; $a<$arg{$flag}{number};$a++)
            {$command="$command $arg{$flag}{$a}";}
	}


#3 list of the output files

$output_file{++$n_output_file}{name}="$run_name"."_result.html";
$output_file{  $n_output_file}{description}="Results";
$output_file{  $n_output_file}{format}="html";


$output_file{++$n_output_file}{name}="$run_name.tc_log";
$output_file{  $n_output_file}{description}="log of t_coffee";
$output_file{  $n_output_file}{format}="ascii";

for ( $a=0; $a< $arg{'-output'}{number}; $a++)
    {
    if ($arg{'-output'}{$a} eq "clustalw"){$extention="aln";}
    else {$extention=$arg{'-output'}{$a}};
	
    $output_file{++$n_output_file}{name}="$run_name.$extention";
    $output_file{  $n_output_file}{description}="Multiple Alignment";
    $output_file{  $n_output_file}{format}="$extention";
    }

$output_file{++$n_output_file}{name}="$run_name.dnd";
$output_file{  $n_output_file}{description}="Dendrogram used for Progressive Alignment";
$output_file{  $n_output_file}{format}="Phyllip";
$output_file{n_output_file}=$n_output_file;
$output_file{command}=$command;

&submit ("$command -out_lib=no -quiet stderr >$run_name.tc_log 2>$run_name.tc_log" );
%output_file=&make_output (%output_file);    


#Make the Output

print "Content-type: text/html\n\n";	
print "<html><body  bgcolor=#FFFFFF><table bgcolor=#FFFFFF  border=8 cellpadding=8 align=center width=60%>";
print "<tr><td><h4>The computation is now finished. Click <a href=$output_file{1}{url}>[HERE]</a>to see the results<\/h4><\/td><\/tr><\/table><\/body><\/html>";


#Clean
for  ( $a=0; $a<$n_tmp_files; $a++ ) 
     {
     unlink ($tmp_files[$a]);
     }


sub submit
    {
    my $command=@_[0];
    
    system $command;
    }

sub make_random_file 
    {
    my $content=@_[0];
    my $name;
    my $cnt;

    if (($content=~/\S/g)==0) 
	{	
	return "no_file";
        }
    $rand_number=int(rand 100000)+1 ;
    $name= &make_random_file_name();
    $name="$name.tcoffee_file";

    $tmp_files[$n_tmp_files++]=$name;    
    open NEW_FILE,">$name";
    
    print NEW_FILE "$content";
    close (NEW_FILE);
    $cnt=chmod 0777, $name;
    return $name;
   }

sub make_random_file_name
    {
    my $name;
    my $rand_number;

    $rand_number=int(rand 100000)+1 ;
    $name="$SCRATCH_AREA/$rand_number"."$$";
    return $name;
    }

sub output_error_page
    {
my $myerror=@_[0];
my $back=@_[1];

print "Content-type: text/html\n\n";
print "<html>\n";
print "<head>\n";
print "<H1><pre>    ERROR   </pre></H1>\n";
print "<H3><pre>$myerror\n\n</pre></H3>\n";
print "<a href=$back>BACK</a></pre>\n";
print "</head>\n";
print "</html>\n";

exit;
}
    
sub make_output
	{
	my %loc_output=@_;
	my $n;

	$n=$loc_output{'n_output_file'};
	for ($a=0; $a<=$n; $a++)
            {
	    if ($loc_output{$a}{name}=~/score_html/)
	       {
	       $new_name=$loc_output{$a}{name};
	       $new_name=~s/score_html/html/;
	       open (IN,  $loc_output{$a}{name});
	       open (OUT,  ">$new_name");
	       while(<IN>){print OUT "$_";}
	       close (IN);
	       close (OUT);
	   
	       unlink ( $loc_output{$a}{name});
	       $loc_output{$a}{name}=$new_name;
	       }
	    chmod 0777, $loc_output{$a}{name};
	    $loc_output{$a}{url}=$loc_output{$a}{name};
	    $loc_output{$a}{url}=~s/$WWW_FILE_HOME/$WWW_HOME/;
	    }
	

	open (RESULT, ">$loc_output{1}{name}");
       
	
	
	print RESULT  "<html><body  bgcolor=#FFFFFF><table bgcolor=#FFFFFF  border=8 cellpadding=8 align=center width=60%>";
        print RESULT  "<tr><td align=center><h1>Ouput of T-Coffee<\/td><\/tr><\/table>";
	
	print RESULT  "<H4> Your data will remain available on this server over the next 3 days. It will then be deleted. Do not forget to bookmark this <a href=$loc_output{1}{url}>URL</a>($loc_output{1}{url}) or save it for further reference </h4>";

        print RESULT  "<table bgcolor=#FFFFFF border=8 cellpadding=8 align=center width=80%>";
        print RESULT  "<tr align=center><td><h1>Type<\/h1><\/td><td><h1>Format<\/h1><\/td><td><h1>Status<\/h1><\/td><\/tr>";


	for ($a=0; $a<=$n; $a++)
            {
	    if  ( -e $loc_output{$a}{name})
	         {print RESULT  "<tr><td><h4><A href=$loc_output{$a}{url}>$loc_output{$a}{description}<\/a><\/td><td>$loc_output{$a}{format}<\/td><td><font color=green> PRODUCED</font><\/tr>";}
	    elsif ($loc_output{$a}{name})
	         {print RESULT  "<tr><td><h4>$loc_output{$a}{description}<\/td><td>$loc_output{$a}{format}<\/td><td><font color=red>NOT PRODUCED</font><\/tr>";}

	    }
 
	print RESULT  "</table>";
	
	print RESULT  "<H1><a href=$CALLING_PAGE>BACK<\/a><\/H1>"; 
	print RESULT  "</body>";
	print  RESULT "</html>";

	close RESULT;
	chmod 0777, $loc_output{1}{name};
	return %loc_output;
	}

sub output_file
	{
	my $file_type    =@_[0];
	my $file_www_ref =@_[1];
	my $file_format  =@_[2];
	
	
	print "<pre>\t<a href=$file_www_ref>$file_type<\/a><\/pre>\n";
	return;
        }

sub html_output 
	{
	if (!$header) 
	    {
	    print "Content-type: text/html\n\n<html>\n<head><H1>T-COFFEE RESULTS</H1></head>"; 
	    $header=1;
	    }
       print "@_[0]\n";
       }

