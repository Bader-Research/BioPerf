#!/usr/bin/env perl
chdir "./example/test";
system "../../bin/t_coffee test.pep -in fast_pair -outfile=new_reference_test.aln >/dev/null";
$string=`diff reference_test.aln new_reference_test.aln|wc`;
@val=($string=~/(\w+)\s*/g);
if ($val[0]<=4 && -e "new_reference_test.aln"){print "\nInstallation of t_coffee Successful\n";}
else {print "\nInstallation of t_coffee Not Successful\n";}
