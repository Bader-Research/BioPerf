#!/bin/bash
source $HOME/.bioperf
echo Please specify your architecture 
echo [A] Alpha
echo [P] PowerPC
echo [X] x86
echo '[H] Custom. (If you have already ran the install-codes.sh successfully' 
echo 'to completion, you can choose this option to run the binaries installed' 
echo 'on your system)'
read arch


case $arch in 
"P" ) 
	arch2=PowerPc;;
"X" ) 
	arch2=x86;;
"H" ) 
	if [ ! -d $BIOPERF/Binaries/$HOSTNAME-Binaries ]; then
	    echo Please make sure install-codes.sh ran successfully 
	    echo There is no directory named $HOSTNAME-Binaries (created by install-codes.sh) in $BIOPERF/Binaries   
	    exit 
	fi
	if [ ! -d $BIOPERF/Scripts/$HOSTNAME-scripts ]; then
	    echo Please make sure install-codes.sh ran successfully 
	    echo There is no directory named $HOSTNAME-Scripts (created by install-codes.sh ) in $BIOPERF/Scripts       
	    exit 
	fi 
	arch2=$HOSTNAME;;
"A" ) 
	if [ -z $SIMDIR ]; then       
	    echo $'\n'
	    echo It appears that you have not installed sim-outorder executable. 
	    echo Download and install the simplesim package "http://www.simplescalar.com/agreement.php3?simplesim-3v0d.tgz" 
	    echo sim-outorder is part of this package
	    echo Set the environment variable SIMDIR to the directory containing sim-outorder
	    echo $'\n'
	    exit
	fi
	if [ -e $SIMDIR/sim-outorder ] 
	    then  echo ""
	else
	    echo $'\n'  
	    echo There is no file named sim-outorder in $SIMDIR
	    echo Please check for sim-outorder again in $SIMDIR
	    echo $'\n'
	    exit
	fi
	cd $BIOPERF/Scripts/Alpha-scripts
	for i in clustalw dnapenny fasta glimmer hmmer predator promlk 
	  do
	  ./alpha-$i.script
	done 
	;;
esac

BLAST=1
CE=1
CLUSTALW=1
FASTA=1
GLIMMER=1
GRAPPA=1
HMMER=1
PHYLIP=1
PREDATOR=1
TCOFFEE=1

echo How would you like to run BioPerf 
echo [C] You want to choose codes you would like to run 
echo [A] Running all codes.
read input3

if [ $input3 = "C" ] ; then  
    echo "Would you like to run Blast (blastn/blastp) [y/n]" 
    read input
    if [ $input = "n" ] ; then 
	BLAST=0 
    fi
    echo 'Would you like to run Ce (ce) [y/n]'
    read input
    if [ $input = "n" ] ; then 
	CE=0 
    fi 
    echo 'Would you like to run Clustalw (clustalw/clustalw-smp) [y/n]' 
    read input
    if [ $input = "n" ] ; then 
	CLUSTALW=0
    fi 
    
    echo 'Would you like to run Fasta (fasta/ssearch: fasta has only one input which takes a long time) [y/n]' 
    read input
    if [ $input = "n" ] ; then 
	FASTA=0  
    fi 
    echo 'Would you like to run Glimmer (glimmer/glimmer-package) [y/n]' 
    read input
    if [ $input = "n" ] ; then 
       GLIMMER=0 
    fi
    echo 'Would you like to run Grappa (grappa) [y/n] '
    read input
    if [ $input = "n" ] ; then 
	GRAPPA=0 
    fi
echo 'Would you like to run Hmmer (hmmsearch/hmmpfam ) [y/n] '
read input
if [ $input = "n" ] ; then 
    HMMER=0 
fi
echo 'Would you like to run Phylip (dnapenny/promlk: dnapenny has only one input which takes a long time) [y/n]' 
read input
if [ $input = "n" ] ; then 
    PHYLIP=0 
fi
echo 'Would you like to run Predator (predator) [y/n]' 
read input
if [ $input = "n" ] ; then 
    PREDATOR=0 
fi
echo 'Would you like to run Tcoffee (tcoffee) [y/n]' 
read input
if [ $input = "n" ] ; then 
    TCOFFEE=0 
fi
fi



echo Please give your input size 
echo [A] class-A
echo [B] class-B
echo [C] class-C
echo [X] Run all
read input 
case $input in 
"A" ) input2="-small"
	if [ $FASTA -eq "1" ] ; then
	    if [ -z $DATABASES ]; then
		echo $'\n'
		echo Fasta requires nr  database for running
		echo Download and install nr from  "ftp://ftp.ncbi.nih.gov/blast/db/nr.tar.gz"
		echo Then set the DATABASES environment variable to the directory containing nr
		echo $'\n'
		echo Would you like to run BioPerf without running Fasta [y/n]
		read user1
		if [ $user1 =  "y" ] ; then
		    FASTA=0
		else 
		    exit
		fi
	    else
		if [ -e $DATABASES/nr ]
		    then  echo ""
		else
		    echo $'\n'  
		    echo There is no file named nr in $DATABASES
		    echo Please check for nr again in $DATABASES
		    echo $'\n'
		    exit
		fi
	    fi
	fi
	;;

"B" ) input2="-medium"
        if [ $BLAST -eq "1" ] ; then 
	    if [ -z $DATABASES ]; then
		echo $'\n'
		echo Blast for medium input requires swissprot.aa database 
		echo Download and install swissprot.aa from  "www.bioperf.org/databases"
		echo Then set the DATABASES environment variable to the directory containing swissprot.aa
		echo $'\n'
		echo Would you like to run BioPerf without running Blast and Fasta[y/n]
		read user1
		if [ "$user1" = "y" ] ; then
		    BLAST=0
		else
		    exit
		fi
	    else 
		if [ -e $DATABASES/swissprot.aa ]
		    then  echo ""
		else
		    echo $'\n'  
		    echo There is no file named swissprot.aa in $DATABASES
		    echo Please check for swissprot.aa  again in $DATABASES
		    echo $'\n'
		    exit
		fi
	    fi
        fi
        if [ $FASTA -eq "1" ] ; then 
	    if [ -z $DATABASES ]; then
		echo $'\n'
		echo Fasta has only one class of inputs which  requires nr  database 
		echo Download and install nr from  "ftp://ftp.ncbi.nih.gov/blast/db/nr.tar.gz"
		echo Then set the DATABASES environment variable to the directory containing nr
		echo $'\n'
		echo Would you like to run BioPerf without running Fasta [y/n]
		read user2
		if [ "$user2" =  "y" ] ; then
		    FASTA=0
		else
		    exit
		fi
	    else
		if [ -e $DATABASES/nr ]
		then  echo ""
		else
		    echo $'\n'  
		    echo There is no file named nr in $DATABASES
		    echo Please check for nr again in $DATABASES
		    echo $'\n'
		    exit
		fi
	    fi	    
	fi
	





 

;;



"C" ) input2="-large"
	if [ $BLAST -eq "1" ] ; then  
	    if [ -z $DATABASES ]; then
		echo $'\n'
		echo Blast for medium input requires swissprot.aa database 
		echo Download and install swissprot.aa from  "www.bioperf.org/databases"
		echo Then set the DATABASES environment variable to the directory containing swissprot.aa
		echo $'\n'
		echo Would you like to run BioPerf without running Blast [y/n]
		read user1
		if [ "$user1" = "y" ] ; then
		    BLAST=0
		else 
		    exit
		fi
	    else 
		if [ -e $DATABASES/swissprot.aa ]
		    then  echo ""
		else
		    echo $'\n'  
		    echo There is no file named swissprot.aa in $DATABASES
		echo Please check for swissprot.aa  again in $DATABASES
		echo $'\n'
		exit
		fi
	    fi
	fi
	if [ $FASTA -eq "1" ] ; then
	    if [ -z $DATABASES ]; then
		echo $'\n'
		echo Fasta has only one class of inputs which  requires nr  database 
		echo Download and install nr from  "ftp://ftp.ncbi.nih.gov/blast/db/nr.tar.gz"
		echo Then set the DATABASES environment variable to the directory containing nr
		echo $'\n'
		echo Would you like to run BioPerf without running Fasta [y/n]
		read user1
		if [ "$user1" =  "y" ] ; then
		    FASTA=0
		else
		    exit
		fi
	    else
		if [ -e $DATABASES/nr ]
		    then  echo ""
		else
		    echo $'\n'  
		    echo There is no file named nr in $DATABASES
		    echo Please check for nr again in $DATABASES
		    echo $'\n'
		    exit
		fi
	    fi	
	fi

;;
"X" ) 		if [ $BLAST -eq "1" ] ; then
	if [ -z $DATABASES ]; then
	    echo $'\n'
	    echo Blast for medium input requires swissprot.aa database 
	    echo Download and install swissprot.aa from  "www.bioperf.org/databases"
	    echo Then set the DATABASES environment variable to the directory containing swissprot.aa
	    echo $'\n'
	    echo Would you like to run BioPerf without running Blast [y/n]
	    read user1
	    if [ "$user1" = "y" ] ; then
		BLAST=0
	    else
		exit
	    fi
	else 
	    if [ -e $DATABASES/swissprot.aa ]
		then  echo ""
	    else
		echo $'\n'  
		echo There is no file named swissprot.aa in $DATABASES
		echo Please check for swissprot.aa  again in $DATABASES
		echo $'\n'
		exit
	    fi
	fi
        fi
			if [ $FASTA -eq "1" ] ; then
	if [ -z $DATABASES ]; then
	    echo $'\n'
	    echo Fasta has only one class of inputs which  requires nr  database 
	    echo Download and install nr from  "ftp://ftp.ncbi.nih.gov/blast/db/nr.tar.gz"
	    echo Then set the DATABASES environment variable to the directory containing nr
	    echo $'\n'
	    echo Would you like to run BioPerf without running Fasta [y/n]
	    read user1
	    if [ "$user1" =  "y" ] ; then
	    FASTA=0
	    else
		exit
	    fi
	else
	    if [ -e $DATABASES/nr ]
		then  echo ""
	    else
		echo $'\n'  
		echo There is no file named nr in $DATABASES
		echo Please check for nr again in $DATABASES
		echo $'\n'
		exit
	    fi
	fi	    
        fi
	;;
esac

for i in Blast Ce Clustalw Fasta Glimmer Grappa Hmmer Phylip Predator Tcoffee
  do 
  case $i in
"Blast") 
	  if [ $BLAST -eq "1" ] ; then
	      echo "*********************************"
	      echo Running $i
	      echo "*********************************"
	      cd $BIOPERF/Scripts/$arch2-scripts/$i/blastn
	      starttime=`date +%s`
	      ./blastn-script  
	      endtime=`date +%s`
	      if [ $? -ne 0 ] ; then
		  echo $i did not run successfully 
	      exit
	      fi
	      blastntotal=`expr $endtime - $starttime`
	      cd $BIOPERF/Scripts/$arch2-scripts/$i/blastp
	    starttime=`date +%s`
	    ./blastp$input2.script
 endtime=`date +%s`
	      if [ $? -ne 0 ] ; then
		  echo $i did not run successfully 
	      exit 
	      fi 
 blastptotal=`expr $endtime - $starttime`
	  fi
	  ;;
"Ce")  if [ $CE -eq "1" ] ; then 
 echo "*********************************"
	      echo Running $i
	      echo "*********************************"
 cd $BIOPERF/Scripts/$arch2-scripts/$i
starttime=`date +%s`
endtime=`date +%s`
	  ./ce-script
	  if [ $? -ne 0 ] ; then                                                                                                                      
	      echo $i did not run successfully                                                                                                          
	      exit                                                                                                                                         
	  fi  
	   cetotal=`expr $endtime - $starttime`
	  fi
	  ;;
"Clustalw") 
 if [ $CLUSTALW -eq "1" ] ; then 
 echo "*********************************"
	      echo Running $i
	      echo "*********************************"
cd $BIOPERF/Scripts/$arch2-scripts/$i/clustalw
starttime=`date +%s`

	  ./clustalw$input2.script
endtime=`date +%s`
	  if [ $? -ne 0 ] ; then
	      echo $i did not run successfully 
exit
	  fi
 clustalwtotal=`expr $endtime - $starttime`
	  cd $BIOPERF/Scripts/$arch2-scripts/$i/clustalw-smp
starttime=`date +%s`
	  ./clustalw-smp$input2.script
endtime=`date +%s`
	  if [ $? -ne 0 ] ; then
	      echo $i did not run successfully                                                                                                          
	      exit
	  fi
 clustalwsmptotal=`expr $endtime - $starttime`
	fi  
	  
	  ;; 
"Fasta")   if [ $FASTA -eq "1" ] ; then 
 echo "*********************************"
	      echo Running $i
	      echo "*********************************"
 cd $BIOPERF/Scripts/$arch2-scripts/$i/fasta
starttime=`date +%s`	
  ./fasta-script
	endtime=`date +%s`  
	  if [ $? -ne 0 ] ; then
echo $i did not run successfully                                                                                                          
exit
	  fi
 fastatotal=`expr $endtime - $starttime`
	  cd $BIOPERF/Scripts/$arch2-scripts/$i/ssearch
starttime=`date +%s`	 
 ./ssearch$input2.script
	endtime=`date +%s`  
	  if [ $? -ne 0 ] ; then
	      echo $i did not run successfully                                                                                                          
exit
	  fi
	   ssearchtotal=`expr $endtime - $starttime`
fi
	  ;;
"Glimmer")  
 if [ $GLIMMER -eq "1" ] ; then 
 echo "*********************************"
	      echo Running $i
	      echo "*********************************"
 cd $BIOPERF/Scripts/$arch2-scripts/$i/glimmer
starttime=`date +%s`
	  ./glimmer$input2.script
endtime=`date +%s`
if [ $? -ne 0 ] ; then
echo $i did not run successfully                                                                                                          
exit

fi
 glimmertotal=`expr $endtime - $starttime`

           cd $BIOPERF/Scripts/$arch2-scripts/$i/glimmer-package
       starttime=`date +%s`
     ./glimmer-package$input2.script
endtime=`date +%s`
if [ $? -ne 0 ] ; then
echo $i did not run successfully                                                                                                          
exit
fi
 glimmerpackagetotal=`expr $endtime - $starttime`
fi

;;
"Grappa")  if [ $GRAPPA -eq "1" ] ; then 
 echo "*********************************"
	      echo Running $i
	      echo "*********************************"
cd $BIOPERF/Scripts/$arch2-scripts/$i/
         starttime=`date +%s`
  ./grappa$input2.script
endtime=`date +%s`

if [ $? -ne 0 ] ; then
echo $i did not run successfully                                                                                                          
exit
fi
 grappatotal=`expr $endtime - $starttime`
fi

;;
"Hmmer")  if [ $HMMER -eq "1" ] ; then 
 echo "*********************************"
	      echo Running $i
	      echo "*********************************"
cd $BIOPERF/Scripts/$arch2-scripts/$i/hmmpfam
        starttime=`date +%s`
   ./hmmpfam$input2.script
endtime=`date +%s`
if [ $? -ne 0 ] ; then
echo $i did not run successfully                                                                                                          
exit
fi
 hmmpfamtotal=`expr $endtime - $starttime`

         cd $BIOPERF/Scripts/$arch2-scripts/$i/hmmsearch
 starttime=`date +%s`
          ./hmmsearch-script
endtime=`date +%s`
if [ $? -ne 0 ] ; then
echo $i did not run successfully                                                                                                          
exit
fi
 hmmsearchtotal=`expr $endtime - $starttime`
fi

;;
"Phylip")  if [ $PHYLIP -eq "1" ] ; then 
 echo "*********************************"
	      echo Running $i
	      echo "*********************************"
cd $BIOPERF/Scripts/$arch2-scripts/$i/dnapenny
    starttime=`date +%s`
       ./dnapenny-script
endtime=`date +%s`
if [ $? -ne 0 ] ; then
echo $i did not run successfully                                                                                                          
exit
fi
 dnapennytotal=`expr $endtime - $starttime`



         cd $BIOPERF/Scripts/$arch2-scripts/$i/promlk
      starttime=`date +%s`
     ./promlk$input2.script
endtime=`date +%s`
if [ $? -ne 0 ] ; then
echo $i did not run successfully                                                                                                          
exit
fi
 promlktotal=`expr $endtime - $starttime`
fi

;;
"Predator")  if [ $PREDATOR -eq "1" ] ; then 
 echo "*********************************"
	      echo Running $i
	      echo "*********************************"
cd $BIOPERF/Scripts/$arch2-scripts/$i/
       starttime=`date +%s`
     ./predator$input2.script
endtime=`date +%s`
if [ $? -ne 0 ] ; then
echo $i did not run successfully                                                                                                          
exit
fi
 predatortotal=`expr $endtime - $starttime`
fi
;;
"Tcoffee")  if [ $TCOFFEE -eq "1" ] ; then 
 echo "*********************************"
	      echo Running $i
	      echo "*********************************"
cd $BIOPERF/Scripts/$arch2-scripts/$i/
            ./tcoffee$input2.script
endtime=`date +%s`
if [ $? -ne 0 ] ; then
echo $i did not run successfully                                                                                                          
exit
fi
 tcoffeetotal=`expr $endtime - $starttime`
fi
;;
esac
done

total=0
for i in Blast Ce Clustalw Fasta Glimmer Grappa Hmmer Phylip Predator Tcoffee

  do 
  case $i in

"Blast") 
if [ $BLAST -eq "1" ] ; then 
echo The time taken for running blastn for $input2 size input is $blastntotal
echo The time taken for running blastp for $input2 size input is $blastptotal
echo $'\n'
total=`expr $total + $blastntotal + $blastptotal`
fi;;

"Ce") 
if [ $CE -eq "1" ] ; then 
echo The time taken for running $i for $input2 size input is $cetotal
echo $'\n'
total=`expr $total + $cetotal`
fi;;



"Clustalw") 
if [ $CLUSTALW -eq "1" ] ; then 
echo The time taken for running clustalw for $input2 size input is $clustalwtotal
echo The time taken for running clustalw-smp for $input2 size input is $clustalwsmptotal
echo $'\n'
total=`expr $total + $clustalwtotal + $clustalwsmptotal`
fi;;

"Fasta") 
if [ $FASTA -eq "1" ] ; then 
echo The time taken for running fasta for $input2 size input is $fastatotal
echo The time taken for running ssearch for $input2 size input is $ssearchtotal

echo $'\n'
total=`expr $total + $fastatotal + $ssearchtotal`
fi;;

"Glimmer") 
if [ $HMMER -eq "1" ] ; then 
echo The time taken for running glimmer for $input2 size input is $glimmertotal
echo The time taken for running glimmer-package for $input2 size input is $glimmerpackagetotal
echo $'\n'
total=`expr $total + $glimmertotal + $glimmerpackagetotal`
fi;;
"Grappa") 
if [ $GRAPPA -eq "1" ] ; then 
echo The time taken for running $i for $input2 size input is $grappatotal
echo $'\n'
total=`expr $total + $grappatotal`
fi;;


"Hmmer") 
if [ $HMMER -eq "1" ] ; then 
echo The time taken for running  hmmpfam for $input2 size input is $hmmpfamtotal
echo The time taken for running hmmsearch for $input2 size input is $hmmsearchtotal
echo $'\n'
total=`expr $total + $hmmpfamtotal + $hmmsearchtotal`
fi;;
"Phylip") 
if [ $PHYLIP -eq "1" ] ; then 
echo The time taken for running dnapenny for $input2 size input is $dnapennytotal
echo The time taken for running promlk for $input2 size input is $promlktotal
echo $'\n'
total=`expr $total + $promlktotal + $dnapennytotal`
fi;;
"Predator") 
if [ $PREDATOR -eq "1" ] ; then 
echo The time taken for running $i for $input2 size input is $predatortotal
echo $'\n'
total=`expr $total + $predatortotal`
fi;;
"Tcoffee")
if [ $TCOFFEE -eq "1" ] ; then 
echo The time taken for running $i for $input2 size input is $tcoffeetotal
echo $'\n'
total=`expr $total + $tcoffeetotal`
fi;;
esac
done


echo  "******************************************************"
echo Total time for running BioPerf is $total seconds
echo  "******************************************************"
