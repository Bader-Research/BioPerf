#!/bin/bash
source $HOME/.bioperf
echo This script will try to install BioPerf codes on your architecture.
echo Previous installation for your architecture with all scripts will be deleted [y/n]
read input 
if [ "$input" = "y" ] ; then 
if [ -d $BIOPERF/Binaries/$HOSTNAME-Binaries ] ; then 
rm -rf $BIOPERF/Binaries/$HOSTNAME-Binaries 
fi
if [ -d $BIOPERF/Scripts/$HOSTNAME-scripts ] ; then 
rm -rf $BIOPERF/Scripts/$HOSTNAME-scripts
fi
fi

if [ "$input" = "n" ] ; then 
echo No installation will be done 
exit 
fi

if [ ! -d $BIOPERF/Binaries/$HOSTNAME-Binaries ] ; then                   
mkdir $BIOPERF/Binaries/$HOSTNAME-Binaries
fi

if [ ! -d $BIOPERF/Scripts/$HOSTNAME-scripts ] ; then                                                                    
mkdir $BIOPERF/Scripts/$HOSTNAME-scripts                                                                                
fi  

for i in "Ce" "Clustalw" "Clustalw_smp" "Fasta" "Glimmer" "Grappa" "Hmmer" "Blast" "Phylip" "Predator" "Tcoffee"
#for i in "Fasta"
do
cd $BIOPERF/Scripts/Compile-scripts/
./compile-$i.script
if [ $? -ne 0 ] ; then
echo $'\n'
echo $i did not compile successfully on $HOSTNAME
echo Please check the installation again by going to directory $BIOPERF/Source-codes/$i
exit
fi
echo $'\n'
echo "*************************************"
echo $'\n'
echo $i successfully compiled on $HOSTNAME
echo $'\n'
#if [ ! -d $BIOPERF/Binaries/$HOSTNAME-Binaries/$i ] ;  then
#mkdir $BIOPERF/Binaries/$HOSTNAME-Binaries/$i                                                               
#fi  

#if [ ! -d $BIOPERF/Scripts/$HOSTNAME-scripts/$i ] ;  then                                                         
#mkdir $BIOPERF/Scripts/$HOSTNAME-scripts/$i                                                                             
#fi  

case $i in
"Ce") 
	if [ ! -d $BIOPERF/Binaries/$HOSTNAME-Binaries/$i ] ;  then
	    mkdir $BIOPERF/Binaries/$HOSTNAME-Binaries/$i                                                               
	fi  
	
	if [ ! -d $BIOPERF/Scripts/$HOSTNAME-scripts/$i ] ;  then                                                         
	    mkdir $BIOPERF/Scripts/$HOSTNAME-scripts/$i                                                                             
	fi  
	
	cd $BIOPERF/Source-codes/$i/
	cp -r CE $BIOPERF/Source-codes/$i/pom $BIOPERF/Binaries/$HOSTNAME-Binaries/$i/
	cd $BIOPERF/Scripts/x86-scripts/$i
	for j in *
	  do  
	  echo $j
	  sed -e 's!x86!$HOSTNAME!g' $j > $BIOPERF/Scripts/$HOSTNAME-scripts/$i/$j
	done 
	echo "";;
    
"Clustalw") 
	if [ ! -d $BIOPERF/Binaries/$HOSTNAME-Binaries/$i ] ;  then
	    mkdir $BIOPERF/Binaries/$HOSTNAME-Binaries/$i                                                               
	fi  
	
	if [ ! -d $BIOPERF/Scripts/$HOSTNAME-scripts/$i ] ;  then                                                         
	    mkdir $BIOPERF/Scripts/$HOSTNAME-scripts/$i                                                                             
	fi  
	
	cd $BIOPERF/Source-codes/$i/
	cp clustalw $BIOPERF/Binaries/$HOSTNAME-Binaries/$i/
	if [ ! -d $BIOPERF/Scripts/$HOSTNAME-scripts/$i/clustalw ] ;  then
	    mkdir $BIOPERF/Scripts/$HOSTNAME-scripts/$i/clustalw
	fi
	
	cd $BIOPERF/Scripts/x86-scripts/$i/clustalw
	for j in *                                                                       
	  do                                                                                                                
	  sed -e 's!x86!$HOSTNAME!g' $j > $BIOPERF/Scripts/$HOSTNAME-scripts/$i/clustalw/$j                                          
	done 
	echo "";;   
    
"Clustalw_smp")
	if [ ! -d $BIOPERF/Binaries/$HOSTNAME-Binaries/Clustalw ] ;  then
	    mkdir $BIOPERF/Binaries/$HOSTNAME-Binaries/Clustalw                                                               
	fi  
	
	if [ ! -d $BIOPERF/Scripts/$HOSTNAME-scripts/Clustalw ] ;  then                                                         
	    mkdir $BIOPERF/Scripts/$HOSTNAME-scripts/Clustalw                                                                              
	fi 
	if [ ! -d $BIOPERF/Scripts/$HOSTNAME-scripts/Clustalw/clustalw-smp ] ;  then 
	    mkdir $BIOPERF/Scripts/$HOSTNAME-scripts/Clustalw/clustalw-smp                                                                             
	fi 
	cd $BIOPERF/Source-codes/$i/
#	if [ ! -d $BIOPERF/Binaries/$HOSTNAME-Binaries/Clustalw/clustalw-smp ] ;  then
	#    mkdir $BIOPERF/Binaries/$HOSTNAME-Binaries/Clustalw/clustalw-smp
	#fi
	    cp clustalw $BIOPERF/Binaries/$HOSTNAME-Binaries/Clustalw/clustalw-smp
	    if [ ! -d $BIOPERF/Scripts/$HOSTNAME-scripts/Clustalw/clustalw-smp ] ;  then
		mkdir $BIOPERF/Scripts/$HOSTNAME-scripts/Clustalw/clustalw-smp
	    fi
	    cd $BIOPERF/Scripts/x86-scripts/Clustalw/clustalw-smp
	    for j in *
	      do
	      sed -e 's!x86!$HOSTNAME!g' $j > $BIOPERF/Scripts/$HOSTNAME-scripts/Clustalw/clustalw-smp/$j
	    done
echo "";;



"Fasta")
	if [ ! -d $BIOPERF/Binaries/$HOSTNAME-Binaries/$i ] ;  then
	    mkdir $BIOPERF/Binaries/$HOSTNAME-Binaries/$i                                                               
	fi  
	
	if [ ! -d $BIOPERF/Scripts/$HOSTNAME-scripts/$i ] ;  then                                                         
	    mkdir $BIOPERF/Scripts/$HOSTNAME-scripts/$i                                                                             
	fi  
	
	cd  $BIOPERF/Source-codes/$i/
	cp ssearch34_t fasta34_t $BIOPERF/Binaries/$HOSTNAME-Binaries/$i/
	if [ ! -d $BIOPERF/Scripts/$HOSTNAME-scripts/$i/fasta ] ;  then
	    mkdir $BIOPERF/Scripts/$HOSTNAME-scripts/$i/fasta
	fi
	
	
	cd $BIOPERF/Scripts/x86-scripts/$i/fasta/
	for j in *
	  do
	  sed -e 's!x86!$HOSTNAME!g' $j > $BIOPERF/Scripts/$HOSTNAME-scripts/$i/fasta/$j
	done
	if [ ! -d $BIOPERF/Scripts/$HOSTNAME-scripts/$i/ssearch ] ;  then
	    mkdir $BIOPERF/Scripts/$HOSTNAME-scripts/$i/ssearch
	fi
	
	
	cd $BIOPERF/Scripts/x86-scripts/$i/ssearch/
	for j in *                                                                                                                         
	  do                                                                                                                                 
	  sed -e 's!x86!$HOSTNAME!g' $j > $BIOPERF/Scripts/$HOSTNAME-scripts/$i/ssearch/$j
	done  
	echo "";;
    
    
"Glimmer") 
	if [ ! -d $BIOPERF/Binaries/$HOSTNAME-Binaries/$i ] ;  then
	    mkdir $BIOPERF/Binaries/$HOSTNAME-Binaries/$i                                                               
	fi  
	
	if [ ! -d $BIOPERF/Scripts/$HOSTNAME-scripts/$i ] ;  then                                                         
	    mkdir $BIOPERF/Scripts/$HOSTNAME-scripts/$i                                                                             
	fi  
	cd $BIOPERF/Source-codes/$i
        build-icm<glimmer.seq>glimmer.icm
	cp build-icm extract get-putative glimmer2 long-orfs run-glimmer2 glimmer.icm $BIOPERF/Binaries/$HOSTNAME-Binaries/$i/
	if [ ! -d $BIOPERF/Scripts/$HOSTNAME-scripts/$i/glimmer ] ;  then
	    mkdir $BIOPERF/Scripts/$HOSTNAME-scripts/$i/glimmer
	fi
	cd $BIOPERF/Scripts/x86-scripts/$i/glimmer/
	for j in *
	  do
	  sed -e 's!x86!$HOSTNAME!g' $j > $BIOPERF/Scripts/$HOSTNAME-scripts/$i/glimmer/$j
	done
	if [ ! -d $BIOPERF/Scripts/$HOSTNAME-scripts/$i/glimmer-package ] ;  then
	    mkdir $BIOPERF/Scripts/$HOSTNAME-scripts/$i/glimmer-package
	fi
	cd $BIOPERF/Scripts/x86-scripts/$i/glimmer-package
	for j in *
	  do
	  sed -e 's!x86!$HOSTNAME!g' $j > $BIOPERF/Scripts/$HOSTNAME-scripts/$i/glimmer-package/$j
	done
	echo "";;





"Grappa")
	if [ ! -d $BIOPERF/Binaries/$HOSTNAME-Binaries/$i ] ;  then
	    mkdir $BIOPERF/Binaries/$HOSTNAME-Binaries/$i                                                               
	fi  
	
	if [ ! -d $BIOPERF/Scripts/$HOSTNAME-scripts/$i ] ;  then                                                         
	    mkdir $BIOPERF/Scripts/$HOSTNAME-scripts/$i                                                                             
	fi  
	
	cd $BIOPERF/Source-codes/$i
	cp grappa $BIOPERF/Binaries/$HOSTNAME-Binaries/$i
	cd $BIOPERF/Scripts/x86-scripts/$i
	for j in *
	  do
	  sed -e 's!x86!$HOSTNAME!g' $j > $BIOPERF/Scripts/$HOSTNAME-scripts/$i/$j
	done
	echo "";;

"Hmmer")
	cd $BIOPERF/Source-codes/$i/src
	cp hmmpfam hmmsearch $BIOPERF/Binaries/$HOSTNAME-Binaries/Hmmer
	if [ ! -d $BIOPERF/Scripts/$HOSTNAME-scripts/$i/hmmpfam ] ;  then
	    mkdir $BIOPERF/Scripts/$HOSTNAME-scripts/$i/hmmpfam
	fi
	cd $BIOPERF/Scripts/x86-scripts/$i/hmmpfam
	for j in *
	  do
	  sed -e 's!x86!$HOSTNAME!g' $j > $BIOPERF/Scripts/$HOSTNAME-scripts/$i/hmmpfam/$j
	done
	if [ ! -d $BIOPERF/Scripts/$HOSTNAME-scripts/$i/hmmsearch ] ;  then
	    mkdir $BIOPERF/Scripts/$HOSTNAME-scripts/$i/hmmsearch
	fi
	cd $BIOPERF/Scripts/x86-scripts/$i/hmmsearch
	for j in *
	  do
	  sed -e 's!x86!$HOSTNAME!g' $j > $BIOPERF/Scripts/$HOSTNAME-scripts/$i/hmmsearch/$j
	done
	echo "";;





"Blast")
	if [ ! -d $BIOPERF/Binaries/$HOSTNAME-Binaries/$i ] ;  then
	    mkdir $BIOPERF/Binaries/$HOSTNAME-Binaries/$i                                                               
	fi  
	
	if [ ! -d $BIOPERF/Scripts/$HOSTNAME-scripts/$i ] ;  then                                                         
	    mkdir $BIOPERF/Scripts/$HOSTNAME-scripts/$i                                                                             
	fi  
	cd $BIOPERF/Source-codes/ncbi
	cp BLOSUM62 build/blastall $BIOPERF/Binaries/$HOSTNAME-Binaries/Blast/
	if [ ! -d $BIOPERF/Scripts/$HOSTNAME-scripts/$i/blastn ] ;  then
	    mkdir $BIOPERF/Scripts/$HOSTNAME-scripts/$i/blastn
	fi
	cd $BIOPERF/Scripts/x86-scripts/$i/blastn
	for j in *
	  do
	  sed -e 's!x86!$HOSTNAME!g' $j > $BIOPERF/Scripts/$HOSTNAME-scripts/$i/blastn/$j
	done
	if [ ! -d $BIOPERF/Scripts/$HOSTNAME-scripts/$i/blastp ] ;  then
	    mkdir $BIOPERF/Scripts/$HOSTNAME-scripts/$i/blastp
	fi
	for j in *
	  do
	  sed -e 's!x86!$HOSTNAME!g' $j > $BIOPERF/Scripts/$HOSTNAME-scripts/$i/blastp/$j
	done
	echo "";;
    
"Phylip") 
	if [ ! -d $BIOPERF/Binaries/$HOSTNAME-Binaries/$i ] ;  then
	    mkdir $BIOPERF/Binaries/$HOSTNAME-Binaries/$i                                                               
	fi  
	
	if [ ! -d $BIOPERF/Scripts/$HOSTNAME-scripts/$i ] ;  then                                                         
	    mkdir $BIOPERF/Scripts/$HOSTNAME-scripts/$i                                                                             
	fi  
	cd $BIOPERF/Source-codes/$i
	cp src/dnapenny src/promlk $BIOPERF/Binaries/$HOSTNAME-Binaries/$i
	if [ ! -d $BIOPERF/Scripts/$HOSTNAME-scripts/$i/dnapenny ] ;  then
	    mkdir $BIOPERF/Scripts/$HOSTNAME-scripts/$i/dnapenny
	fi
	cd $BIOPERF/Scripts/x86-scripts/$i/dnapenny
	for j in *
	  do
	  sed -e 's!x86!$HOSTNAME!g' $j > $BIOPERF/Scripts/$HOSTNAME-scripts/$i/dnapenny/$j
	done
	if [ ! -d $BIOPERF/Scripts/$HOSTNAME-scripts/$i/promlk ] ;  then
	    mkdir $BIOPERF/Scripts/$HOSTNAME-scripts/$i/promlk
	fi
	for j in *
	  do
	  sed -e 's!x86!$HOSTNAME!g' $j > $BIOPERF/Scripts/$HOSTNAME-scripts/$i/promlk/$j
	done 
	echo "";;
    
"Predator")
	if [ ! -d $BIOPERF/Binaries/$HOSTNAME-Binaries/$i ] ;  then
	    mkdir $BIOPERF/Binaries/$HOSTNAME-Binaries/$i                                                               
	fi  
	
	if [ ! -d $BIOPERF/Scripts/$HOSTNAME-scripts/$i ] ;  then                                                         
	    mkdir $BIOPERF/Scripts/$HOSTNAME-scripts/$i                                                                             
	fi  
	cd $BIOPERF/Source-codes/$i
	cp predator $BIOPERF/Binaries/$HOSTNAME-Binaries/$i
	cd $BIOPERF/Scripts/x86-scripts/$i
	for j in *
	  do
	  sed -e 's!x86!$HOSTNAME!g' $j > $BIOPERF/Scripts/$HOSTNAME-scripts/$i/$j
	done
	echo "";;
    
"Tcoffee")
	if [ ! -d $BIOPERF/Binaries/$HOSTNAME-Binaries/$i ] ;  then
	    mkdir $BIOPERF/Binaries/$HOSTNAME-Binaries/$i                                                               
	fi  
	
	if [ ! -d $BIOPERF/Scripts/$HOSTNAME-scripts/$i ] ;  then                                                         
	    mkdir $BIOPERF/Scripts/$HOSTNAME-scripts/$i                                                                             
	fi  
	cd $BIOPERF/Source-codes/$i
	cp bin/t_coffee bin/lalign2list $BIOPERF/Binaries/$HOSTNAME-Binaries/$i
	cd $BIOPERF/Scripts/x86-scripts/$i
	for j in *
	  do
	  sed -e 's!x86!$HOSTNAME!g' $j > $BIOPERF/Scripts/$HOSTNAME-scripts/$i/$j
	done
	echo "";;
    
esac
echo "*************************************"
echo $'\n'
done
chmod -R 755 $BIOPERF/Binaries/$HOSTNAME-Binaries/ $BIOPERF/Scripts/$HOSTNAME-scripts

