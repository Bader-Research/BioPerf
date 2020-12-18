#!/bin/bash
source $HOME/.bioperf
echo This script deletes all output files in $BIOPERF/Outputs
echo "Do you want to proceed [y/n]" 
read "input"

if [ "$input" =  "y" ] ;  then 
#echo The outputs deleted
for i in Blast Ce Clustalw Fasta Glimmer Grappa Hmmer Phylip Predator Tcoffee
do 
rm -rf $BIOPERF/Outputs/$i/*
done 
echo Outputs Deleted
fi 
