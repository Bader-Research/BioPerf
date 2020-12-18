#!/bin/bash
source $HOME/.bioperf
echo You can do the following: 
echo '[R] Run BioPerf' 
echo '[I] Install BioPerf on  your architecture'
echo  ' (if your architecture is not PowerPc, x86)'
echo '[C] Clean outputs in '$BIOPERF/Outputs  
echo '[D] Display all versions of the installed codes'
read option
case $option in 
"C" ) $BIOPERF/Scripts/Run-scripts/CleanOutputs.sh;;
"D" ) $BIOPERF/Scripts/Run-scripts/display-versions.sh;;
"I" ) $BIOPERF/Scripts/Run-scripts/install-codes.sh;;
"R" ) $BIOPERF/Scripts/Run-scripts/run-bioperf.sh;;
esac
