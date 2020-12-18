#!/bin/bash
echo Please specify the installation directory for BioPerf [$HOME/BioPerf]
read BIOPERF
#echo This is BIOPERF $BIOPERF
if [ -z $BIOPERF ]; then 
BIOPERF=$HOME/BioPerf
fi
echo "#!/bin/bash" $'\n'"export BIOPERF=$BIOPERF" $'\n'"if [ ! -d $BIOPERF ] ; then" $'\n'"echo The path to BioPerf is set incorrectly" $'\n'"echo Please set the correct path to your BioPerf installation in  $HOME/.bioperf" $'\n'"exit" $'\n'"fi" $'\n' > $HOME/.bioperf
if [ $? != 0 ] ; then
echo Unable to set the environement variable 
exit
fi
echo The environment variable BIOPERF is set to $BIOPERF in $HOME/.bioperf.
echo Please do not change this file
if [ ! -d $BIOPERF ] ; then 
echo Making installation directory $BIOPERF
mkdir -p $BIOPERF
fi 
if [ $? != 0 ] ; then 
echo Unable to create Installation directory $BIOPERF 
echo Installation aborted
exit 
fi
echo Installation proceeding in $BIOPERF
mv BioPerf-codes.tar.gz install-BioPerf.sh $BIOPERF
cd $BIOPERF
zcat BioPerf-codes.tar.gz | tar xvf -
if [ $? = 0 ] ; then 
echo Installation of all codes done in $BIOPERF/
$BIOPERF/Scripts/Run-scripts/display-versions.sh
echo $'\n' 
echo $'\n'
echo The path to BioPerf is set to $BIOPERF in $HOME/.bioperf
echo 'Always make sure the path to BioPerf is set correctly in' $HOME/.bioperf 
else 
echo Installation could not be completed in BIOPERF. 
fi
