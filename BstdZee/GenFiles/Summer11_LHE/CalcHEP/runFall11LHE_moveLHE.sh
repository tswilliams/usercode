#!/bin/bash

##########################################
# This is a bash shell script used for automation of moving already produced 
# LHE files during private production of high pT DY(ee)+Jets samples.
# 
# Author: T. Williams
# Date: Feb 2012
###########################################

sampleName=$1
sampleBaseDir=$2
origLHEDir=$3

lheDir=${sampleBaseDir}/lhe/

echo runFall11LHE_moveLHE.sh script being run with the following arguments:
echo "   Arg 1 (sample name):         " $1 
echo "   Arg 2 (sample store dir):    " $2
echo "   Arg 3 (store dir - orig lhe):" $3
echo 

#####################################
# 0. Check that output directory exists

if [ ! -d "$sampleBaseDir" ]; then
   echo "Sample base directory does not exist. Exiting script now ..."
   exit 1
fi

if [ ! -d "$lheDir" ]; then
   echo "LHE directory does not exist. Creating it now ..."
   mkdir $lheDir
   sleep 15
fi

#####################################
# 1. Copy over the LHE file(s) 

cp -r $origLHEDir/*$sampleName* $lheDir/
gzip -d $lheDir/*.lhe.gz

echo "Have reached the end of runFall11LHE_moveLHE.sh script."