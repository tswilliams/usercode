#!/bin/bash

sampleBaseDir=${3}
wrapperLogFile=${sampleBaseDir}/commandHistory.txt

echo >> $wrapperLogFile
echo ___________________________________________________________________________________ >> $wrapperLogFile
echo =================================================================================== >> $wrapperLogFile
echo Date/Time: `date` >> $wrapperLogFile
echo Command: $* >> $wrapperLogFile
echo "   ------------------------------    " >> $wrapperLogFile

$* 2>&1 | tee -a $wrapperLogFile

echo 
echo 
echo "    *** This output has been written to the file: "${wrapperLogFile}" ***"