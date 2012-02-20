#!/bin/bash

##########################################
# This is a bash shell script used for automating the hadronisation of 
# LHE files during private production of high pT DY(ee)+Jets samples.
# 
# Author: T. Williams
# Date: Feb 2012
###########################################

qsub_sh=/opt/ppd/newscratch/williams/Summer11MG/Summer11_LHE/mgScripts/qsub_wrapper.sh

sampleName=$1
sampleBaseDir=$2
inputFileName=$3
nEvtsPerLHE=$4
nEvtsPerJob=$5

inputFileDir=${sampleBaseDir}/lhe
targetDir=${sampleBaseDir}/GEN_SIM

echo runFall11LHE_hadroniseLHE.sh script being run with following arguments:
echo "   Arg 1 (sample name):      " $1 
echo "   Arg 2 (sample store dir): " $2
echo "   Arg 3 (input file name):  " $3
echo "   Arg 4 (nEvts per LHE):    " $4
echo "   Arg 5 (num evts per job): " $5
echo

##############
# 0. Check that directories are present
if [ ! -d "$sampleBaseDir" ]; then
   echo "Directory '" $sampleBaseDir "' does not exist"
   echo "Exiting script early"
   exit 1
fi

if [ ! -d ${inputFileDir} ]; then
   echo "Directory '" $inputFileDir "' does not exist"
   echo "Exiting script early"
   exit 1
fi

if [ ! -d ${targetDir} ]; then 
   echo "Target directory '" $targetDir "' does not exist ..."
   echo "Will now create this directory ..."
   mkdir $targetDir
   mkdir ${targetDir}/jobLogs
   sleep 20
fi

##############
# 1. Submit the batch jobs for specified input file
jobIdsFile=${targetDir}/jobids_`date "+%y%m%d_%H%M%S"`.txt
echo Job ids file is: $jobIdsFile  

lheFileName=${inputFileDir}/${inputFileName}

# Checking that input file exists  
if [ ! -e ${lheFileName} ]; then
   echo "  ** WARNING: The LHE file " ${lheFileName} " doesn't exist. "
   echo "          Skipping submission of hadronisation jobs for this file. "
   exit    
fi

# Adding comment line to job id file  
echo "#Hadronise:: "${sampleName}", "${inputFileName}", "${nEvtsPerLHE}", "${nEvtsPerJob} >> $jobIdsFile

# Setting i=1, since i used to be used for running over multiple input files ...
i=1
# Splitting up processing of each LHE file into a number of GEN-SIM jobs (to keep job times down to acceptable levels)
nEvtsSkipped=0
for ((j=0; $(($j * $nEvtsPerJob))<$nEvtsPerLHE; j++)); do
   nEvtsSkipped=$(($nEvtsPerJob * $j))
   gensimFileName=${targetDir}/7TeV_${sampleName}_run${i}p${j}_GENSIM.root
   jobLogFileName=${targetDir}/jobLogs/7TeV_${sampleName}_run${i}p${j}_GENSIM.log
   echo "Job no. " ${i}"."${j}" being submitted ..."
   echo "       InputFile: " $lheFileName
   echo "       OutputFile: " $gensimFileName
   echo "       nEvts Skipped: " $nEvtsSkipped " , nEvts Processed: " $nEvtsPerJob
   echo "#==> "${gensimFileName} >> $jobIdsFile
   $qsub_sh -q prod  -j oe -o ${jobLogFileName} -l cput=72:00:00,walltime=96:00:00,mem=1gb -N ${sampleName}_GENSIM_${i}p${j} \
        -v LHEtoGENSIM_inFile=${lheFileName},LHEtoGENSIM_outFile=${gensimFileName},LHEtoGENSIM_runNum=${i},LHEtoGENSIM_numEvtsProc=${nEvtsPerJob},LHEtoGENSIM_numEvtsSkip=${nEvtsSkipped} \
         batch_LHEtoGENSIM.sh >> ${jobIdsFile}
done

echo "Have reached the end of runFall11LHE_hadronise.sh script."