#!/bin/bash

######################################
# Purpose: Script that submits batch jobs that run BstdZeeFirstAnalyser over various different sets of input files
#             - Arguments are:
#                  $2 => ANR_OUTPUT_DIR  : Directory that output will be stored in.
#                  $1 => ANR_CONFIG_FILE : Location of config file that details of what to run over will be read in from 
#             - 
# Author: Tom S Williams
# Date: June 2012
######################################

if [ ! $# == 2 ]; then
   echo " ERROR : Incorrect number of aguments given for script $0 - there should be 2."
   echo "            * Arg 1 : Location of config file specifying what to run over"
   echo "            * Arg 2 : Base directory for output to be stored in"
   echo "         Script $0 now exiting early !"
   exit 1
fi

ANR_OUTPUT_DIR=$2
ANR_CONFIG_FILE=$1

echo 
echo "   ANR_OUTPUT_DIR  ="$ANR_OUTPUT_DIR
echo "   ANR_CONFIG_FILE ="$ANR_CONFIG_FILE

DATE_TIME_STRING=`date +%Y%m%d_%H%M%S`


######################################
# 0. Copy my CMSSW libs & Analyser binary to /opt/newscratch in order to freeze them, and adapt LD_LIBRARY_PATH accordingly
echo
echo "  ** 0: Copying over libs & Analyser to frozen dir"

if [ -z "$CMSSW_BASE" ]; then
   echo "  ERROR : Env var CMSSW_BASE is not set ; please issue cmsenv command before running this script."
   echo "          Script exiting early !" 
   exit 1
fi

THE_CMSSW_LIB_DIR=$CMSSW_BASE/lib/$SCRAM_ARCH
THE_CMSSW_BIN_DIR=$CMSSW_BASE/bin/$SCRAM_ARCH

FROZEN_BASE_DIR=/opt/ppd/newscratch/williams/BatchJobs/BstdZee_FrozenLibs/$DATE_TIME_STRING
FROZEN_LIB_DIR=$FROZEN_BASE_DIR/lib
FROZEN_BIN_DIR=$FROZEN_BASE_DIR/bin

## Copying over to freeze
mkdir -p $FROZEN_LIB_DIR
cp $THE_CMSSW_LIB_DIR/* $FROZEN_LIB_DIR
mkdir $FROZEN_BIN_DIR
cp $THE_CMSSW_BIN_DIR/BstdZeeFirstAnalyser $FROZEN_BIN_DIR

## Adapting LD_LIBRARY_PATH accordingly
LD_LIBRARY_PATH="${LD_LIBRARY_PATH/${THE_CMSSW_LIB_DIR}:/}"
LD_LIBRARY_PATH=$FROZEN_LIB_DIR:$LD_LIBRARY_PATH

echo "        Frozen lib dir is: "$FROZEN_LIB_DIR
echo "          ... containing:"
echo "`ls $FROZEN_LIB_DIR`"
echo

echo "        Frozen bin dir is: "$FROZEN_BIN_DIR
echo "          ... containing:"
echo `ls $FROZEN_BIN_DIR`
echo

echo "        And now, LD_LIBRARY_PATH is:"
echo $LD_LIBRARY_PATH
echo 

###########################################
# 1. Parse config file to produce BstdZeeFirstAnalyser runtime options with python script 
echo
echo "  ** 1: Parsing the config file:"

declare PARSED_CONFIG_FILE=$(scripts/parse_batch_config_file.py $ANR_CONFIG_FILE $ANR_OUTPUT_DIR)

echo
echo "$PARSED_CONFIG_FILE"

###########################################
# 2. Submit the batch jobs
echo 
echo 
echo "  ** 2: Submit the batch jobs"

mkdir -p $ANR_OUTPUT_DIR

QSUB_CMD=scripts/qsub.sh

## Run over each line in PARSED_CONFIG_FILE
echo "$PARSED_CONFIG_FILE" | while read line; do

   ## Split PARSED_CONFIG_FILE line into batch job params
   read THE_JOB_NAME THE_JOB_LOG_FILE THE_ANR_ARGS <<<$(IFS=";"; echo $line)

   export BSTDZEE_ANR_ARGS="$THE_ANR_ARGS"
   export BSTDZEE_ANR_OUT_DIR="$ANR_OUTPUT_DIR"
   
   echo 
   echo "  + About to submit job '"$THE_JOB_NAME"' with following parameters:"
   echo "       - jobLogFile ="$THE_JOB_LOG_FILE
   echo "       - anrOptions ="$THE_ANR_ARGS
   echo "       - outputDir  ="$BSTDZEE_ANR_OUT_DIR

   $QSUB_CMD -q prod -j oe -o$THE_JOB_LOG_FILE -l cput=12:00:00,walltime=12:00:00,mem=1gb -N $THE_JOB_NAME \
          -v LD_LIBRARY_PATH,BSTDZEE_ANR_ARGS,BSTDZEE_ANR_OUT_DIR,BSTDZEE_ANR_BIN_DIR=$FROZEN_BIN_DIR scripts/batch_run_analyser.sh >& __tmp_jobNo.txt
   echo "       * Job submitted ("`cat __tmp_jobNo.txt`")"
   rm __tmp_jobNo.txt

done

echo
echo "End of script reached. Goodbye!"
echo 
echo 
