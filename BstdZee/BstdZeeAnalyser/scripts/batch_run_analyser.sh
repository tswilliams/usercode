#!/bin/bash

########################################################
# Script that is run in analysis exe batch jobs
# Env vars passed in by qsub:
#       * LD_LIBARY_PATH      : Modified LD_LIBRARY_PATH that points to frozen version of my libs & boost/root libs etc
#       * BSTDZEE_ANR_BIN_DIR : Directory containing the BstdZeeAnalyser exe
#       * BSTDZEE_ANR_ARGS    : Arguments to use with BstdZeeAnalyser
#       * BSTDZEE_ANR_OUT_DIR : Location to copy output files to
# 
# Author: Tom S Williams
# Date:   June 2012
########################################################

echo "Entering batch_run_analyser.sh script ... "
echo 
echo "The values of the env vars that qsub should have passed to me are"
echo "  * LD_LIBARY_PATH      ="$LD_LIBRARY_PATH
echo "  * BSTDZEE_ANR_BIN_DIR ="$BSTDZEE_ANR_BIN_DIR
echo "  * BSTDZEE_ANR_ARGS    ="$BSTDZEE_ANR_ARGS
echo "  * BSTDZEE_ANR_OUT_DIR ="$BSTDZEE_ANR_OUT_DIR


# STAGE 0. cd into the data (i.e. large capacity) dir for this job
echo
JOB_SCRATCH_DIR=/scratch/${PBS_JOBID}
echo "About to cd into job's scratch dir: "$JOB_SCRATCH_DIR
mkdir -p $JOB_SCRATCH_DIR
cd $JOB_SCRATCH_DIR


# STAGE 1. Run the analyser
echo
echo " ================ Now running the analyser =================== "
echo

$BSTDZEE_ANR_BIN_DIR/BstdZeeFirstAnalyser  $BSTDZEE_ANR_ARGS

echo
echo " ================ FINSHED RUNNING ANALYSER =================== "
echo


# STAGE 2. Copy the analyser output to /opt/ppd/newscratch

echo "Current directory: "`pwd`
echo " * Contents:"
ls -lh
echo
echo "Copying ROOT files to "$BSTDZEE_ANR_OUT_DIR
cp *.root $BSTDZEE_ANR_OUT_DIR/


echo
echo
echo "... reached end of batch_run_analyser.sh script. Goodbye! "
echo

