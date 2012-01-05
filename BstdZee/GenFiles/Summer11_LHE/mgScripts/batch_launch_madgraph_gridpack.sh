#!/bin/bash

set -o verbose

echo "Starting job on " `date`
echo "Running on " `uname -a`
echo "System release " `cat /etc/redhat-release`


####################################
# 0. Sorting out preliminaries - environment vars etc.
# N.B: Two environment variables must be set before running script
#         a) 'name' - The name of the gridpack that will be used (w/o _grid.tar.gz)
#         b) 'rnum' - The random number seed for event generation

unset PBS_O_WORKDIR
echo " name = " ${name}
echo " rnum = " ${rnum}

#export name=${1} #Put the name corresponding to the needed gridpack in the official repository (without _grid.tar.gz)
#export rnum=$2
#export name=W4jets
export RELEASE=CMSSW_4_1_6
export STOREHOME=/opt/ppd/newscratch/williams/Summer11MG/centralGridPack_DYJetsToLL_incl/batchTests/2012-01-05/${name}
export nevt=10000

#define the random number generator seed from lumi section for this run
#rnum=`grep 'JobID="'${jobid}'"' $RUNTIME_AREA/arguments.xml | awk -F "\"" '{print $4}' | head -1 `

#setup CMS env
. /gridsoft/SL5/cms/cmsset_default.sh

export JOB_SCRATCH_DIR=/scratch/${PBS_JOBID} 
mkdir $JOB_SCRATCH_DIR


####################################
# 1. Retrieve the wanted gridpack from the official/local repository 
cd $JOB_SCRATCH_DIR
cp /afs/cern.ch/cms/generators/www/slc5_ia32_gcc434/madgraph/V5_1.1/7TeV_Summer11/Gridpacks/${name}_gridpack.tar.gz ./

export SOURCE=`pwd`/${name}_gridpack.tar.gz


####################################
# 2. Setup the working area for event generation
scram project -n ${name}_${rnum} CMSSW ${RELEASE} 
cd ${name}_${rnum} 
mkdir -p work  
cd work
eval `scram runtime -sh`


# force the f77 compiler to be the CMS defined one
ln -s `which gfortran` f77
ln -s `which gfortran` g77
export PATH=`pwd`:${PATH}

cp ${SOURCE} . 
tar xzf ${name}_gridpack.tar.gz  
rm -f ${name}_gridpack.tar.gz 
cd madevent

# run the production stage
./bin/compile
./bin/clean4grid

# Cleaning and recompilation of the CERNLIB for possible left-over in Madgraph5
rm -f lib/libcernlib.a
rm -f Source/CERNLIB/*.o
cd Source/CERNLIB/
make
cd ../..

####################################
# 3. Now actually generate the events
cd ..
./run.sh ${nevt} ${rnum}

####################################
# 4. Un-gzip the LHE file, and copy to the storage directory 
gunzip events.lhe.gz
rfmkdir ${STOREHOME}
rfcp events.lhe ${STOREHOME}/events_${name}_nevt${nevt}_seed${rnum}.lhe

rm -r ${JOB_SCRATCH_DIR}/*

echo "End of job on " `date`
exit 0

