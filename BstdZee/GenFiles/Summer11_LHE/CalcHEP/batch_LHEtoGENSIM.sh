#!/bin/bash

set -o verbose

echo "Starting job on " `date`

echo
echo Env vars passed into this job are:
echo "    LHEtoGENSIM_inFile      = " $LHEtoGENSIM_inFile
echo "    LHEtoGENSIM_outFile     = " $LHEtoGENSIM_outFile
echo "    LHEtoGENSIM_runNum      = " $LHEtoGENSIM_runNum
echo "    LHEtoGENSIM_numEvtsSkip = " $LHEtoGENSIM_numEvtsSkip
echo "    LHEtoGENSIM_numEvtsProc = " $LHEtoGENSIM_numEvtsProc
echo 

###################################################################
# 1) Setting up the environment on the node ...

export JOB_SCRATCH_DIR=/scratch/${PBS_JOBID}
mkdir $JOB_SCRATCH_DIR
cd $JOB_SCRATCH_DIR

. /gridsoft/SL5/cms/cmsset_default.sh # To setup CMS env vars
export SCRAM_ARCH=slc5_amd64_gcc434
scram project CMSSW CMSSW_4_1_7_patch1
cd CMSSW_4_1_7_patch1/src
eval `scram runtime -sh`
export CVSROOT=:pserver:anonymous@cmscvs.cern.ch:/cvs/CMSSW
cvs -d :pserver:anonymous:98passwd@cmscvs.cern.ch:/cvs/CMSSW login
cvs co -r V01-00-86 Configuration/GenProduction 
cvs co -r V01-00-86 Configuration/Generator
cvs logout
scramv1 b --j 4

###################################################################
# 2) Copying over input LHE files, and doing LHE -> GEN conversion
cp ${LHEtoGENSIM_inFile} inputLHEfile.lhe
cmsDriver.py MCDBtoEDM --conditions START311_V2::All -s NONE --eventcontent RAWSIM --datatier GEN --filein=file:inputLHEfile.lhe -n -1 --no_exec --python_filename=orig_LHEtoGEN_cfg.py
sed 's/"LHESource",/"LHESource", firstRun=cms.untracked.uint32\('${LHEtoGENSIM_runNum}'\), /g' orig_LHEtoGEN_cfg.py > LHEtoGEN_cfg.py
grep "LHESource" LHEtoGEN_cfg.py
cmsRun LHEtoGEN_cfg.py >& logFile_LHEtoGEN.txt

###################################################################
# 3) Now doing GEN -> GEN-SIM conversion (i.e. PYTHIA hadronisation + SIM)
cmsDriver.py Configuration/GenProduction/python/Hadronizer_TuneZ2_7TeV_madgraph_cff.py \
  --step GEN --beamspot Realistic7TeV2011Collision \
  --conditions START311_V2::All \
  --pileup NoPileUp \
  --datamix NODATAMIXER \
  --eventcontent RAWSIM \
  --datatier GEN-SIM \
  --filein=file:MCDBtoEDM_NONE.root --fileout=outputGENSIM.root -n ${LHEtoGENSIM_numEvtsProc} --no_exec --python_filename=orig_GENtoGENSIM_cfg.py

sed 's/"PoolSource",/"PoolSource", skipEvents=cms.untracked.uint32\('${LHEtoGENSIM_numEvtsSkip}'\), /g' orig_GENtoGENSIM_cfg.py > GENtoGENSIM_cfg.py
grep "PoolSource" GENtoGENSIM_cfg.py
cmsRun GENtoGENSIM_cfg.py >& logFile_GENtoGENSIM.txt
####################################################################
# And finally copy output files & clean up
mv outputGENSIM.root ${LHEtoGENSIM_outFile}
mv logFile_LHEtoGEN.txt ${LHEtoGENSIM_outFile}.log1
mv logFile_GENtoGENSIM.txt ${LHEtoGENSIM_outFile}.log2

rm -r ${JOB_SCRATCH_DIR}/*

echo "End of job on " `date`
exit 0
