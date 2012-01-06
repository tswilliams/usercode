#!/bin/bash

set -o verbose

echo "Starting job on " `date`

#LHEtoGENSIM_inFile=/opt/ppd/newscratch/williams/Summer11MG/centralGridPack_DYJetsToLL_incl/batchTests/CMSSW_4_1_7_patch1/src/events_zjets_smzerobmass_nevt5000_seed1_qcut10_mgPostv2.lhe
#LHEtoGENSIM_outFile=/opt/ppd/newscratch/williams/Summer11MG/centralGridPack_DYJetsToLL_incl/batchTests/CMSSW_4_1_7_patch1/src/events_zjets_smzerobmass_nevt5000_seed1_qcut10_mgPostv2_batch2.root

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
cmsDriver.py MCDBtoEDM --conditions START311_V2::All -s NONE --eventcontent RAWSIM --datatier GEN --filein=file:inputLHEfile.lhe -n -1 >& logFile_LHEtoGEN.txt

###################################################################
# 3) Now doing GEN -> GEN-SIM conversion (i.e. PYTHIA hadronisation + SIM)
cmsDriver.py Configuration/GenProduction/python/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph_tauola_cff.py \
  --step GEN,SIM --beamspot Realistic7TeV2011Collision \
  --conditions START311_V2::All \
  --pileup NoPileUp \
  --datamix NODATAMIXER \
  --eventcontent RAWSIM \
  --datatier GEN-SIM \
  --filein=file:MCDBtoEDM_NONE.root --fileout=outputGENSIM.root -n -1 >& logFile_GENtoGENSIM.txt

####################################################################
# And finally copy output files & clean up
mv outputGENSIM.root ${LHEtoGENSIM_outFile}
mv logFile_LHEtoGEN.txt ${LHEtoGENSIM_outFile}.log1
mv logFile_GENtoGENSIM.txt ${LHEtoGENSIM_outFile}.log2

rm -r ${JOB_SCRATCH_DIR}/*

echo "End of job on " `date`
exit 0
