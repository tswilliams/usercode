#
# Python script for copying multiple NTuple files from dcache to my /opt/ppd/scratch area ...
##

# Imports ...
import os

print '*** Python script moveFilesAndMerge.py is now running ...'

####################################
## Define the source directories, the dccp target directories, and the directories/filenames for the merged NTuple files ...
sourceDir = []
targetDir = []
mergedNTupleFname = []

#sourceDir.append("/pnfs/pp.rl.ac.uk/data/cms/store/user/tsw/WWTo2L2Nu_TuneZ2_7TeV_pythia6_tauola/NTuple-42Xv1b_WWTo2L2Nu_2011-07-17/c996fd5ccd1ef011a4c1a6d202aa2c27")
#targetDir.append("/opt/ppd/scratch/williams/Datafiles/NTuples/42Xv1b/2011-07-17/src/WWTo2L2Nu/")
#mergedNTupleFname.append("/opt/ppd/scratch/williams/Datafiles/NTuples/42Xv1b/2011-07-17/fullMC/WWTo2L2Nu_NTuple-42Xv1b_2011-07-17.root")

#sourceDir.append("/pnfs/pp.rl.ac.uk/data/cms/store/user/tsw/WZTo3LNu_TuneZ2_7TeV_pythia6_tauola/NTuple-42Xv1b_WZTo3LNu_2011-07-17/c996fd5ccd1ef011a4c1a6d202aa2c27")
#targetDir.append("/opt/ppd/scratch/williams/Datafiles/NTuples/42Xv1b/2011-07-17/src/WZTo3LNu/")
#mergedNTupleFname.append("/opt/ppd/scratch/williams/Datafiles/NTuples/42Xv1b/2011-07-17/fullMC/WZTo3LNu_NTuple-42Xv1b_2011-07-17.root")

#sourceDir.append("/pnfs/pp.rl.ac.uk/data/cms/store/user/tsw/ZZTo2L2Nu_TuneZ2_7TeV_pythia6_tauola/NTuple-42Xv1b_ZZTo2L2Nu_2011-07-17/c996fd5ccd1ef011a4c1a6d202aa2c27")
#targetDir.append("/opt/ppd/scratch/williams/Datafiles/NTuples/42Xv1b/2011-07-17/src/ZZTo2L2Nu/")
#mergedNTupleFname.append("/opt/ppd/scratch/williams/Datafiles/NTuples/42Xv1b/2011-07-17/fullMC/ZZTo2L2Nu_NTuple-42Xv1b_2011-07-17.root")

#sourceDir.append("/pnfs/pp.rl.ac.uk/data/cms/store/user/tsw/Photon/NTuple-42Xv1b_Photon_2011-07-20/c996fd5ccd1ef011a4c1a6d202aa2c27")
#targetDir.append("/opt/ppd/scratch/williams/Datafiles/NTuples/42Xv1b/2011-07-20/src/data-Photon/")
#mergedNTupleFname.append("/opt/ppd/scratch/williams/Datafiles/NTuples/42Xv1b/2011-07-20/data/data-Photon_0-5invfb_NTuple-42Xv1b_2011-07-20.root")

#sourceDir.append("/pnfs/pp.rl.ac.uk/data/cms/store/user/tsw/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/NTuple-42Xv1b_WJetsToLNu_2011-07-17/c996fd5ccd1ef011a4c1a6d202aa2c27")
#targetDir.append("/opt/ppd/scratch/williams/Datafiles/NTuples/42Xv1b/2011-07-17/src/WJetsToLNu/")
#mergedNTupleFname.append("/opt/ppd/scratch/williams/Datafiles/NTuples/42Xv1b/2011-07-17/fullMC/WJetsToLNu_NTuple-42Xv1b_2011-07-17.root")

#sourceDir.append("/pnfs/pp.rl.ac.uk/data/cms/store/user/tsw/QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6/NTuple-42Xv1b_QCD-dblEM_2011-07-17/c996fd5ccd1ef011a4c1a6d202aa2c27")
#targetDir.append("/opt/ppd/scratch/williams/Datafiles/NTuples/42Xv1b/2011-07-17/src/QCD-dblEM/")
#mergedNTupleFname.append("/opt/ppd/scratch/williams/Datafiles/NTuples/42Xv1b/2011-07-17/mc/QCD-dblEM_NTuple-42Xv1b_2011-08-03.root")
#
#sourceDir.append("/pnfs/pp.rl.ac.uk/data/cms/store/user/tsw/DYJetsToLL_PtZ-100_TuneZ2_7TeV-madgraph-tauola/NTuple-42Xv1b_DYJetsToLL-ZpT100_2011-08-01/c996fd5ccd1ef011a4c1a6d202aa2c27")
#targetDir.append("/opt/ppd/scratch/williams/Datafiles/NTuples/42Xv1b/2011-08-01/src/DYJetsToLL-ZpT100/")
#mergedNTupleFname.append("/opt/ppd/scratch/williams/Datafiles/NTuples/42Xv1b/2011-08-01/mc/DYJetsToLL-ZpT100_NTuple-42Xv1b_2011-08-03.root")
#
#sourceDir.append("/pnfs/pp.rl.ac.uk/data/cms/store/user/tsw/DYJetsToLL_PtZ-100_TuneZ2_7TeV-madgraph-tauola/NTuple-42Xv1b_DYJetsToLL-ZpT100-ee_2011-08-01/c996fd5ccd1ef011a4c1a6d202aa2c27")
#targetDir.append("/opt/ppd/scratch/williams/Datafiles/NTuples/42Xv1b/2011-08-01/src/DYJetsToLL-ZpT100-ee/")
#mergedNTupleFname.append("/opt/ppd/scratch/williams/Datafiles/NTuples/42Xv1b/2011-08-01/mc/DYJetsToLL-ZpT100-ee_NTuple-42Xv1b_2011-08-03.root")
#
#sourceDir.append("/pnfs/pp.rl.ac.uk/data/cms/store/user/tsw/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/NTuple-42Xv1b_DYJetsToLL-ee_2011-08-01/c996fd5ccd1ef011a4c1a6d202aa2c27")
#targetDir.append("/opt/ppd/scratch/williams/Datafiles/NTuples/42Xv1b/2011-08-01/src/DYJetsToLL-ee/")
#mergedNTupleFname.append("/opt/ppd/scratch/williams/Datafiles/NTuples/42Xv1b/2011-08-01/mc/DYJetsToLL-ee_NTuple-42Xv1b_2011-08-03.root")
#
#sourceDir.append("/pnfs/pp.rl.ac.uk/data/cms/store/user/tsw/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/NTuple-42Xv1b_DYJetsToLL-TauTau_2011-08-01/c996fd5ccd1ef011a4c1a6d202aa2c27")
#targetDir.append("/opt/ppd/scratch/williams/Datafiles/NTuples/42Xv1b/2011-08-01/src/DYJetsToLL-TauTau/")
#mergedNTupleFname.append("/opt/ppd/scratch/williams/Datafiles/NTuples/42Xv1b/2011-08-01/mc/DYJetsToLL-TauTau_NTuple-42Xv1b_2011-08-03.root")

sourceDir.append("/pnfs/pp.rl.ac.uk/data/cms/store/user/tsw/MuEG/NTuple-42Xv1b_MuEG_2011-07-20/c996fd5ccd1ef011a4c1a6d202aa2c27")
targetDir.append("/opt/ppd/scratch/williams/Datafiles/NTuples/42Xv1b/2011-07-20/src/data-MuEG/")
mergedNTupleFname.append("/opt/ppd/scratch/williams/Datafiles/NTuples/42Xv1b/2011-07-20/data/data-MuEG_0-5invfb_NTuple-42Xv1b_2011-08-05.root")

########################################
# Do the work for each dataset ...
for datasetIdx in range(len(sourceDir)):
   print " "
   print " * ---------------------------------------------- *"
   print " * New dataset ..."
   print "  dcache dir:    " + sourceDir[datasetIdx]
   print "  dccp target:   " + targetDir[datasetIdx]
   print "  merged NTuple: " + mergedNTupleFname[datasetIdx]
   # Copy the job output files into the src directory in /opt/ppd/scratch area
   listOfFiles = os.listdir(sourceDir[datasetIdx])
   counter = 0
   for fname in listOfFiles:
      if counter==30:
         break
      os.system('dccp ' + sourceDir[datasetIdx]+'/'+fname + ' ' + targetDir[datasetIdx])
      counter+=1
   
   # Merge the first 5 of the source files in the src directory into the fullMC directory ROOT file ...
   listOfSrcFiles = os.listdir(targetDir[datasetIdx])
   idx=0
   filesToMerge_str = ''
   while (idx<30 and idx<len(listOfSrcFiles)):
      filesToMerge_str = filesToMerge_str + ' ' + targetDir[datasetIdx] + listOfSrcFiles[idx]
      idx+=1
   os.system('hadd ' + mergedNTupleFname[datasetIdx] + filesToMerge_str)
