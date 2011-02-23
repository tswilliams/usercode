
import os

#Setting which hlt menu config file will be run, and what data/mc datafiles it will be run on...
hltMenuConfigFile = "hlt_mc_11-02-05.py"
signalMCfnames = []
hltResultsFnames = []

#signalMCfnames.append("DataFiles/MC/BoostedZeeLFNs_0-75TeVexcitedQuark_11-02-08.txt")
#hltResultsFnames.append("hltResults/11-02-13/hltResults_0-75TeV_unifiedThres.log")

#signalMCfnames.append("DataFiles/MC/BoostedZeeLFNs_1-00TeVexcitedQuark_11-02-08.txt")
#hltResultsFnames.append("hltResults/11-02-13/hltResults_1-00TeV_unifiedThres.log")

#signalMCfnames.append("DataFiles/MC/BoostedZeeLFNs_2-00TeVexcitedQuark_11-02-08.txt")
#hltResultsFnames.append("hltResults/11-02-13/hltResults_2-00TeV_unifiedThres.log")

signalMCfnames.append("DataFiles/MC/BoostedZeeBkgdLFNs_SMZee_11-02-19.txt")
hltResultsFnames.append("hltResults/11-02-19/hltResults_std5e32onSMZee.log")

# Running the HLT menu config file over all of the datasets specified above, and summarising the standard output obtained...
for i in range(len(signalMCfnames)):
	#Running the HLT menu config file over the i'th dataset... 
   cmsRunCommand = 'cmsRun ' + hltMenuConfigFile + " " + signalMCfnames[i] + " >& " + hltResultsFnames[i]
   print "Running the command: " + cmsRunCommand
   os.system(cmsRunCommand)
   print "  ...done."
	
	#Opening up the hltResults file and writing summary information to the hltResultsSummary file...
   hltResultsSummFname = hltResultsFnames[i].replace("/hltResults","/hltResultsSummary")
   print "Writing summary information to", hltResultsSummFname
   hltResultsFile = open(hltResultsFnames[i],"r")
   hltResultsLines = hltResultsFile.readlines()
   hltResultsFile.close()
   hltResultsSummLines = []
   hltResultsSummLines.append("A summary of the HLT results resulting from the config file " + hltMenuConfigFile + " ...\n")
   hltResultsSummLines.append(" ... running over the first events from the datafiles list in " + signalMCfnames[i] + "\n")
   hltResultsSummLines.append("\n")
   copyLines = False
   for j in range(len(hltResultsLines)):
      if hltResultsLines[j].startswith("TrigReport ---------- Event  Summary ------------")==True:
         copyLines = True
      elif hltResultsLines[j].startswith("TrigReport -------End-Path   Summary ------------")==True:
         break
      if copyLines == True:
         hltResultsSummLines.append(hltResultsLines[j])

   hltResultsSummFile = open(hltResultsSummFname,"w")
   hltResultsSummFile.writelines(hltResultsSummLines)
   hltResultsSummFile.close()
   print "  ...done."
   
   #Opening up the hltResultsSummary file and writing just the electron and photon trigger path data to the hltResultsElePhotSummary file...
   hltResultsElePhotSummFname = hltResultsFnames[i].replace("/hltResults","/hltResultsElePhotSummary")
   print "Writing summary information for electron and photon triggers to", hltResultsElePhotSummFname
   hltResultsElePhotSummLines = []
   headerInfo = True
   for j2 in range(len(hltResultsSummLines)):
      elePhotTrigger = False
      if (hltResultsSummLines[j2].find("HLT_Ele")>-1) or (hltResultsSummLines[j2].find("HLT_Phot")>-1):
         elePhotTrigger=True
      elif (hltResultsSummLines[j2].find("HLT_DoubleEle")>-1) or (hltResultsSummLines[j2].find("HLT_DoublePhot")>-1):
         elePhotTrigger=True
      if (hltResultsSummLines[j2].find("Jet")>-1) or (hltResultsSummLines[j2].find("Tau")>-1):
         elePhotTrigger=False
      if headerInfo==True or elePhotTrigger==True:
         hltResultsElePhotSummLines.append(hltResultsSummLines[j2])
      if hltResultsSummLines[j2].startswith("TrigReport  Trig Bit#        Run     Passed     Failed      Error Name")==True: #Checks whether header info has ended
         headerInfo=False
   
   hltResultsElePhotSummFile = open(hltResultsElePhotSummFname,"w")
   hltResultsElePhotSummFile.writelines(hltResultsElePhotSummLines)
   hltResultsElePhotSummFile.close()
   print "  ...done."

