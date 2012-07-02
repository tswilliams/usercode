#!/usr/bin/python

#############################################################
#
#   2 arguments - config file & output dir
#
# Author: Tom Williams
# Date: June 2012
##############################################################

# Imports
import sys, re

if len(sys.argv)!=3:
   print " ERROR : Incorrect number of arguments given. Script exits early ..."
   print "         There should be two arguments:"
   print "           + arg 1 : Location of config file"
   print "           + arg 2 : dir that ouput files will be copied to"
   sys.exit()

configFile_name = sys.argv[1]
baseOutputDir = sys.argv[2]

configFile = open(configFile_name)
configFile_lines = configFile.readlines()
configFile.close()

for ithLine in configFile_lines:
   # Skip lines beginning with comment symbol '#' 
   if ithLine[0]=="#":
      continue
   elif ithLine=="\n":
      continue
   ithLine = ithLine.rstrip("\n")
   ithLine = ithLine.rstrip(" ")
   
   # Split line into different config params
   listOfCfgParams = re.split("\s*;\s*", ithLine)
   numCfgParams = len(listOfCfgParams)
   
   if numCfgParams<3 or numCfgParams>5 : 
      print " ERROR : Incorrect number of config parameters on current line. Script exits early ..."
      print "         Line was '"+ithLine+"'"
      sys.exit()
      
   # Check that first param is mc or data
   mcOrData = listOfCfgParams[0]
   mcOrData = mcOrData.lower()
   if mcOrData=="mc":
      isMC = True
   elif mcOrData=="data":
      isMC = False
   else:
      print " ERROR : Invalid 1st config param on current line. Script exits early ..."
      print "         Line was '"+ithLine+"'"
      sys.exit()
      
   # Parse the other params
   anrArgs_inputFileName=listOfCfgParams[1]
   
   anrArgs_outputFileTag=listOfCfgParams[2]
   if anrArgs_outputFileTag.endswith(".root"):
      anrArgs_outputFileTag = anrArgs_outputFileTag.replace(".root", "")
   if anrArgs_outputFileTag.find("/")!=-1:
      print " ERROR : Invalid 3rd config param on current line."
      print "         outputFileTag should NOT contain a `/`"
      print "         Line= `"+ithLine+"`"
      print "         outputFileTag ... parsed value: '"+anrArgs_outputFileTag+"'"
      print "         Script exiting early !"
      sys.exit()
      
   if numCfgParams>3:
      anrArgs_numJobs = int(listOfCfgParams[3])
      if anrArgs_numJobs<=0:
         print " ERROR : Invalid 4th config param on current line (Total number of jobs '"+str(anrArgs_numJobs)+"' cannot be <=0)."
         print "         Line was '"+ithLine+"'"
         print "         Script exiting early !"
         sys.exit()
   else:
      anrArgs_numJobs = 1

      
   if numCfgParams>4:
      anrArgs_maxEvents = int(listOfCfgParams[4])
   else:
      anrArgs_maxEvents = -1
   
   # Construct the arguments for BstdZeeFirstAnalyser
   if isMC:
      anrArgsString = " --mc"
   else:
      anrArgsString = " "
   
   anrArgsString += " --from-list -i "+anrArgs_inputFileName
   anrArgsString += " --maxEvents "+str(anrArgs_maxEvents)
   
   
   # Print out jobName ; jobLogFile ; Anr options
   jobName = "BstdZAna_"+anrArgs_outputFileTag
   jobLogFile = baseOutputDir+"/"+anrArgs_outputFileTag+".log"
   
   for subSampleNo in range(1,anrArgs_numJobs+1):
      ithJobName            = jobName
      ithJobLogFile         = jobLogFile
      ithAnrArgs_outFileTag = anrArgs_outputFileTag
      ithAnrArgsString      = anrArgsString
      
      if anrArgs_numJobs != 1:
         ithJobName += ("_"+str(subSampleNo))
         ithJobLogFile = ithJobLogFile.replace(".log","_"+str(subSampleNo)+".log")
         ithAnrArgs_outFileTag += ("_"+str(subSampleNo))
         ithAnrArgsString += " -o "+ithAnrArgs_outFileTag+" --tot-num-jobs="+str(anrArgs_numJobs)+" --num-this-job="+str(subSampleNo)
      else:
         ithAnrArgsString += " -o "+ithAnrArgs_outFileTag
         
      lineToPrint = ithJobName+" ; "+ithJobLogFile+" ; "+ithAnrArgsString
      print lineToPrint
   
   
      
   
      
   