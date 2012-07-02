#!/usr/bin/python

###########################################
# Arguments:
#      * Directories of CRAB tasks
# When run:
#      * Produces file containing a non-duplicate list of output files from this CRAB task
#          ('non-duplicate' => if any one of the jobs in a CRAB task is run multiple times, 
#                              then only the output file from the latest 'running' of that 
#                              job is written to the list. Thus, the analysis code does not
#                              run over the same events multiple times.)
#      * Checks how many output files there should be and reports if any are missing
#      * ... 
# Author: T. S. Williams
# Date: June 2012
############################################

import subprocess, time
import re
import datetime

# Define functions 

def dateStamp(dateFormat='%Y%m%d'):
    return datetime.datetime.now().strftime(dateFormat)

def timeStamp(timeFormat='%H%M%S'):
    return datetime.datetime.now().strftime(timeFormat)

def dateTimeStamp():
   return dateStamp()+'-'+timeStamp()

def runProcess(shellCmd):
   maxExecTime = 300.0
   subproc = subprocess.Popen(shellCmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
   t0 = time.time()
   while( (time.time()-t0)<maxExecTime):
      time.sleep(5.0)
      retcode = subproc.poll() #returns None while subprocess is running
      if (retcode is not None):
         return subproc.stdout.readlines()
   print " ERROR: Running command '"+exe+"' in function runProcess took longer than",maxExecTime,"sec"
   print "        An empty list will now be returned !"
   return []

def GetNumJobsInCrabTask(crabTaskDir):
   args_xml_file = open(crabTaskDir+"/share/arguments.xml")
   args_xml_contents = args_xml_file.read()
   args_xml_file.close()
   return args_xml_contents.count('</Job>')
   
def GetCrabTaskOutputDir(crabTaskDir):
   crab_cfg_file = open(crabTaskDir+"/share/crab.cfg")
   crab_cfg_contents = crab_cfg_file.readlines()
   crab_cfg_file.close()
   for ithLine in crab_cfg_contents:
      if ithLine.count("user_remote_dir")>0:
         return ithLine.split("=",1)[1].rstrip('\n')

def GetCrabTaskOutputFileName(crabTaskDir):
   crab_cfg_file = open(crabTaskDir+"/share/crab.cfg")
   crab_cfg_contents = crab_cfg_file.readlines()
   crab_cfg_file.close()
   for ithLine in crab_cfg_contents:
      if ithLine.count("output_file")>0:
         return ithLine.split("=",1)[1]

def WriteCrabOutputFileList(crabTaskDir):
   print "Crab task:", crabTaskDir
   
   # Check from CRAB directory how many output files there should be  
   numJobs = GetNumJobsInCrabTask(crabTaskDir)
   print "   * Total num of jobs =", numJobs
   
   # Retrieve output directory for this task
   outputDir = GetCrabTaskOutputDir(crabTaskDir)
   outputDir = "/pnfs/pp.rl.ac.uk/data/cms/store/user/tsw/"+outputDir
   print "   * Output dir =", outputDir
   
   # Get list of file in output directory
   outputFileNames = runProcess("ls "+outputDir)
   outputFileNames = [eachLine.rstrip("\n") for eachLine in outputFileNames] # Strip newline character from end of each line using a list comprehension
   
   # Setup corresponding list of job no, 'job resubmit no' tuples 
   thePattern = re.compile('\S+_(\d+)_(\d+)_\w\w\w.root')
   outputFilesInfo = [ [int(thePattern.match(eachLine).group(1)), int(thePattern.match(eachLine).group(2))] for eachLine in outputFileNames]
   
   # Tie these lists together and sort by (jobNo, jobResubNo) tuples
   combinedList = zip(outputFilesInfo, outputFileNames)
   combinedList.sort()
   combinedList.reverse()
   
   # Now go through this list and extract the non-duplicate output file list
   nonDuplOutputFiles  = [] 
   jobsWithOutputFiles = []
   prevJobNo = -1
   for ithItem in combinedList:
      fileName     = ithItem[1]
      jobNoForFile = ithItem[0][0]
#      print ithItem
      if (jobNoForFile != prevJobNo):
#         print " * Will add this file to my list"
         nonDuplOutputFiles.append(fileName)
         jobsWithOutputFiles.append(jobNoForFile)
      prevJobNo = jobNoForFile
   nonDuplOutputFiles.reverse()
#   print "The list of non-duplicate files is:"
#   for fileName in nonDuplOutputFiles:
#      print "    +", fileName

   numOutFilesMissing = numJobs - len(nonDuplOutputFiles)
   jobsWithoutOutputFiles = set(range(1,numJobs+1)) - set(jobsWithOutputFiles)
   jobsWithoutOutputFiles = list(jobsWithoutOutputFiles)
   jobsWithoutOutputFiles.sort()
   if numOutFilesMissing>0:
      print "   * There are",numOutFilesMissing,"jobs without any output file yet"
      print "        -> Job numbers:", jobsWithoutOutputFiles
   elif numOutFilesMissing==0:
      print "   * Output file found for each job"
   else:
      print "  ERROR: In function WriteCrabOutputFileList (crabTaskDir='"+crabTaskDir+"')"
      print "         numOutFilesMissing(="+str(numOutFilesMissing)+") is NEGATIVE"
      print "         This does not make sense - function will return BEFORE writing non-duplicate filelist to disk" 
      return
   if len(jobsWithoutOutputFiles) != numOutFilesMissing:
      print "  ERROR: In function WriteCrabOutputFileList (crabTaskDir='"+crabTaskDir+"')"
      print "         Length of list jobsWithoutOutputFiles (="+str(len(jobsWithoutOutputFiles))+") does *NOT* equal numOutFilesMissing(="+str(numOutFilesMissing)+")"
      print "         This does not make sense - function will return BEFORE writing non-duplicate filelist to disk"
      return
   
   # Now, write nonDuplOutputFiles to a file
   strDateTime = dateTimeStamp()
   headerLines = ['## nonDuplFileList\n']
   headerLines.append('## Created by python function WriteCrabOutputFileList with arg crabTaskDir="'+crabTaskDir+'"\n')
   headerLines.append('## Date: '+strDateTime+'\n')
   headerLines.append('## numFilesMissing = '+str(numOutFilesMissing)+'\n')
   
   fileListFileName = crabTaskDir+'/nonDuplFiles_'
   if numOutFilesMissing==0:
      fileListFileName += 'all'
   else:
      fileListFileName += str(numOutFilesMissing)+'JobsMissing'
   fileListFileName +='_'+strDateTime+'.txt' 
   theFile = open(fileListFileName, "w")
   theFile.writelines(headerLines)
   
   fileListLines = [ (outputDir+'/'+eachLine+'\n') for eachLine in nonDuplOutputFiles ]
   theFile.writelines(fileListLines)
   
   theFile.close()
   
   print "   * Non-duplicate list of output files written to:",fileListFileName
   
   

   
   
#################################################
## Main script part

# Imports
import sys, os

# Write file list for each output file
for arg in sys.argv[1:]:
   print
   if( os.path.isdir(arg) ):
      WriteCrabOutputFileList(arg)
   else:
      print "  WARNING: Script detected that argument '"+arg+"' is NOT a directory"
      print "           It will now be skipped"
   
