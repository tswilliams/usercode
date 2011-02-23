
# Setting which hltResultsSummary files to compare...
hltResultsSummFnames = []

hltResultsSummFnames.append("hltResults/11-02-12/hltResultsElePhotSummary_0-75TeV.log")
hltResultsSummFnames.append("hltResults/11-02-12/hltResultsElePhotSummary_1-00TeV.log")
hltResultsSummFnames.append("hltResults/11-02-12/hltResultsElePhotSummary_2-00TeV.log")

# Extracting the trigger path info (i.e. bit#, name, #events run over, #events passing) from each file...
trigPathBitNumCompList = []
numEvtsRunOverCompList = []
trigPathEffiCompList   = []
trigPathNameCompList   = []

for i in range(len(hltResultsSummFnames)):
	#Reading the contents of the i'th hltResultsSummary---.log file...
	hltResultsSummFile = open(hltResultsSummFnames[i],"r")
	fileContents = hltResultsSummFile.readlines()
	hltResultsSummFile.close()
	#Reducing this contents into a list of just the information for each of the trigger paths...
	convertLine=False
	trigPathBitNumList = [] #List of strings
	numEvtsRunOverList = [] #List of integers
	trigPathEffiList   = [] #List of floats
	trigPathNameList   = [] #List of strings
	for j in range(len(fileContents)):
		#Converting each table line into a list...
		if convertLine==True:
			lineCellContents = fileContents[j].split(None)
			#Extracting the four pieces of info that I want from this list, and appending them to the list of info from this hltResultsSummary file....
			trigPathBitNumList.append(lineCellContents[2])
			numEvtsRunOverList.append(lineCellContents[3])
			trigPathEffiList.append(str(float(lineCellContents[4])/float(lineCellContents[3])))
			trigPathNameList.append(lineCellContents[7])
			
		if fileContents[j-1].startswith("TrigReport ---------- Path   Summary ------------"):
			convertLine=True
	
	trigPathBitNumCompList.append(trigPathBitNumList)
	numEvtsRunOverCompList.append(numEvtsRunOverList)
	trigPathEffiCompList.append(trigPathEffiList)
	trigPathNameCompList.append(trigPathNameList)

for i2 in range(len(hltResultsSummFnames)):
	if i2!=0:
		for j2 in range(len(trigPathBitNumCompList[i2]):
			if trigPathBitNumCompList[i2][j2]!=trigPathBitNumCompList[i2-1][j2]:
				print "***Error!***"
				break
			if trigPathNameCompList[i2][j2]!=trigPathNameCompList[i2-1][j2]:
				print "***Error!***"
				break


ouputTableLineContents = []
firstLine = "Trigger\tBit#\t"
for i in range(len(hltResultsSummFnames)):
	fileTag = hltResultsSummFnames[i].
ouputTableLineContents.append(firstLine)
for j3 in range(len(trigPathBitNumCompList[0]))
	for i3 in range(len(hltResultsSummFnames))
		ouputTableLineContents
