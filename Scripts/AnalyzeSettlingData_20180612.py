import sys
import os

Directory = "/Users/chris/Google Drive/Dunham Lab/Postdoc_Brewing/SettlingAssay/MinimumValues_20180621/"
FileNames = os.listdir(Directory)
WriteFile = open("/Users/chris/Google Drive/Dunham Lab/Postdoc_Brewing/SettlingAssay/SettlingRatios_20180621.txt", 'w+')

def MinimumValue(Array):
	MinValue = 999999999.9
	for i in Array:
		if MinValue > float(i[1]):
			MinValue = float(i[1])
	return MinValue
	
def MaximumValue(Array):
	MaxValue = 0.0
	for i in Array:
		if MaxValue < float(i[1]):
			MaxValue = float(i[1])
	return MaxValue	
def NameInfo(Name):
	SampleName = Name[0:Name[Name.find("_")+1:].find("_") + Name.find("_") + 1]
	Time = Name[Name[Name.find("_")+1:].find("_") + Name.find("_") + 2:Name.find(".")]
	Replicate = Name[Name.find("JPG") + 3:].replace(".txt","")
	return [SampleName, Time, Replicate]

for Name in FileNames:
	if Name == ".DS_Store":
		continue
	else:
		NameList = NameInfo(Name)
   		OpenFile = open(Directory + Name, 'r')
    	FileList = []
    	for Line in OpenFile:
    		FileList.append(Line.strip().split("\t"))
    	MaxValue = MaximumValue(FileList)
    	MinValue = MinimumValue(FileList)
    	LastLine = 0.0
    	PastMiddle = False
    	for Line in FileList:
    		if float(Line[1]) == MinValue:
    			PastMiddle = True
    		if (MaxValue/2) < float(Line[1]) and (MaxValue/2) > LastLine and PastMiddle: # Checks that we are past the minimum value, then checks if we are halfway through
    			WriteFile.write(NameList[0] + "\t" + NameList[1] + "\t" + NameList[2] + "\t" + str(float(Line[0])/len(FileList)) + "\n")
    			PastMiddle = False
    		LastLine = float(Line[1])
    		
    	OpenFile.close()