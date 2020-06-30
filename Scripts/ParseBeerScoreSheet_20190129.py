import os
import sys

beerScoreSheetFile = open(sys.argv[1],'r')

beerScoreSheetLists = []
''' Makes a set of lines that if it begins with, it skips'''
linePrefixesEliminate = ["Beer", "Aroma", "Appearance", "Flavor", "Flaws", "Mouthfeel"]
linePrefixesTwoRow = ["Malt","Hops","Esters","Phenols","Alcohol","Sweetness","Acidity","Other", "Aspect", "Clarity", "Head Size", "Head Retention", "Head Texture", "Sweetness", "Bitterness", "Harshness", "Fault"]
lineThreePrefixesTwoRow = ["Brown", "Black"]
linePrefixesThreeRow = ["Acetaldehyde","Hot/Alcoholic","Astringent","Diacetyl","DMS","Estery","Grassy","Light-struck","Medicinal","Metallic","Musty","Oxidized","Plastic","Solvent/Fusel","Sour/Acidic","Smoky","Spicy","Sulfur","Vegetal","Vinegary","Yeasty"]

for line in beerScoreSheetFile:
	beerScoreSheetLists.append(line.strip().split('\t'))

judgeNumber = 1
sys.stdout.write("Judge\tBeer\tCharecteristic\tScore_1\tScore_2\tScore_3\n")
for line in beerScoreSheetLists:
	#print line
	if line[0] == 'Judge:':
	## Set judge number on judge line
		judgeNumber = line[1].strip()
		rowNumber = 0
		continue
	elif line[0].startswith(tuple(linePrefixesEliminate)):
	## Skip the aforementioned lines
		continue
	elif line[0].startswith("Aspect") and line[2].startswith("Malt traits?"):
	## Skip the aspect line if on the same line as Malt traits?
		continue
	elif len(line) < 6:
	## Skips the empty liness
		continue
	sys.stdout.write(judgeNumber + '\t' + '1' + '\t')
	columnCounter = 0
	for columnNumber in range(len(line)): 
	## Able to carriagereturn the 2nd and 3rd beer lines
		if 24 > columnNumber >= 12:
			beerNumber = "2"
		elif columnNumber >= 24:
			beerNumber = "3"
		else:
			beerNumber = "1"
		if columnNumber == 12 or columnNumber == 24:
			sys.stdout.write('\n' + judgeNumber + '\t' + beerNumber+ '\t')
			columnCounter = 0
		## Outputs either every other line with something depending on what the line starts with
		if line[0].startswith(tuple(linePrefixesThreeRow)) and columnCounter > 3 and columnNumber % 2 > 0:
			sys.stdout.write(line[columnNumber] + '\n' + judgeNumber + '\t' + beerNumber + '\t')
		elif line[2].startswith(tuple(lineThreePrefixesTwoRow)) and columnNumber % 2 > 0:
			sys.stdout.write(line[columnNumber] + '\n' + judgeNumber + '\t' + beerNumber + '\t')
		elif line[0].startswith(tuple(linePrefixesTwoRow)) and columnNumber % 2 > 0:
			sys.stdout.write(line[columnNumber] + '\n' + judgeNumber + '\t' + beerNumber + '\t')
		else:
			sys.stdout.write(line[columnNumber] + '\t')
		columnCounter += 1
	if len(line) > 0:
		sys.stdout.write('\n')
	rowNumber += 1

	