'''

Convert vcf to allele frequency
Chris Large
Take a vcf from samtools and converts it into the ratio of reference and alternate allele frequencies.

python calc_het_wholegenome.py input.vcf output.vcf

'''

import sys

vcfFile = open(sys.argv[1], 'r')
outFile = open(sys.argv[2], 'w')

chrDict = {'I':['1',0,230218],
	'II':['2',230218,813184],
	'III':['3',1043402,316620],
	'IV':['4',1360022,1531933],
	'V':['5',2891955,576874],
	'VI':['6',3468829,270161],
	'VII':['7',3738990,1090940],
	'VIII':['8',4829930,562643],
	'IX':['9',5392573,439888],
	'X':['10',5832461,745751],
	'XI':['11',6578212,666816],
	'XII':['12',7245028,1078177],
	'XIII':['13',8323205,924431],
	'XIV':['14',9247636,784333],
	'XV':['15',10031969,1091291],
	'XVI':['16',11123260,948066],
	'M':''}

chromosome = 1

for VCFLine in vcfFile:
	if(VCFLine[0] =='#'): # Skips header
		continue
	else: 
		tabVCFLine = VCFLine.split('\t')
		if tabVCFLine[0] == 'chrM':
			continue
		print tabVCFLine[0] + str(int(chrDict[tabVCFLine[0][tabVCFLine[0].index('chr',0)+3:]][2]) - int(tabVCFLine[1]))
		if int(tabVCFLine[1]) < 30000 or (int(chrDict[tabVCFLine[0][tabVCFLine[0].index('chr',0)+3:]][2]) - int(tabVCFLine[1]) < 30000):
			print "skip"
			continue
		
		VCFposition = int(tabVCFLine[1]) + int(chrDict[tabVCFLine[0][tabVCFLine[0].index('chr',0)+3:]][1])
		
		outFile.write(str(chrDict[tabVCFLine[0][tabVCFLine[0].index('chr',0)+3:]][0]) + '\t' + str(VCFposition) + '\t') # Converts chrI to 1
		
		infoDelim = tabVCFLine[7][4 + tabVCFLine[7].find('DP4='):tabVCFLine[7].find('MQ=') - 1] # Finds DP4 in INFO tab
		DP4Delim = infoDelim.split(', ')
		numerator = float(DP4Delim[0]) + float(DP4Delim[1])
		denominator = numerator + float(DP4Delim[2]) + float(DP4Delim[3])
		
		if denominator == 0:
			outFile.write('1' + '\n')
		else:
			outFile.write(str(numerator/denominator) + '\n')



