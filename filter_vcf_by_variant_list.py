#written by noah friedman
#given a vcf, call bcftools to filter a vcf to only include the variants included in the list
import subprocess
import pandas as pd
import sys

def filter_vcf_by_variant_list(variantTextFile, inputVcf, outputFileName):
	#first we need to bgzip the file
	#check if the file is not compressed, if it is compress it
	vcf = inputVcf
	if '.gz' not in vcf:
		cmd = 'bgzip {iVcf}'.format(iVcf = inputVcf)
		subprocess.Popen(cmd, shell=True).wait()
		vcf = inputVcf + '.gz' 
		cmd = 'tabix {iVcf}'.format(iVcf = vcf)
		subprocess.Popen(cmd, shell=True).wait()
	cmd = 'bcftools view {iVcf} -R {variantTFile} -o {outputFName}'.format(
			iVcf = vcf,
			variantTFile = variantTextFile,
			outputFName = outputFileName
		)
	print 'bcftools command: ', cmd
	subprocess.Popen(cmd, shell=True).wait()
	print outputFileName
	return outputFileName

#utility function to write a xls to a variant list
#hey alert, we probaly want to write this file to the GC work file area
def write_xls_to_variant_list(xlsName, udnId):
	xls = pd.ExcelFile(xlsName)
	df = xls.parse(xls.sheet_names[0]) 
	listToBeWritten = []
	for index, row in df.iterrows():
		print row
		listToBeWritten.append((row['Chromosome'], row['Position']))
	outputFName = udnId + '_variants.txt'
	listToBeWritten = sorted(listToBeWritten, key=lambda x: x[0]) #fix the order
	with open(outputFName, 'w') as f: #write to output
		for line in listToBeWritten:
			f.write(str(line[0]) + '\t' + str(line[1]) + '\n')
	return outputFName
