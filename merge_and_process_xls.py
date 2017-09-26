#written by Noah Friedman
"""
Tools for processing the final xls for use by GCs and visualization
"""

import pandas as pd
import numpy as np
import sys
import os

#a dict mapping the column name we want for the final xls to [columnNamePipelineXls, columnNameUserXls]
columnMappings = {'CHROM': ['CHROM', 'Chromosome'], 'POS': ['POS', 'Position']}

#I do this with inefficient looping over both spreadsheets
#and looping over DFs
#this will not work if we have big spreadsheets, in that case we need to implement a pandas merge
def merge_and_add_columns(pipelineDf, userDf, userColumnsToAdd):
	colNames = pipelineDf.columns.tolist() + userColumnsToAdd
	numRows = len(pipelineDf.index)
	returnMergedDf = pd.DataFrame(np.empty((numRows, len(colNames)), dtype=str), columns=colNames) 
	for index, row in pipelineDf.iterrows():
		key = str(row[columnMappings['CHROM'][0]]) + ':' + str(row[columnMappings['POS'][0]])
		for idx, r in userDf.iterrows():
			userKey = str(r[columnMappings['CHROM'][1]]) + ':' + str(r[columnMappings['POS'][1]])
			if key == userKey:
				#copy over values from pipeline df
				for column in pipelineDf.columns.tolist():
					returnMergedDf.set_value(index, column, row[column])
				#copy over values from userDf
				for column in userColumnsToAdd:
					returnMergedDf.set_value(index, column, r[column])
	return returnMergedDf

#utility function to convert a value to float if an only if it can be
def convert_to_float_or_zero(i):
	try: 
		float(i)
		return float(i)
	except:
		return 0

def add_allele_freq_summary_column(df):
	alleleFreqCols = ['AF_EAS', 'AF_NFE', 'AF_SAS', 'AF_AMR', 'AF_AFR']
	df['GNOMAD_Max_Allele_Freq'] = np.empty(len(df.index))
	for index, row in df.iterrows():
		freqs = [convert_to_float_or_zero(row[i]) for i in alleleFreqCols]
		population = alleleFreqCols[freqs.index(max(freqs))]
		freq = max(freqs)
		#ALERT todo: include the population from which the max freq comes
		#df.set_value(index, 'Max_Allele_Freq', str(population) + ':' + str(freq))
		df.set_value(index, 'GNOMAD_Max_Allele_Freq', freq)
	#print df

def add_tier_column():
	return 0

#read xls sheets and create a dictionary mapping sheet names to dictionaries
#it also returns a list of sheetNames
def read_xls_sheets(xlsName):
	xls = pd.ExcelFile(xlsName)
	sheetNames = xls.sheet_names
	sheetDict = dict()
	for sheet in sheetNames:
		sheetDict[sheet] = xls.parse(sheet)
	return sheetDict, sheetNames

def merge_columns_across_spreadsheets(spreadSheetPipeline, spreadSheetUser, outputDir):
	sheetDictPipeline, sheetDictPipelineNames = read_xls_sheets(spreadSheetPipeline)
	sheetDictUser, sheetDictUserNames = read_xls_sheets(spreadSheetUser)
	#We expect the user's xls to be just a single sheet output from ingenuity.  If its not, that breaks our code and we exit
	if len(sheetDictUser) != 1:
		print 'error we except an excel sheet from the user with a single sheet'
		sys.exit()
	mergedDf = merge_and_add_columns(sheetDictPipeline[sheetDictPipelineNames[1]], sheetDictUser[sheetDictUserNames[0]], [ #indicies indicate where we can find the two actual data spreadsheets
	'Transcript ID', 'Transcript Variant', 'Protein Variant', 'Gene Region', 'Gene Symbol']) #list of columns to add from the user uploaded columns
	add_allele_freq_summary_column(mergedDf)

	#save everything to an xlsx
	outputXlsxName = os.path.join(outputDir, 'outputTestMerged.xlsx')	 
	writer = pd.ExcelWriter(outputXlsxName)
	sheetDictPipeline[sheetDictPipelineNames[0]].to_excel(writer, 'Column Descriptions', index = False)
	mergedDf.to_excel(writer,'Sheet1', index = False)
	writer.save()
	return outputXlsxName



