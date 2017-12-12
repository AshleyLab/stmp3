#written by Noah Friedman
"""
Tools for processing the final xls for use by GCs and visualization
"""

import pandas as pd
import numpy as np
import sys
import os


reload(sys)
sys.setdefaultencoding('utf-8')
sys.setdefaultencoding("latin-1")  #ingenuity exports are encoded this way I think

#a dict mapping the column name we want for the final xls to [columnNamePipelineXls, columnNameUserXls]
columnMappings = {'CHROM': ['CHROM', 'Chromosome'], 'POS': ['POS', 'Position']}

#I do this with inefficient looping over both spreadsheets
#and looping over DFs
#this will not work if we have big spreadsheets, in that case we need to implement a pandas merge
def merge_and_add_columns(df1, df2, userColumnsToAdd, idxForColMappingsdf1, idxForColMappingsdf2):
	#now do the merge
	colNames = df1.columns.tolist() + userColumnsToAdd
	numRows = len(df1.index)
	returnMergedDf = pd.DataFrame(np.empty((numRows, len(colNames)), dtype=str), columns=colNames) 
	for index, row in df1.iterrows():
		df1key = str(row[columnMappings['CHROM'][idxForColMappingsdf1]]) + ':' + str(row[columnMappings['POS'][idxForColMappingsdf1]])
		for idx, r in df2.iterrows():
			df2Key = str(r[columnMappings['CHROM'][idxForColMappingsdf2]]) + ':' + str(r[columnMappings['POS'][idxForColMappingsdf2]])
			#print df2Key
			#print '__________'
			if df1key == df2Key:
				#copy over values from pipeline df
				for column in df1.columns.tolist():
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
	df['GNOMAD_Max_Allele_Freq_POP'] = np.empty(len(df.index))
	df['GNOMAD_Max_Allele_Freq_POP'] = df['GNOMAD_Max_Allele_Freq_POP'].astype(str) #this will be a column of strings
	for index, row in df.iterrows():
		freqs = [convert_to_float_or_zero(row[i]) for i in alleleFreqCols]
		population = alleleFreqCols[freqs.index(max(freqs))][3:]
		freq = max(freqs)
		#ALERT todo: include the population from which the max freq comes
		#df.set_value(index, 'Max_Allele_Freq', str(population) + ':' + str(freq))
		df.set_value(index, 'GNOMAD_Max_Allele_Freq', freq)
		df.set_value(index, 'GNOMAD_Max_Allele_Freq_POP', population)
	#print df

def fix_ref_or_alt_column_for_indels(df, refKey, altKey):
	for index, row in df.iterrows():
		if type(row[altKey]) == float: 
			df.set_value(index, altKey, '_')
		if type(row[refKey]) == float: 
			df.set_value(index, refKey, '_')

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

def sort_sheets(df):
	print df
	sys.exit()
	print df.sort([''])

#arbitrary function for fixing issues with STMP/ingenuity inconsistency
#basically if there is a stmp pos that is close to or equal to ingenuity pos overwrite the ingenuity pos with the stmp pos
def find_closest_match_for_pos(stmpDf, ingenuityChrom, ingenuityPos):
	minDistance = 100000 #default distance is arbitrary big number
	newChrom, newPos, newRef, newAlt = -1, -1, -1, -1
	for idx, r in stmpDf.iterrows():
		stmpChrom, stmpPos = str(r['CHROM']), str(r['POS'])
		stmpRef, stmpAlt = str(r['REF']), str(r['ALT'])
		if stmpChrom == ingenuityChrom:
			distance = abs(int(stmpPos) - int(ingenuityPos))
			if distance < minDistance:
				minDistance = distance
				newChrom = stmpChrom
				newPos = stmpPos
				newRef = stmpRef
				newAlt = stmpAlt
	return newChrom, newPos, newRef, newAlt

#in order to properly merge with the ingenuity XLS, we need to right align our STMP output
#not generic function iterates over both dfs
def realign_ingenuity_sheet(ingenuityDf, stmpDf):
	for index, row in ingenuityDf.iterrows():
		ingenuityChrom, ingenuityPos, ingenuityRef, ingenuityAlt = str(row['Chromosome']), str(row['Position']), str(row['Reference Allele']), str(row['Sample Allele'])
		newChrom, newPos, newRef, newAlt = find_closest_match_for_pos(stmpDf, ingenuityChrom, ingenuityPos)
		if newChrom > 0: #if we found a match please set it up
			ingenuityDf.set_value(index, 'Chromosome', newChrom)
			ingenuityDf.set_value(index, 'Position', newPos)
			ingenuityDf.set_value(index, 'Reference Allele', newRef)
			ingenuityDf.set_value(index, 'Sample Allele', newAlt)
	return ingenuityDf

def merge_columns_across_spreadsheets(spreadSheetPipeline, spreadSheetUser, outputDir, udnId):
	sheetDictPipeline, sheetDictPipelineNames = read_xls_sheets(spreadSheetPipeline)
	sheetDictUser, sheetDictUserNames = read_xls_sheets(spreadSheetUser)

	#fix_ref_or_alt_column_for_indels(sheetDictUser[sheetDictUserNames[0]], 'Reference Allele', 'Sample Allele') #this code fixes the syntax of indels for the user inputed sheet
	#fix_ref_or_alt_column_for_indels(sheetDictPipeline[sheetDictPipelineNames[1]], 'REF', 'ALT') #this code fixes the syntax of indels for the pipeline inputed sheet
	sheetDictUser[sheetDictUserNames[0]] = realign_ingenuity_sheet(sheetDictUser[sheetDictUserNames[0]], sheetDictPipeline[sheetDictPipelineNames[1]])

	#We expect the user's xls to be just a single sheet output from ingenuity.  If its not, that breaks our code and we exit
	if len(sheetDictUser) != 1:
		print 'error we expect an excel sheet from the user with a single sheet'
		sys.exit()
	mergedDfPipeline = merge_and_add_columns(sheetDictPipeline[sheetDictPipelineNames[1]], sheetDictUser[sheetDictUserNames[0]], [ #indicies indicate where we can find the two actual data spreadsheets
	'Transcript ID', 'Transcript Variant', 'Protein Variant', 'Gene Region', 'Gene Symbol'], 0, 1) #list of columns to add from the user uploaded columns
	
	colsToAddToDf = ['AF_EAS', 'AF_NFE', 'AF_SAS', 'AF_AMR', 'AF_AFR', 'NC', 'NI', 'NA', 'ESP_AF_POPMAX', 'KG_AF_POPMAX', 'SD', 'SF', 'QUAL', 'ID', 'FILTER', 'GT', 'NJ', 'SX', 'GI', 'AN_AFR', 'AN_AMR', 'AN_ASJ', 'AN_EAS', 'AN_FIN', 'AN_NFE', 'AN_OTH', 'AN_SAS', 'clinvar_pathogenic', 'KG_AF_GLOBAL', 'KG_AC', 'POPMAX', 'AN_POPMAX', 'AC_POPMAX', 'AF', 'AN', 'AN_Female', 'AN_Male']
	colsFinal = []
	for x in colsToAddToDf: #fix the columns in case they are missing
		if x in sheetDictPipeline[sheetDictPipelineNames[1]].columns.tolist():
			colsFinal.append(x)

	mergedDfUser = merge_and_add_columns(sheetDictUser[sheetDictUserNames[0]], sheetDictPipeline[sheetDictPipelineNames[1]], colsFinal, 1, 0) # the 0s /1s relate to is it the GC xls first or the pipeline xls first
	add_allele_freq_summary_column(mergedDfUser)

	if False:
		sort_sheets(mergedDfUser)

	#save everything to an xlsx
	outputXlsxName = os.path.join(outputDir, udnId + '_merged.xlsx')

	writer = pd.ExcelWriter(outputXlsxName,options={'encoding':'latin-1'})

	#sheetDictPipeline[sheetDictPipelineNames[0]].to_excel(writer, 'Column Descriptions', index = False)
	mergedDfUser.to_excel(writer,'Sheet1', index = False)

	writer.save()
	return outputXlsxName


renameDict = {'SX': 'SwissProtExpression', 'GI': 'ProximalGeneInfo', 'SD': 'SwissProtDiseaseAssociation', 'mTaster': 'MutationTaster', 'SF': 'SwissProtFunction', 'phylop': 'phyloP', 'NI': 'MutationTasterPVal', 'sift': 'Sift', 'GNOMAD_Max_Allele_Freq': 'GNOMADMaxAlleleFreq', 'POPMAX': 'ExacPopmax', 'AN_POPMAX': 'ExacANPopmax', 'AC_POPMAX': 'ExacACPopmax', 'AF':'ExacAf', 'AN':'ExacAn', 'AN_Female': 'GNOMAD_AN_FEMALE', 'AN_Male': 'GNOMAD_AN_MALE'}
#utility function we use to make the columns of the xls human readable as needed
def make_cols_human_readable(df): 
	colsToRename = renameDict
	#make sure we dont try to rename a column that isnt actually there
	for key, value in colsToRename.items():
		if key not in df.columns.tolist():
			del colsToRename[key]
	df = df.rename(columns=colsToRename)
	return df

#check each cell.  If the cell has a comma and the value before the comma is a float, just take the first half
def clean_up_comma_separated_values(df):
	cols = df.columns
	for index, row in df.iterrows():
		for col in cols:
			commaIdx = str(row[col]).find(',')
			if commaIdx > 0:
				val = row[col][:commaIdx]
				try:
					float(val)
					df.set_value(index, col, val)
				except:
					pass
	return df

def clean_cell_values(df):
	df = clean_up_comma_separated_values(df)
	valsToMarkAsEmpty = ['no value', '.']
	emptyValue = ''
	for col in df.columns:
		if col == 'CADD Score':
			df[col] = df[col].replace(['< 10'], 1)  #clean CADD score
		df[col] = df[col].replace(valsToMarkAsEmpty, emptyValue) #normalize empty columns 
		#clean up commas separated values
	return df

orderedCols = ['Chromosome', 'Position', 'Reference Allele', 'Sample Allele', 'Variation Type', 'Gene Region', 
'Gene Symbol', 'Transcript ID', 'Transcript Variant', 'Protein Variant', 'Translation Impact',
'GNOMAD_Max_Allele_Freq'
]
def sort_cols(df):
	#inefficient way to order the columns
	firstCols = [] #all columns whose order we care about
	lastCols = [] #put all columns we dont specify specifically at the end
	for col in orderedCols:
		if col in df.columns.tolist():
			firstCols.append(col)
	for col in df.columns.tolist():
		if col not in firstCols:
			lastCols.append(col)
	sortedCols = firstCols + lastCols
	df = df.reindex_axis(sortedCols, axis=1)
	return df


#function to rename columns, reorder columns etc in xls data
def improve_legibility_of_xls(xlsName):
	xls = pd.ExcelFile(xlsName)
	df = xls.parse(xls.sheet_names[0]) 
	#to imporve legibility we do two things: make columns readable and sort columns
	df = make_cols_human_readable(df)
	df = sort_cols(df)
	df = clean_cell_values(df)
	outputXlsxName = xlsName
	writer = pd.ExcelWriter(outputXlsxName,options={'encoding':'latin-1'})
	df.to_excel(writer,'Sheet1', index = False)
	writer.save()
	return outputXlsxName



