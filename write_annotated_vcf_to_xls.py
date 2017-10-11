#written by noah friedman
import sys
import pandas as pd
import re
import numpy as np
import os

reload(sys)
sys.setdefaultencoding("latin-1")
#sys.setdefaultencoding('utf-8')

#VCF indicies: CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT_TAGS_NAMES	FORMAT_TAGS_VALUES	UDN78139
#indexDict = {'CHROM':0,'POS':1,'ID':2,'REF':3,'ALT':4,'QUAL':5,'FILTER':6,'INFO':7,'FORMAT':8,'UDNHeader':9}

def parse_info_header(header, infoFieldsDict):
	header = header[11:-2] #string off the leading ##INFO=< and trailing >
	headerVals = header.split(',') #split by comma
	#we do a bit of extra processing for the information in the description tag because this will become the column names for our XLS
	headerVals[3] = headerVals[3].strip('Description=').strip('\"')
	infoFieldsDict[headerVals[0]] = headerVals[1:] 
	#infoFieldsList.append(headerVals)
	#return infoFieldsDict

#initializes the df that is the first sheet of our xls which explains what column headers mean 
def initialize_info_df(infoFieldsDict, formatTags):
	otherColNames = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER']
	numRows = len(infoFieldsDict) + len(formatTags) + len(otherColNames) #figure out how many rows will be in the df
	df = pd.DataFrame(np.empty((numRows, 2), dtype=str), columns=['Column header', 'Description']) #initialize df and all cells strings
	cntr = 0
	for c in otherColNames:
		df.set_value(cntr, 'Column header', c)
		df.set_value(cntr, 'Description', c)
		cntr += 1
	for f in formatTags: 
		df.set_value(cntr, 'Column header', f)
		df.set_value(cntr, 'Description', f)
		cntr += 1
	for key, value in infoFieldsDict.items():
		df.set_value(cntr, 'Column header', key)
		df.set_value(cntr, 'Description', value[2])
		cntr += 1
	return df

#initializes the pandas df object we will temporarily be storing all of our data in
def initialize_df(infoFieldsDict, formatTags):
	colNames = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER']
	for f in formatTags: colNames.append(f)
	for key, value in infoFieldsDict.items():
		#print key, infoFieldsDict[key]
		colNames.append(key) #append key as the column
	#add chrom pos etc
	df = pd.DataFrame(columns=colNames)
	df = df.loc[:,~df.columns.duplicated()] #remove duplicate column names--this happens to us.  Solution from: https://stackoverflow.com/questions/14984119/python-pandas-remove-duplicate-columns
	return df

def set_info_fields_for_row_of_df(df, addIdx, infoFieldsRow, infoFieldsDict):
	for v in infoFieldsRow:
		iFieldDictKey = v[:v.find('=')] #we get the key by substringing the entry from beginning to the equals sign
		value = v[v.find('=') + 1:]
		key = iFieldDictKey
		df.set_value(addIdx, key, value) #set the cell of the current row of the df

#sets chrom, pos, ref, alt etc for the df
def set_core_cols_for_row_of_df(df, addIdx, coreFieldsRow):
	cntr = 0
	for key in ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER']:
		value = coreFieldsRow[cntr]
		df.set_value(addIdx, key, value)
		cntr += 1

#sets the cells of the df corresponding to the format tags (not if the format tags are not the same for every record this will break!! alert)
def set_format_tags_for_row_of_df(df, addIdx, formatTagsNames, formatTagsRow):
	cntr = 0
	for key in formatTagsNames:
		value = formatTagsRow[cntr]
		df.set_value(addIdx, key, value)
		cntr += 1

def set_row_of_df(df, colNames, splitLineFields, infoFieldsDict):
	#first add an empty row to the df:
	addIdx = len(df.index)
	coreFieldsRow = splitLineFields[:7]
	infoFieldsRow = splitLineFields[7].split(';')
	formatTagsNames = splitLineFields[8].split(':')
	formatTagsRow = splitLineFields[9].split(':')
	df.loc[addIdx] = ['no value' for i in range(len(colNames))] #initialize the row to say 'no value' everywhere
	set_core_cols_for_row_of_df(df, addIdx, coreFieldsRow)
	set_info_fields_for_row_of_df(df, addIdx, infoFieldsRow, infoFieldsDict)
	set_format_tags_for_row_of_df(df, addIdx, formatTagsNames, formatTagsRow)

def reorder_dataframe_columns():
	return 0
	cols = df.columns.tolist()

def only_keep_specified_columns(df, specifiedColumns):
	colsToDrop = set(df.columns) - set(specifiedColumns) #take the set difference to get columns to drop
	df = df.drop(colsToDrop, 1)
	return df

#for use on the info df to only keep specific rows we want
def only_keep_specified_rows(df, specifiedRows):
	#assumes that the rows are indexed by integers
	idxsToRemove = []
	cntr = 0
	for index, row in df.iterrows():
		if row['Column header'] not in specifiedRows: idxsToRemove.append(cntr)
		cntr += 1
	df = df.drop(idxsToRemove, 0)
	return df

#global variable: which columns we want to include
columnsToInclude = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'DP', 'GT', 'NC', 'ENSGN', 'ARIC_AA', 'TAMR', 'NA', 'AC_AFR', 'AC_EAS', 'AC_AMR', 'ARIC_EA', 'NG', 'NI', 'ENS', 'MT', 'AF_AFR', 'AF_AMR', 'AF_ASJ', 'AF_EAS', 'AF_FIN', 'AF_NFE', 'AF_OTH', 'AF_SAS', 'KG_AF_POPMAX', 'ESP_AF_POPMAX', 'NE', 'SX', 'CLNSG']
#only keep specified columns in our output
onlyKeepSpecifiedColumns = False

#main loop through the vcf file
#goes through each line, first gets columns from the info fields then headers which it exports
def vcf_to_xls(vcfFilePath, outputDir, udnId):
	print 'writing ', vcfFilePath, ' to xls format'
	df = None #declare df variable in global scope
	infoDf = None #declare infoDf variable in global scope
	with open(vcfFilePath) as f:
		firstRecordFlag = True #flag to trigger the initialization of the df
		lines = f.readlines()
		infoFieldsDict = {} #a dict indexed by key containing the info field information: NUMBER?, TYPE, DESCRIPTION 
		for line in lines:
			#print line
			if line[0] == '#':
				#read in the info fields
				if line[0:6] == '##INFO':
					parse_info_header(line, infoFieldsDict)
			else:
				splitLineFields = line.split('\t')
				if firstRecordFlag: #do special initialization of the df for the first row 
					formatTags = splitLineFields[8].split(':')
					infoDf = initialize_info_df(infoFieldsDict, formatTags)
					df = initialize_df(infoFieldsDict, formatTags)
					firstRecordFlag = False
				set_row_of_df(df, df.columns, splitLineFields, infoFieldsDict) #sets the current row 
	if onlyKeepSpecifiedColumns:
		infoDf = only_keep_specified_rows(infoDf, columnsToInclude)
		df = only_keep_specified_columns(df, columnsToInclude)
	
	#and finally write the dataframe to an xls	
	outputXlsxName = os.path.join(outputDir, udnId + '_stmpAnnotatedOutput.xlsx')	 
	writer = pd.ExcelWriter(outputXlsxName)
	print 'I am a cat'
	infoDf.to_excel(writer, 'Column Descriptions', index = False)
	df.to_excel(writer,'Sheet1', index = False)  
	writer.save()
	print 'xlsx data written to ', outputXlsxName
	return outputXlsxName

#just for testing
#alert add processing for proper gzipping etc
#vcf = '/home/noahfrie/common/udn/gateway/data/UDN781395/747009-UDN781395-P_HLYC7BCXY-2-ID07._ccP_rhP_smA_rmD_chP_final_preprocessed_fbL_testVCFANNO.vcf'
#vcf_to_xls(vcf)


