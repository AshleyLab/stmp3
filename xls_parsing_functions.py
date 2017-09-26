#Noah Friedman noahfrie
#provides a class of xls parsing functions for various other scripts 

import pandas as pd

#iterates through all the sheets of a b excel file and 
def read_xls_sheets(xlsName):
	xls = pd.ExcelFile(xlsName)
	sheetNames = xls.sheet_names
	sheetDict = dict()
	for sheet in sheetNames:
		sheetDict[sheet] = xls.parse(sheet)
	return sheetDict

#parses what the GC's wrote for IGV and returns the proper class label
#two modes: passFailMode where the labels are just pass or fail
#multi class mode where there are more labels
def parse_igv_comment(comment, passFailMode = True):
	if passFailMode:
		if 'pass' in comment:
			return 'PASS'
		elif 'fail' in comment:
			return 'FAIL'
		else:
			#if the words pass or fail arent in the comment we cant the variant so we return a sentinel
			return -1
	else:
		print 'Warning other modes not supported yet'

#method for getting xls values with error checking for missing values
def get_xls_value(row, key, missingVal = 'value missing'):
	if key in row:
		return str(row[key])
	else:
		return missingVal

#automatically tells us what to exclude from our 
def exclude_sheets_based_on_missing_columns(sheetDict):
	sheetsWithMissingColumns = []
	for key, value in sheetDict.items():
		curDf = sheetDict[key]
		#print curDf.columns
		if 'CHROM' not in curDf.columns:
			sheetsWithMissingColumns.append(key)
	#print 'king rhoam'
	#print sheetsWithMissingColumns
	return sheetsWithMissingColumns


