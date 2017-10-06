#written by noah friedman noahfrie
#a program aimed to provide annotations by searching the web
from py_ms_cognitive import PyMsCognitiveWebSearch
import sys
import pandas as pd
import numpy as np

bingWebSearchKey ='960301165233425c8cb9bdd8d1d718a4' #note this is a free trial key and will eventually expire or be exhausted

#mapping the column names to the function that fills them
searchFunctions = {		
 	'omim': lambda row: find_omim_link(row)		
 	#'Gene_Summary': lambda x: gene_summary_to_val(x),		
 	#'Function_Summary': lambda x: function_summary_to_val(x)
 }

def safely_convert_val_to_str(val):
	return str(val)

def find_first_correct_result_url(firstFiftyResults, correctWebpageString):
 	for result in firstFiftyResults:
 		if correctWebpageString in result.url: return result.url 
 	return 'no url found'

def find_omim_link(row):
	geneName = safely_convert_val_to_str(row['Gene Symbol']) #search based on the HGMD gene name
	searchTerm = 'omim ' + geneName
	search_service = PyMsCognitiveWebSearch(bingWebSearchKey, searchTerm)
	firstFiftyResults = search_service.search(limit=50, format='json')
	url = find_first_correct_result_url(firstFiftyResults, 'omim')
	return url

def add_cols_for_df(df, searchKeys):
	for key, function in searchFunctions.items():
		if key in searchKeys:
			df[key] = np.empty(len(df.index))
			df[key] = df[key].astype(str) #alert change based on column--we need to set columns to default to be strings

#based on an input list of search keys we create a 
def annotate_from_searches(searchKeys, xlsName):
	xls = pd.ExcelFile(xlsName)
	df = xls.parse(xls.sheet_names[0]) #this is idx 0 because our merged xls files will always have their first sheet be a key for values
	add_cols_for_df(df, searchKeys)
	for index, row in df.iterrows():
		#first iterate over all entries in the df
		#then iterate over all search functions that we are on
		for key, function in searchFunctions.items():
			if key in searchKeys: #only seach the keys that the user actually wants to search
				function = searchFunctions[key]
				value = function(row)
				df.set_value(index, key, value)

	#write the df and return it
	outputXlsxName = xlsName
	writer = pd.ExcelWriter(outputXlsxName,options={'encoding':'latin-1'})
	df.to_excel(writer,'Sheet1', index = False)
	writer.save()
	return outputXlsxName


