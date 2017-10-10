#written by noah friedman noahfrie
#a program aimed to provide annotations by searching the web
from py_ms_cognitive import PyMsCognitiveWebSearch
import sys
import pandas as pd
import numpy as np
import urllib
import requests
import re

bingWebSearchKey ='960301165233425c8cb9bdd8d1d718a4' #note this is a free trial key and will eventually expire or be exhausted

#mapping the column names to the function that fills them
searchFunctions = {		
 	'omim': lambda row: find_omim_link(row),	
 	'rvis': lambda row: find_rvis_value(row),
 	'fathmm': lambda row: parse_varsome_for_col(row, 'fathmm_score'),
 	'mgi': lambda row: parse_mgi(row)
 }

def safely_convert_val_to_str(val):
	return str(val).decode('utf-8', 'ignore')
	#return str(val)

#helper function to find the first result in a search that is the correct url
def find_first_correct_result_url(firstFiftyResults, correctWebpageString):
 	for result in firstFiftyResults:
 		if correctWebpageString in result.url: return result.url 
 	return 'no url found'

 #helper function for getting raw html of a webpage
def get_raw_html(url):
	f = urllib.urlopen(url)
	return f.read()

#code that finds a regex, and cleans off the parts surrounding the target
def find_and_clean_regex(regex1, regexTarget, regex2, text):
	regexPattern = r"" + regex1 + regexTarget + regex2
	result = ''
	try:
		pattern = re.compile(r"" + regexPattern)
		matches = re.findall(pattern, text)
		if len(matches) > 0: result= matches[0].strip(regex1).strip(regex2)
	except:
		pass
	return result


def find_omim_link(row):
	geneName = safely_convert_val_to_str(row['Gene Symbol']) #search based on the ingenuity assigned gene name
	searchTerm = 'omim ' + geneName
	search_service = PyMsCognitiveWebSearch(bingWebSearchKey, searchTerm)
	firstFiftyResults = search_service.search(limit=50, format='json')
	url = find_first_correct_result_url(firstFiftyResults, 'omim')
	return url

#finds the rvis value of a gene (note this is not resilient to changes in the RVIS website)
def find_rvis_value(row):
	geneName = safely_convert_val_to_str(row['Gene Symbol']) #search based on the ingenuity assigned gene name
	#note the regex pattern is brittlely dependent on spacing for indentation
	firstHalfRegexPattern = geneName + '\n                </td>\n                <td>\n                    '
	actualValueRegexPattern = '.*?'
	secondHalfRegexPattern = '\n'
	baseRVISurl = "http://genic-intolerance.org/Search?query=" + geneName
	html = get_raw_html(baseRVISurl)
	rvisScore = find_and_clean_regex(firstHalfRegexPattern, actualValueRegexPattern, secondHalfRegexPattern, html)
	return rvisScore

#gets a single column value from varsome (called multiple times)
def parse_varsome_for_col(row, field):
	#get the key for the varsome search
	chromosome = safely_convert_val_to_str(row['Chromosome'])
	position = safely_convert_val_to_str(row['Position'])
	ref = safely_convert_val_to_str(row['Reference Allele'])
	alt = safely_convert_val_to_str(row['Sample Allele'])
	prefixVarsomeUrl = 'https://api.varsome.com/lookup/'
	suffixVarsomeUrl = '?add-all-data=1'
	varsomeUrl = prefixVarsomeUrl + ':'.join([chromosome,position,ref,alt]) + suffixVarsomeUrl #properly format the url
	response = requests.get(varsomeUrl)
	val = ''
	try:
		jsonResponse = response.json()
		val = jsonResponse['dbnsfp'][0][field][0]
	except:
		pass
	return val

#parses the MGI information for the gene in question
#we get the html for a text file associated with the phenotype report
def parse_mgi(row):
	geneName = safely_convert_val_to_str(row['Gene Symbol'])
	prefixMGIUrl = 'http://www.informatics.jax.org/allele/report.txt?phenotype=&nomen='
	suffixMGIUrl = '&chromosome=any&cm=&coordinate=&coordUnit=bp'
	returnVal = ''
	MGIUrl = prefixMGIUrl + geneName + suffixMGIUrl
	#url = 'http://www.informatics.jax.org/allele/report.txt?phenotype=&nomen=TKT&chromosome=any&cm=&coordinate=&coordUnit=bp'
	html = get_raw_html(MGIUrl)
	#parse this text file to get what I want
	phenotypeIndex = 8
	for line in html.split('\n')[1:]: #skip the header line
		splitByTab = line.split('\t')
		if len(splitByTab) > 8:
			returnVal += splitByTab[8] #the value I want is at line 8
	return returnVal

def add_cols_for_df(df, searchKeys):
	for key, function in searchFunctions.items():
		if key in searchKeys:
			df[key] = np.empty(len(df.index))
			df[key] = df[key].astype(str) #alert change based on column--we need to set columns to default to be strings

#use this as a global variable so we dont have to do edge cases for how functions from our dictionary are set up
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


