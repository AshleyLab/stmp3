#written by noah friedman noahfrie
#a program aimed to provide annotations by searching the web
from py_ms_cognitive import PyMsCognitiveWebSearch
import sys
import pandas as pd
import numpy as np
import urllib
import requests
import re
#alert please change!!!
sys.path.append('/home/noahfrie/noahfrie/devCode/stmp3/variant-api-client-python/') #the varsome api code
from variantapi.client import VariantAPIClient
api_key = 'pC4Kh2WCojBfs!*9SY?XcwXYdE6T*4870eZQvsXR' #THIS KEY EXPIRES ON NOVEMEBER 6th!

bingWebSearchKey ='960301165233425c8cb9bdd8d1d718a4' #note this is a free trial key and will eventually expire or be exhausted

#mapping the column names to the function that fills them
searchFunctions = {		
 	'omim': lambda row: find_omim_link(row),	
 	'rvis': lambda row: find_rvis_value(row),
 	'fathmm': lambda row: parse_varsome_for_col(row, 'fathmm_score'),
 	'mTaster': lambda row: parse_varsome_for_col(row, 'mutationtaster_converted_rankscore'),
 	'sift': lambda row: parse_varsome_for_col(row, 'sift_score'),
 	'phylop': lambda row: parse_varsome_for_col(row, 'phylop100way_vertebrate_rankscore'),
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

#corrects ref/alt for varsomes way of dealing with indels
#based on the assumption of left alignment
def correct_ref_alt_pos_for_varsome_indels(ref, alt, pos):
	if len(ref) > len(alt):
		newRef = ref[len(newRef): ]
		dif = len(ref) - len(newRef)
		position = int(pos) + dif #if ref if longer than alt we are basically deleting the front
	else:
		return 0


#gets a single column value from varsome (called multiple times)
'gnomad_exomes'
def parse_varsome_for_col(row, field):
	api = VariantAPIClient(api_key)
	#get the key for the varsome search
	chromosome = safely_convert_val_to_str(row['Chromosome'])
	position = safely_convert_val_to_str(row['Position'])
	ref = safely_convert_val_to_str(row['Reference Allele'])
	alt = safely_convert_val_to_str(row['Sample Allele'])
	#if len(ref) != len(alt): position, ref, alt = correct_ref_alt_pos_for_varsome_indels(ref, alt, position)
	varsomeLookup = 'chr' + '-'.join([chromosome, position, ref, alt])
	val = ''
	try:
		result = api.lookup(varsomeLookup, params={'add-source-databases': 'all'}, ref_genome='hg19')
		try:
			#gnomad data is formatted differently than dbnsfp data
			if 'an_' in field:
				val = result['gnomad_exomes'][0][field]
			else:
				val = result['dbnsfp'][0][field][0]
		except:
			pass
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
	#alwaysAnnotateFields = ['an_amr','an_afr','an_asj','an_eas','an_fin','an_nfe','an_oth','an_sas'] #things we always take from varsome
	alwaysAnnotateFields = []
	xls = pd.ExcelFile(xlsName)
	df = xls.parse(xls.sheet_names[0]) #this is idx 0 because our merged xls files will always have their first sheet be a key for values
	add_cols_for_df(df, searchKeys)
	cntr = 0
	for index, row in df.iterrows():
		print cntr
		#first iterate over all entries in the df
		#then iterate over all search functions that we are on
		for key, function in searchFunctions.items():
			if key in searchKeys: #only seach the keys that the user actually wants to search
				function = searchFunctions[key] 
				value = function(row)
				df.set_value(index, key, value)
		#now do the annotations we always do
		for field in alwaysAnnotateFields:
			value = parse_varsome_for_col(row, field)
			if not str(value).isdigit(): value = float('nan')  #set the value to be nan if there is no value
			df.set_value(index, field, value)

		cntr += 1

	#write the df and return it
	outputXlsxName = xlsName
	writer = pd.ExcelWriter(outputXlsxName,options={'encoding':'latin-1'})
	df.to_excel(writer,'Sheet1', index = False)
	writer.save()
	return outputXlsxName


