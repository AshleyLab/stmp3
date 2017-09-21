"""written by noah friedman
The intent of this script is to serve as a master high level pipeline that executes all analysis from sequencing file to final analysis (xls or visualization)

usage: analysis_pipeline_master_script.py pipeline_control_file.tsv
Pipeline_control_file is a tsv.  Rach row is a filed of input arguments.
The first column of each is the master name for the input arguments. Then an arbitrary number of arguments follow in each column 
Steps in the code use these arguments in order to perform pipeline tasks
A dictionary of argument abbreviations can be found below or on github

Note that higher level utility scripts will generate these control tsv files on their own

########################################################
Control.tsv specs:
each argument that is included in the filename is named with a unique three letter identifier
the first two letters are a lower case version of the argument, and the third letter is an uppercase letter identifying which class of arguments it falls under
###############Analysis direction arguments#################
line 1 preprocessing: 
	mfV (merge family vcfs)
	rhP (reheader)
	smA (split multi-allelic)
	ccP (concat snp/indels)
	chP (string chromosome prefix)
	rmD (removeDups)
line 2 calling: rgC (call rtg) 
line 3 filtering: afF (allele frequency) cvF (clinvar filtering) sgF (segregation filtering) fbL (filter by variant list)
line 4 annotation: cvA (clinvar annotation) exA (exac annotation)

#*^^^^^*#Note: if the ccP preprocessing argument is specified we will search the directory the input vcfs are in an merge them with the complementary snp or indel vcf
###############File direction arguments#########################
line 5 inputOrProbandVcf: path to the proband's vcf, or if we are running it on a database vcf, that vcf
inputOrProbandBam
inputOrProbandFastq
line 6 familyVcfs: path to vcfs for each family member.  The script will verify that they are indeed in the same family
familyBams
familyFastqs
line 7 seqFilesForRecalling: path to relevant files for recalling
line 8 variantListToFilterOn: path to list of specific variants to filter on

###############Program arguments###############################
line 9 debugArguments: saveFiles (keep all files in scratch) printSqlQueries (print sql queries)
line 10 runArchitecture: parallel (submit individual slurm jobs for subprocesses in parallel)
line 11 ped file 

"""

##############IMPORTS
import sys
import os
import subprocess
#Other scripts to import
import filter_vcf_by_variant_list
import prepare_vcfanno_conf
import write_annotated_vcf_to_xls
import merge_and_process_xls
#import segregation_util

#where this script is (useful for tools like vcf anno where we need to cd in and out of the directory)
scriptPath = '/home/noahfrie/noahfrie/devCode/stmp2/code'

#parses the control tsv specified at the outset and returns a dictionary with controls
#ALERT the code will break if the tsv doesnt have all required fields
#todo: add spot check to ensure inclusion of all required fields
def parse_control_tsv(tsvPath):
	controlParamDict = dict()
	with open(tsvPath) as f:
		lines = f.readlines()
		for line in lines:
			if line[0] != '#': #treat lines beginning with # as comments in the tsv
				vals = line.strip('\n').split('\t') #get rid of new lines and split line by tabs
				#fill the values for the control parameters dictionary
				#first remove and possible empty strings in vals[1:]
				newVals = []
				for val in vals[1:]: 
					if not val.isspace() and val != '': newVals.append(val)  
				controlParamDict[vals[0]] = newVals
	return controlParamDict

#performs a series of checks designed to ensure that we have a valid set of sequencing files with which to perform our analyses
def check_for_presenece_of_valid_seq_files(controlParamDict):
	#check 1: do we have a valid input seq file
	return 0
	"""elif: len(controlParamDict['snpAndIndelVcf']) > 0:
		if len(controlParamDict['snpAndIndelVcf']) != 2:
			sys.exit('error: you must specify precisely two files for snp and indel vcf')
		#TODO: check to ensure that snp and indel vcfs have the words snp and indel in them
		#TODO: validate that they are good vcfs
		print 'required checks are unimplemented'
		"""

	"""elif: len(controlParamDict['familyVcfs']) > 0:
		#TODO: validate each family vcf to make sure they are good vcfs
		#TODO: validate that each family vcf comes from the same family?
		print 'required checks are unimplemented'

	elif: len(controlParamDict['seqFilesForRecalling']) > 0:
		#TODO: validate these files and make sure other parameters are in concordance
		print 'required checks are unimplemented'"""


#checks to ensure that if a filtering argument is specified the required annotation is included
def check_filter_annotation_concordance(controlParamDict):
	#clinvar
	return 
	if 'cvF' in controlParamDict['filtering']:
		if 'cvA' not in controlParamDict['annotation']:
			sys.exit('error: you must perform clinvar annotations in order to filter on them')
	#frequency
	if 'afF' in controlParamDict['filtering']:
		#TODO: complete if statment: basically if not of the freq databases are annotated we cant filter on freq
		if 'exA' not in controlParamDict['annotation'] and 'wfA' not in controlParamDict['annotation']:
			sys.exit('error: necessary freq databases arent included ergo we cant do filtering')

	if len(controlParamDict['inputOrProbandVcf']) > 0:
		if len(controlParamDict['inputOrProbandVcf']) > 1:
			print controlParamDict['inputOrProbandVcf']
			sys.exit('error more than one input or proband vcf specified')
	else:
		sys.exit('no proband specified')
	#TODO add more


#function to check the coherence of specified input arguments
#a series of if statements check to make sure the arguments specified by the user are valid
def check_argument_coherence(controlParamDict):
	#check 1: do we have a valid input seq file
	check_for_presenece_of_valid_seq_files(controlParamDict)
	#print controlParamDict
	check_filter_annotation_concordance(controlParamDict)
	#TODO add more

#utility to add a suffix to a vcf file
def add_suffix_to_vcf(filename, suffix):
	return filename.replace('.vcf', '_' + suffix + '.vcf')

#utility function to get the directory of the file we are working with
def get_directory_of_file(filePath):
	directory, filename = os.path.split(filePath)
	return directory

#find the final preprocessed vcf which is the vcf with the suffix final_preprocessed
def find_final_preprocessed_vcf(dirPath):
	files = os.listdir(dirPath)
	for f in files:
		if 'final' in f and 'tbi' not in f:
			return f
	print 'error no final preprocessed vcf found'

#given the base directory of the file and either the snp or indel file it finds the other and returns both
def get_snp_and_indel_files(fileDirectory, specifiedFile):
	filesInDir = os.listdir(fileDirectory)
	snpFile = ''
	indelFile = ''
	for f in filesInDir:
		#we consider a file the snp/indel file if it contains that string and ends in 'vcf' or vcf.gz
		if 'SNP' in f and f[-3:] == 'vcf': snpFile = os.path.join(fileDirectory, f)
		if 'INDEL' in f and f[-3:] == 'vcf': indelFile = os.path.join(fileDirectory, f)
	if snpFile == '' or indelFile == '':
		print 'error no snp/indel file found'
		sys.exit()
	return snpFile, indelFile
#____________________________________________________________________________________________________________#

#ARGUMENT PARSING AND COHERENCY CHECKS###################################
controlParamDict = parse_control_tsv(sys.argv[1])
check_argument_coherence(controlParamDict) #make sure the user didn't ask us to perform an impossible pipeline
currentWorkingVcf = None
currentWorkingXls = None

#BEGIN PIPELINE#####################################################

############---------CALLING--------------##################
if len(controlParamDict['calling']) > 0:
	print 'executing calling arguments: ', controlParamDict['calling']
	#if there are arguments for calling, run calling pipelines
	#run rtg, scotch etc
##########################################################################################


#depending on whether family vcfs have been specified or not, the preprocessing either operates on a single vcf or a compendium of family vcfs
if len(controlParamDict['inputOrProbandVcf']) > 0:  #only do the following if an input vcf is specified
	vcfs = []
	vcfs.append(controlParamDict['inputOrProbandVcf'][0])
	if len(controlParamDict['familyVcfs']) > 0:
		for v in controlParamDict['familyVcfs']:
			vcfs.append(v)
	cntr = 0
	for v in vcfs: 
		fileDirectory = get_directory_of_file(v)
		#for now we store all output in a separate directory that we call 
		outputDir = os.path.join(fileDirectory, 'analysisPipelineOutput')
		#use a linux command to create the output directory
		cmd = 'mkdir "{d}"'.format(d = outputDir)
		print cmd
		#subprocess.Popen(cmd, shell=True).wait()

		#if family mode is specified vcfs = all family memebers as identified by ped file
		#else vcfs = proband only

		############---------PREPROCESSING--------------##################
		if len(controlParamDict['preprocessing']) > 0:
			print 'executing the preprocessing arguments: ', controlParamDict['preprocessing']
			
			#Step 1: check if we are in SNP/indel mode
			snpVcf = ''
			indelVcf = ''
			if 'ccP' in controlParamDict['preprocessing']:
				snpVcf, indelVcf = get_snp_and_indel_files(fileDirectory, 'ALERT PICK THIS PARAMETER TO MAKE IT VCF')

			#now we run general preprocessing dot py
			#ALERT! some of these parameters should maybe live in the pipeline controls file
			preprocessingScriptPath = '/scratch/users/noahfrie/devCode/stmp2/code/general_preprocessing.py'

			#file abbreviations from analysis_pipeline_master
			#smA (split multi-allelic)
			#chP (string chromosome prefix)
			#rhP (reheader)
			#ccP (concat snp/indels)
			#rmD (removeDups)

			#set some of the parameters for 
			splitMultiallelic = ''; reheaderVcf = ''; concat = ''; stripChrPrefix = ''; removeDups = ''
			if 'smA' in controlParamDict['preprocessing']: splitMultiallelic = 'smA'
			if 'chP' in controlParamDict['preprocessing']: stripChrPrefix = 'chP'
			if 'rhP' in controlParamDict['preprocessing']: reheaderVcf = 'rhP'
			if 'ccP' in controlParamDict['preprocessing']: concat = 'ccP'
			if 'rmD' in controlParamDict['preprocessing']: removeDups = 'rmD'
			cmd = 'python {preprocessingScript} {iVcf} {d} {sMAllelic} {sChrPrefix} {reheadVcf} {ccat} {rDups} {dIFiles} {snp} {indel}'.format(
					preprocessingScript = preprocessingScriptPath,
					iVcf = v,
					d = outputDir,
					sMAllelic = splitMultiallelic,
					sChrPrefix = stripChrPrefix,
					reheadVcf = reheaderVcf,
					ccat = concat,
					rDups = removeDups,
					dIFiles = False, #ALERT this needs to be changed
					snp = snpVcf,
					indel = indelVcf
					)
			print cmd
			subprocess.Popen(cmd, shell=True).wait()
			#finally we reset the value in our vcfs array to be the final preprocessed vcf
			vcfs[cntr] = os.path.join(fileDirectory, find_final_preprocessed_vcf(fileDirectory))
			cntr += 1

	#if we have specified mfV (merge family vcfs) we merge all vcfs in the vcfs array
	#ALERT temporary hack
	currentWorkingVcf = vcfs[0]
"""
if 'mfV':
	print 'noah'
	#/scratch/PI/euan/common/udn/gateway/data/UDN639487/682358-UDN639487-M_H2VY7BCXY-1-ID09._snp_indel.reheader_final_preprocessed.vcf.gz
	currentWorkingVcf = ''
	#TODO: figure out where the files are and what they are called etc 
	#ALERT add in code to call an
	#call merge vcf script
else:
	currentWorkingVcf = vcfs[0]"""

##########################################################################################

#ALERT! #design choice-- what directory do we want files going into if we have family mode and merging

############---------PRE-ANNOTATION FILTERING--------------##################
if len(controlParamDict['filtering']) > 0:
	if 'sgF' in controlParamDict['filtering']:
		print 'ALERT segregation filtering not set up yet'
	if 'fbL' in controlParamDict['filtering']:
		outputFileName = add_suffix_to_vcf(currentWorkingVcf, 'fbL')
		outputFileName = outputFileName.strip('.gz') #we should strip the .gz because filter vcf by variant list outputs to a uncompressed file
		'ALERT! PROBLEM BUG HERE THE OUTPUT SHOULDNT BE BGZIPED'
		#alert fix the hack on variant list to filter on
		currentWorkingVcf = filter_vcf_by_variant_list.filter_vcf_by_variant_list(controlParamDict['variantListToFilterOn'][0], currentWorkingVcf, outputFileName)
		#run segregation filtering script
		#segregation_util.filter_by_segregation(currentWorkingVcf, pedFile, outputDir, segregationModelType)
##########################################################################################

############---------ANNOTATION--------------##################
#Paths for annotation files
exacPath = '/share/PI/euan/apps/stmp/data/datasets/datasetsForVcfanno/ExAC.r0.3.1.sites.vep.vcf.gz'
caddPath = '/scratch/users/noahfrie/devCode/stmp2/vcfanno/annotationDataFiles/cadd_v1.3.vcf.gz'
gnomadPath = '/share/PI/euan/apps/stmp/data/datasets/datasetsForVcfanno/gnomad.exomes.r2.0.1.sites.vcf.gz'
clinvarPath = '/share/PI/euan/apps/stmp/data/datasets/datasetsForVcfanno/clinvar_20170905.vcf.gz'


#ALERT todo: move all annotation files to a common location instead of inside my directory
if len(controlParamDict['annotation']) > 0:
	#run annotation scripts (range annotation, point annotation etc)
	#TODO: change these to be adjusted by if/else logic
	myTestConfDict = {}
	#TODO need to sucessfully export the whole thing
	if 'cdD' in controlParamDict['annotation']: #do cadd score annotation
		myTestConfDict[caddPath] = ['raw', 'phred']
	if 'exA' in controlParamDict['annotation']: #exac
		#do all exac annotations
		myTestConfDict[exacPath] = ['KG_AF_POPMAX', 'ESP_AF_POPMAX']
	if 'gnA' in controlParamDict['annotation']: #gnomad
		myTestConfDict[gnomadPath] = ['AF_AFR', 'AF_AMR', 'AF_ASJ', 'AF_EAS', 'AF_FIN', 'AF_NFE', 'AF_OTH', 'AF_SAS']
	if 'clV' in controlParamDict['annotation']: #clinvar  alert unclear if it is working
		myTestConfDict[clinvarPath] = ['CLNSIG']
		#myTestConfDict['/scratch/users/noahfrie/devCode/stmp2/vcfanno/annotationDataFiles/common_no_known_medical_impact_20170905.vcf.gz'] = ['CLNSIG']

	#ALERT/ NOTE!  you must provide absolute paths for vcf anno to work!
	confFileName = os.path.join(os.getcwd(), 'myTestConfFile.toml')
	prepare_vcfanno_conf.write_conf_file(confFileName, myTestConfDict)
	vcfannoPath = '/home/noahfrie/noahfrie/devCode/stmp2/vcfanno/'
	outputVcfPath = add_suffix_to_vcf(currentWorkingVcf, 'testVCFANNO')
	#We need to cd into the vcfanno directory, run it, then cd back into our current directory
	os.chdir(vcfannoPath)

	cmd = './vcfanno_linux64' + ' -p 4 -lua /home/noahfrie/noahfrie/devCode/stmp2/vcfanno/example/custom.lua ' + confFileName + ' ' + currentWorkingVcf + ' > ' + outputVcfPath
	print cmd 
	subprocess.Popen(cmd, shell=True).wait()

	#go back to where we were before
	os.chdir(scriptPath)
	currentWorkingVcf = outputVcfPath
##########################################################################################


############---------POST ANNOTATION FILTERING--------------##################
#again go an run filters, but this time the ones that occur after annotation
#ALERT PLEASE CHANGe
#tiered_output_xls = stmp_tiering_util.tier_real(args, joined_outfile, yaml_commands)

if len(controlParamDict['filtering']) > 0:
	#call to stmp tiering
	#TODO allow us to specify the filtering that we do
	pass
##########################################################################################

###########-----------XLS CREATION AND PROCESSING-------------------######################
#CREATE the tiered xls--Now we work with excel sheets etc
if currentWorkingVcf != None:
	currentWorkingXls = write_annotated_vcf_to_xls.vcf_to_xls(currentWorkingVcf)

if len(controlParamDict['alreadyGeneratedXls']) > 0:  #set the current working vcf to be what the user specified if they specified something
	currentWorkingXls = controlParamDict['alreadyGeneratedXls'][0]

if len(controlParamDict['gcXls']) > 0:  #if a gc (genetic counselor) xls is included, go and perform the spreadsheet merging script
	gcXls = controlParamDict['gcXls'][0]
	merge_and_process_xls.merge_columns_across_spreadsheets(currentWorkingXls, gcXls)


print 'final working vcf: ', currentWorkingVcf
print 'stmp completed'
print 'you ran stmp with the following parameters:', controlParamDict





