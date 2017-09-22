#author: Noah Friedman
import sys

sys.path.append("/home/noahfrie/noahfrie/devCode/stmp2/code/trio")
import pedigreeUtils
import general_utils
import vcfUtils
import vcf
import os
import logging
import subprocess

#Initialization of a logger
logger = logging.getLogger("segregation")
FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
logging.basicConfig(level=logging.INFO, filemode='w', filename='segregation_util.log')
logger.setLevel(logging.INFO)
logger.info("\n\nLog File for segregation")

#a dictionary that tells you what index you need for accessing a certain argument value from the command line.  In my opinion this solution is simpler than using an argument parser especially since this is a mini utility. 
args_dict = {"ped": 0, "vcf": 1, "udnid": 2, 'modelType': 3, "dir": 4}

#------------------------------------------------------------------------------------
#UNIT TESTS

#function to ensure that the matricies are consistent with the order we ought to expect
def consistentcy_checks(familyMembersList, familyMembersMatrix, vcfMatrix, inputVcf, affectedDict, verbose=False, n = 100):
	#confirm that familyMembersMatrix is consistent with affected dict
	print "verifying that familyMembersMatrix is consistent with affected dict"
	for i in range(len(familyMembersList)):
		if verbose: print familyMembersList[i], familyMembersMatrix[i], affectedDict[familyMembersList[i]]
		if familyMembersMatrix[i] != affectedDict[familyMembersList[i]]:
			print "error inconsistency detected"
			sys.exit()
		else:
			if verbose: print "affection status for ", familyMembersList[i], " is consistent"
	print "passed"

	#confirm that order is consistent between familyMembers list and samples in the first n lines of the vcf
	print "verifying that order is consistent between familyMembers list and samples in the first ", n, " lines of the vcf"
	vcf_reader = vcf.Reader(open(inputVcf, 'r'))
	cntr = 0
	for record in vcf_reader:
		for i in range(len(record.samples)):
			if verbose: print record.samples[i].sample, familyMembersList[i]
			if record.samples[i].sample != familyMembersList[i]:
				print "error inconsistency in ordering of samples"
				sys.exit()
		cntr += 1
		if cntr >= n: break
	print "passed"

	vcf_reader2 = vcf.Reader(open(inputVcf, 'r'))
	cntr = 0
	for record in vcf_reader2:
		for i in range(len(record.samples)):
			gt = record.samples[i]['GT']
			if gt != vcfMatrix[cntr][i]: 
				print "error matrix and vcf are inconsistent"
			else:
				if verbose: print gt, vcfMatrix[cntr][i]
		cntr += 1
		if cntr >= n: break
	print "passed"

	print "passed all tests"

#------------------------------------------------------------------------------------------------

def get_arg(args, arg):
	return args[args_dict[arg]]

def get_gt_matrix_from_vcf(inputVcf):
	#iterates through a vcf and gives back the gts 
	vcf_reader = vcf.Reader(open(inputVcf, 'r'))
	m = []
	cntr = 0
	for record in vcf_reader:
		row = []
		#alert added for testing throw out
		cntr += 1
		#if cntr > 20: return m
		for sample in record.samples:
			row.append(sample['GT'])
		m.append(row)
	return m

#parse the first vcf record to figure out which family members show up in the vcf
#note: it is important that we do this instead of simply using the pedigree in case the pedigree lists family members who aren't in the vcf
def get_family_members_in_vcf(inputVcf):
	vcf_reader = vcf.Reader(open(inputVcf, 'r'))
	#just check the first record (I know this is not the prettiest method but still..)
	cntr = 0
	familyMembers = []
	for record in vcf_reader:
		for s in record.samples:
			familyMember = s.sample
			familyMembers.append(familyMember)
		cntr += 1
		if cntr >= 1:
			return familyMembers

#take the full affection status dict from trio tools and specifcially pulls out info relevant to the family in question
def get_fam_affection_status(affectedDict, familyMembers):
	tempDict = dict()
	for key, value in affectedDict.items():
		if key in familyMembers:
			tempDict[key] = value
	fMatrix = []
	#return this information as a list which conforms to the order we found in the vcf!  maintaining this order is important!
	for familyMember in familyMembers:
		fMatrix.append(tempDict[familyMember])
	return fMatrix

#IMPORTANT NOTE--our code makes the assumption that the order that the samples occur in the vcf is consistent

#checks to see if anyone in the family is missing a genotype for the spcified line
def no_missing_genotypes(gtLine):
	for gt in gtLine:
		if gt == './.': return False
	return True

#returns true if all genotypes for the current alelle in the vcf are not the same
#NOTE: IN The future this can (should) be adjusted so that it only takes into account parents and affecteds.  If we have an unaffected sibling whose genotype is different from unaffected mom and dad we will include the variant even if everyone else is the same
def all_genotypes_are_not_the_same(gtLine):
	gt1 = gtLine[0]
	#if any genotype is different than the first genotype we can conclude that all genotypes are not the same
	#as a fun exercise, the reader can prove this by induction:)
	for gt in gtLine:
		if gt != gt1:
			return True
	return False

#returns true if the genotypes of at least one of the affecteds is not 0/0 (homozygous ref)
def no_homozygous_ref_in_all_affected(gtLine, aIdxs):
	for idxA in aIdxs:
		print gtLine[idxA]
		if gtLine[idxA] != '0/0':
			return True
	return False

#returns true if there are no homozygous alternate genotypes in the unaffecteds that are composed of alternate alleles the proband has
#some examples: proband is 1/1, mom is 1/1, dad is 0/1--throw out variant cant be disease causing 
#proband is 0/1, mom is 1/1, dad is 0/0--throw out, variant cant be implicated in compound het model
#proband is 0/1, mom is 1/2, dad is 0/2, unaffected sister is 2/2--keep sister's homozygous alt isnt disase causing 
#proband is 1/2, mom is 2/2, dad is 0/1--throw out
#proband is 1/2, mom is 0/2, dad is 1/1--throw out
#proband is 1/2, mom is 0/2, dad is 0/1--KEEP
def no_homozygous_alt_in_any_unaffecteds(gtLine, uIdxs, aIdxs):
	#get the genotype of the proband
	gtProband = gtLine[aIdxs[0]]
	#isolate their alleles 
	pAllele1 = gtProband[0]
	pAllele2 = gtProband[2]
	#create a list of the 1 or 2 homozygous alts that can be formed from the probands alelles 
	pHomozygotes = []
	if pAllele1 == '1' or pAllele2 == '1': pHomozygotes.append('1/1')
	if pAllele1 == '2' or pAllele2 == '2': pHomozygotes.append('2/2')
	if pAllele1 == '3' or pAllele2 == '3': pHomozygotes.append('3/3')

	for idxU in uIdxs:
		if gtLine[idxU] in pHomozygotes:
			return False
	return True

#DEPRECATED CURRENTLY
def segregation_inconsistent_between_afffected_and_unaffected(gtLine, aIdxs, uIdxs):
	#we can't filter on lines with incomplete data
	#print gtLine, aIdxs, uIdxs
	#if len(gtLine) != len(uIdxs): return True
	for idxU in uIdxs:
		gtU = gtLine[idxU]
		for idxA in aIdxs:
			gtA = gtLine[idxA]
			if gtU == gtA:
				return False
	return True

#DEPRECATED CURRENTLY
def segregation_consistent_across_affecteds(gtLine, aIdxs):
	#we can't filter on lines with incomplete data
	if len(gtLine) != len(aIdxs): return True
	#by default segregation is consistent if there aren't two affecteds
	if len(aIdxs) < 2: return True
	else:
		gt1 = gtLine[aIdxs[0]]
		for idx in aIdxs:
			gt = gtLine[idx]
			if gt1 != gt: return False
	return True

#to be consistent with a dominant model of inheritance, all affecteds must have the alt, and not unaffecteds can have the alt
#note it permits missing records
def consistent_with_dominant_model_of_inheritance(gtLine, aIdxs, uIdxs):
	for idxA in aIdxs:
		gtA = gtLine[idxA]
		if gtA == '0/0': return False
	for idxU in uIdxs:
		gtU = gtLine[idxU]
		if gtU == '0/1' or gtU == '1/1': return False
	return True

#to be consistent with a recessive model of inheritance, all affecteds must be homozygous recessive and their parents must be heterozygous carriers
#again it is permissive with missing records
def consistent_with_recessive_model_of_inheritance(gtLine, aIdxs, uIdxs):
	#check that affecteds are homozygous alt
	for idxA in aIdxs:
		gtA = gtLine[idxA]
		if gtA == '0/1' or gtA == '0/0':
			return False
	#check that unaffecteds are heterozygous 
	for idxU in uIdxs:
		gtU = gtLine[idxU]
		if gtU == '0/0' or gtU == '1/1':
			return False
	return True


#returns the indicies of affecteds and unaffecteds in list form
def get_idxs(fMatrix):
	cntr = 0
	aIdxs = []
	uIdxs = []
	for f in fMatrix:
		if f == '2': aIdxs.append(cntr)
		elif f == '1': uIdxs.append(cntr)
		cntr += 1
	return aIdxs, uIdxs

#main function called from tiering that checks ped file to see if the pedigree is consistent
def filter_vcf_by_segregation(gtMatrix, fMatrix, inputVcf, modelType):
	aIdxs, uIdxs = get_idxs(fMatrix)
	inputVcf = general_utils.open_compressed_or_regular(inputVcf,"r")
	lines = inputVcf.readlines()
	linesToWrite = []
	rowNum = 0
	for line in lines:
		segregationFailReason = ''
		if line[0] == '#':
			print 'Line has been added because it begins with a #'
			linesToWrite.append(line)
		else:
			#by deafult we automatically include a variant if anyone in the family is missing a genotype 
			if not no_missing_genotypes(gtMatrix[rowNum]):
				#print 'Line has been added because there are missing genotypes'
				#linesToWrite.append(line)
				pass
				#Do nothing right now

			#Core Segregation Logic	
			else:
				if all_genotypes_are_not_the_same(gtMatrix[rowNum]):
					if no_homozygous_ref_in_all_affected(gtMatrix[rowNum], aIdxs):
						if no_homozygous_alt_in_any_unaffecteds(gtMatrix[rowNum], uIdxs, aIdxs):
							
							#the dominant model type requires a different test to see if it fits a dominant model of inheritance
							if modelType == 'Dominant':
								if consistent_with_dominant_model_of_inheritance(gtMatrix[rowNum], aIdxs, uIdxs):
									print 'Line has been added because it is consistent with a dominant model of inheritance'
									linesToWrite.append(line)
								else:
									segregationFailReason = 'inconsistent with a dominant model of inheritance'
							#the recessive model type requires us to do an additional test to see if it fits a recessive model of inheritance
							elif modelType == 'Recessive':
								if consistent_with_recessive_model_of_inheritance(gtMatrix[rowNum], aIdxs, uIdxs):
									print 'Line has been added because it is consistent with a recessive model of inheritance'
									linesToWrite.append(line)
								else:
									segregationFailReason = 'inconsistent with a recessive model of inheritance'									

							#deafult model is the most permissive
							elif modelType == 'Default':
								print 'Line has been added because it is consistent with the default model of inheritance'
								linesToWrite.append(line)
						
						else: 
							segregationFailReason = 'homozygous alt present in an unaffected'
							#print segregationFailReason
					else:
						segregationFailReason = 'homozygous reference present in all affecteds'
						#print segregationFailReason
				else: 
					segregationFailReason = 'All genotypes are the same'
					#print segregationFailReason

			#if segregation failed
			if segregationFailReason != '':
				logger.info('excluded line for failing segregation: ' + line)
				logger.info('failed because ' + segregationFailReason)
				logger.info('\n')
			
			rowNum += 1
	return linesToWrite

#----------------------------------------------------------------------------------
def write_vcf(linesToWrite, filename, outputDir):
	vcf_name = vcfUtils.getSampleName(filename)
	#change this please please NOAH
	outfilepath = os.path.join(outputDir, vcf_name+'_segreation_filtered_TEST_IN_JANUARY'+ '.vcf')
	f = open(outfilepath, 'w')
	for line in linesToWrite:
		f.write(line)
	f.close()

	cmd = 'bgzip -f {vcf1}'.format(vcf1=outfilepath)
	proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
	proc.wait()
    
	outfilepath += '.gz'

	cmd2 = 'bcftools index {vcf1}'.format(vcf1=outfilepath)
	proc = subprocess.Popen(cmd2, stdout=subprocess.PIPE, shell=True)
	proc.wait()

	return outfilepath

#main function for segregation filtering that is called by other scripts
def filter_by_segregation(inputVcf, inputPed, outputDir, modelType):
	#munge the ped file using the traverse function in pedigree utils
	#note that this provides a lot of extra information that we don't currently use
	#i just wanted to reuse the function instead of reinventing the wheel
	founders, parent_dict, trios, quartets, unique_trios, affectedDict = pedigreeUtils.traverse(inputPed)
	#probandUDNID = get_arg('udnid')

	familyMembersList = get_family_members_in_vcf(inputVcf)
	fMatrix = get_fam_affection_status(affectedDict, familyMembersList)
	print "generating gt_matrix"
	m = get_gt_matrix_from_vcf(inputVcf)

	#consistentcy_checks(familyMembersList, fMatrix, m, inputVcf, affectedDict)

	print "filtering by segregation"
	
	lineToWrite = filter_vcf_by_segregation(m, fMatrix, inputVcf, modelType)
	print "writing vcf"

	outputPath = write_vcf(lineToWrite, inputVcf, outputDir)
	print outputPath
	return outputPath

def main():
	print "usage segregation_util.py inputVcf udnid inputPed modelType outputDir (note it is order sensitive); model type can be Dominant, Default or Recessive"
	args = sys.argv[1:]
	inputVcf = get_arg(args, 'vcf')
	inputPed = get_arg(args, 'ped')
	print(inputPed)
	modelType = get_arg(args, 'modelType')
	outputDir = get_arg(args, 'dir')
	#NOTE THE TEMPORARY HACK BEING USED FOR FILEPATH HERE
	filter_by_segregation(inputVcf, inputPed, outputDir, modelType)

if __name__ == '__main__':
	main()




