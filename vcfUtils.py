#Rick Dewey 5.22.12, Prag Batra 2015
#module for general vcf utilities
#usage: utils module for vcf parsing
#!/usr/bin/env python

import os
import re
import sys
import vcfHeaders
import yaml_keys
import yaml_utils
import subprocess
import logging

import sys
import linecache
import sqlite3
import os
from subprocess import *
import subprocess
import mmap
import datetime
import gzip
import multiprocessing
from multiprocessing import Pool

# third-party modules
import vcf
import yaml

# our internal imports
import yaml_keys
import stmp_consts
import vcfHeaders
import vcfUtils
import tsv_utils
import yaml_utils
#import general_utils
#import expand_intersectedBeds
#import stmp_annotation_checker
#import db_utils
#import logging_utils
#import logging
#import snpeff_annotation
#import annovar_annotation
#import stmp_range_annotation
#import point_annotation
#import file_join_util
#import stmp_annotation_util


# gets name of the headers for the sample columns (in the order in which they are listed in the VCF)
# returns array with 1 header at each index
def get_sample_headers(vcf_file_loc):
    logger = logging.getLogger(__name__)
    cmd = 'bcftools query -l {vcf}'.format(vcf=vcf_file_loc)
    logger.debug('cmd to get vcf sample col headers: ' + str(cmd))
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    output=proc.stdout.read()
    return output.split('\n')

def getNumSampleCols(vcf_file_loc):
    logger = logging.getLogger(__name__)
    cmd = 'bcftools query -l {vcf}|wc -l'.format(vcf=vcf_file_loc)
    logger.debug('cmd to get # of vcf sample cols: ' + str(cmd))
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    output = proc.stdout.read()
    return int(output)

#gets proper head from vcf or tsv file, ignoring all other stuff
def get_head(infile):
    logger = logging.getLogger(__name__)
    f1 = open(infile, "r")
    if(not infile.endswith('.vcf')):
        head = f1.readline().replace('#', '')
        f1.close()
        return head
    #else
    while True:
        line = f1.readline()
        if "#CHROM" in line:
            head = line.replace("#", "")
            break
        if not line:
            logger.error("Error in vcfUtils.get_head - End of file reached prior to finding vcf header: vcf header is likely malformed", exc_info=True)
            raise ValueError("Error in vcfUtils.get_head - End of file reached prior to finding vcf header: vcf header is likely malformed")
#             print >> sys.stderr, "Error in vcfUtils.get_head - End of file reached prior to finding vcf header: vcf header is likely malformed"
#             exit(1)
    f1.close()
    return head

#gets indexes from head of annotation file for downstream filtering/annotation
def get_listindex(head_list, viplist):
    indexlist = []
    tempindex = 0
    for item in head_list:
        if item in viplist:
            indexlist.append(tempindex)
        tempindex+=1
    return indexlist

# (from trio vcfUtils)
#parses format field of vcf file to find location of AVR, DP, RQ
def parse_format(formatstats):
    logger = logging.getLogger(__name__)
    temp2 = formatstats.split(':')
    counter = 0
    avr_num = 0 
    dp_num = 0
    rq_num = 0
    gq_num = 0
    found_avr = False
    found_dp = False
    found_rq = False
    found_gq = False
    for infos in temp2:
        if (str(infos) == "AVR"): 
            avr_num = int(counter)
            found_avr = True
        if (str(infos) == "DP"): 
            dp_num = int(counter)
            found_dp = True
        if (str(infos) == "RQ"): 
            rq_num = int(counter)
            found_rq = True
        if (str(infos) == "GQ"): 
            gq_num = int(counter)
            found_gq = True
        counter+=1
    if(not found_avr):
         logger.warn('warning could not find AVR in FORMAT stats ' + str(formatstats))
        #None
    if(not found_dp):
         logger.warn('warning could not find DP in FORMAT stats ' + str(formatstats))
        #None
    if not found_rq:
         logger.warn('warning could not find RQ in FORMAT stats ' + str(formatstats))
        #None
    if not found_gq:
         logger.warn('warning could not find GQ in FORMAT stats ' + str(formatstats))
        #None
    return avr_num, dp_num, rq_num, gq_num, found_avr, found_dp, found_rq, found_gq

# (from trio vcfUtils)
#parses vcf format genotypes into allele codes, works for bialleleic positions, need to update for multiallelic positions, works for GATK standard vcf format only so far
def allele_coder(GT, d, alt, q, p, avr_n, dp_n, rq_n, gq_n, data_type):
    logger = logging.getLogger(__name__)
    logger.debug('GT: ' + str(GT))
    
    gt_list = GT.split(':')
    genotype = gt_list[0]
    if(genotype == '0/1' or genotype == '0|1' or genotype == '1|0'):
        logger.debug('found 0/1 gt in list ' + str(GT))
        alleles = '1'
    elif(genotype == '1/1' or genotype == '1|1'):
        logger.debug('found 1/1 gt in list ' + str(GT))
        alleles = '2'
    elif genotype == '0/0' or genotype == '0|0':
        logger.debug('found 0/0 gt in list ' + str(GT))
        alleles = '0'
    else:
        logger.warn('warning: unknown genotype ' + str(genotype) + ' FORMAT line ' + str(GT)) 
        alleles = "4"
    return alleles
    # TODO implement proper filtering criteria to check validity of allele calls
#     #3 for alleles not meeting filtering criteria
#     else:
#         alleles = '3'
#     #4 for alleles not called
#     else:
#         alleles = "4"
    
    #deals with illumina format vcf files
    if (data_type == "ILLUMINA") or (data_type == "ILL"):
        if ("./." in GT) == 0:
            gt_list = GT.split(":")
            likelihoodvec = gt_list[4].split(",")
            logger.debug("\nlike vec "+str(likelihoodvec))
            gq = float(gt_list[3])
            temp = gt_list[1].split(",")
            if not ".:." in GT:
            #if int(temp[1]) != ".":
                al = int(temp[1])
                totalDepth = int(temp[0])+int(temp[1])
            # handles imputed genotypes from GATK phase by transmission tool which may impute genotypes without read information
            if ".:." in GT:
                al = 0
                totalDepth = 0
            
            if (gt_list[0] == "0/1" or gt_list[0] == "0|1" and totalDepth >= d and gq >= q and al >= alt and int(likelihoodvec[1]) <= p):
                alleles = "1"
            elif (gt_list[0] == "1|0" and totalDepth >= d and gq >= q and int(temp[0]) >= alt and int(likelihoodvec[0]) <= p):
                alleles = "1"
            elif (gt_list[0] == "1/1" or gt_list[0] == "1|1" and totalDepth >= d and gq >= q and al >= alt and int(likelihoodvec[2]) <= p):
                alleles = "2"
            elif (gt_list[0] == "0/0" or gt_list[0] == "0|0" and totalDepth >= d and gq >= q and int(likelihoodvec[0]) <= p):
                alleles = "0"

            #3 for alleles not meeting filtering criteria
            else:
                alleles = "3"

        #4 for alleles not called
        else:
            alleles = "4"
        return alleles

    #deals with complete genomics format vcf files
    if (data_type == "COMPLETE") or (data_type == "CG"):
        if ("./." in GT) == 0:
            gt_list = GT.split(":")
            if gt_list[0] == "0/1" or gt_list[0] == "0|1" or gt_list[0] == "1|0":
                alleles = "1"
            elif gt_list[0] == "1/1" or gt_list[0] == "1|1":
                alleles = "2"
            elif gt_list[0] == "0/0" or gt_list[0] == "0|0":
                alleles = "0"
            else:
                alleles = "3"

        #4 for alleles not called
        else:
            alleles = "4"
        return alleles
        
    #deals with RTG vcf files
    #AVR definition: adaptive variant rescoring value. It is a value between 0 and 1 that represents the probability that the variant for the sample is correct.
    #In the RTG version of the allele coder, likelihoodvec stores the AVR value, the --post_hoc_filter_cutoff represents the lower limit for inclusion
    #Cutoffs occur if the likelihoodvec is < p or pl
    if (data_type == "RTG") or (data_type == "rtg"):
        logger.debug("\nGT "+str(GT))
        if ("./." in GT) == 0:
            gt_list = GT.split(":")
            try:
                likelihoodvec = float(gt_list[avr_n])
            except ValueError:
                logger.debug('gt list: ' + str(gt_list))
                logger.debug('avr_n: ' + str(avr_n))
                logger.debug('gt_list[avr_n]: ' + str(gt_list[avr_n]))
                raise
            gq = float(gt_list[gq_n])
            totalDepth = int(gt_list[dp_n])
            logger.debug("\n p "+str(p)," AVR "+str(likelihoodvec))
            if (gt_list[0] == "1/1" or gt_list[0] == "1|1"):
                al = totalDepth * 0.5
            if totalDepth != 0:
                #In RTG formatted vcf files, the alternate allele depth not reported in the stats field
                al = int(99)
                totalDepth = int(temp[0])+int(temp[1])
                logger.debug("\ntotaldepth "+str(totalDepth), "\nal "+str(al), "\nlikelihoodvec "+str(likelihoodvec), "\ngq "+str(gq))
            # handles imputed genotypes from RTG which may impute genotypes without read information
            if totalDepth == 0:
                al = 0
                logger.debug("\ntotaldepth "+str(totalDepth), "\nal "+str(al), "\nlikelihoodvec "+str(likelihoodvec), "\ngq "+str(gq))
            
            if (gt_list[0] == "0/1" or gt_list[0] == "0|1" and totalDepth >= d and gq >= q and al >= alt and likelihoodvec >= p):
                alleles = "1"
                logger.debug("\ntotaldepth "+str(totalDepth), "\nal "+str(al), "\nlikelihoodvec "+str(likelihoodvec), "\ngq "+str(gq))
            elif (gt_list[0] == "1|0" and totalDepth >= d and gq >= q and al >= alt and likelihoodvec >= p):
                alleles = "1"
                logger.debug("\ntotaldepth "+str(totalDepth), "\nal "+str(al), "\nlikelihoodvec "+ str(likelihoodvec), "\ngq "+str(gq))
            elif (gt_list[0] == "1/1" or gt_list[0] == "1|1" and totalDepth >= d and gq >= q and al >= alt and likelihoodvec >= p):
                alleles = "2"
                logger.debug("\ntotaldepth "+str(totalDepth), "\nal "+str(al), "\nlikelihoodvec "+str(likelihoodvec), "\ngq "+str(gq))
            elif (gt_list[0] == "0/0" or gt_list[0] == "0|0" and totalDepth >= d and gq >= q and likelihoodvec >= p):
                logger.debug("\ntotaldepth "+str(totalDepth), "\nal "+str(al), "\nlikelihoodvec "+str(likelihoodvec), "\ngq "+str(gq))
                alleles = "0"
                logger.debug("\ntotaldepth "+str(totalDepth), "\nal "+str(al), "\nlikelihoodvec "+str(likelihoodvec), "\ngq "+str(gq))
            #3 for alleles not meeting filtering criteria
            else:
                alleles = "3"

        #4 for alleles not called
        else:
            alleles = "4"
        return alleles

    else:
        print >> sys.stderr, "Invalid input - unknown or unsupported data format: "+data_type
        exit(1)

#modeled off of "is_rare"--returns 0 for a variant that is not the reference, returns 1 for 0/0 or ./. variables---also it has sentimnetal value for noah because its the first real function he has written in stmp
def is_non_ref(templine, bp_indexlist):
    l = templine.split('\t')
    gt = l[bp_indexlist[0]]
    if gt != '0/0' and gt != './.':
        return 1
    else: return 0

#boolean, returns 1 for a variant with frequency below threshold (for all background populations consider) or novel or unknown allele frequency, 0 otherwise 
def is_rare(templine, freq_thresh, bp_indexlist, yaml_commands):
    logger = logging.getLogger(__name__)
    logger.debug('bp_indexlist: ' + str(bp_indexlist))
    multimatch_delimiter = yaml_utils.get_dataset_defaults(yaml_commands)[yaml_keys.kDMultimatchDelimiter]
    templinelist = templine.split("\t")
    rare_flag = 0 # not rare
    #note the logic has been changed so that we only keep a variant if there is proof that it is rare 
    for i in bp_indexlist:
        if templinelist[i] != "":
            templinelistElts = templinelist[i].split(multimatch_delimiter)
            # should throw an exception if templinelist[i] isn't a float
            logger.debug('templinelist[{i}]: '.format(i=i) + str(float(templinelist[i])))
            logger.debug('freq thresh: ' + str(float(freq_thresh)))
            for elt in templinelistElts:
                if(float(elt) < float(freq_thresh)):
                    rare_flag = 1 # rare
    return rare_flag

# returns true if the specified text is found in any of the specified columns
# NOTE: columns must already be converted using yaml_utils (so they need to be absolute column names in the annotated output file)
def contains_text(text, templine, columns, header, yaml_commands, case_sensitive=False):
    lineSubset = ''
    lineContents = templine.rstrip("\n").split("\t")
    for column in columns:
        lineSubset += ' ' + lineContents[header.index(column)]
    
    if((case_sensitive and text in lineSubset) or (not case_sensitive and text.lower() in lineSubset.lower())):
        return True
    #else
    return False


#boolean, returns 1 for a variant with user-defined functional properties as annotated by annovar
def is_functional(templine, function_list, functional_columns, header):
    lineFunctionalSubset = ''
    lineContents = templine.rstrip("\n").split("\t")
    for column in functional_columns:
        lineFunctionalSubset += ' ' + lineContents[header.index(column)]
    
    templine = lineFunctionalSubset
    
    functional = 0
    if ("stoploss" in templine) and ("stoploss" in function_list):
        functional=1
    elif ("stopgain" in templine) and ("stopgain" in function_list):
        functional=1
    elif ("splicing" in templine) and ("splicing" in function_list):
        functional=1
    elif ("frameshift" in templine) and (("nonframeshift" in templine) == 0) and ("frameshift" in function_list):
        functional=1
    elif ("nonframeshift" in templine) and ("nonframeshift" in function_list):
        functional=1
    elif ("nonsynonymous" in templine) and ("nonsynonymous" in function_list):
        functional=1
    elif ("exonic" in templine) and ("exonic" in function_list):
        functional=1
    elif ("intronic" in templine) and ("intronic" in function_list):
        functional=1
    elif ("UTR5" in templine) and ("UTR5" in function_list):
        functional=1
    elif ("UTR3" in templine) and ("UTR3" in function_list):
        functional=1
    elif ("ncRNA" in templine) and ("ncRNA" in function_list):
        functional=1
    elif ("upstream" in templine) and ("upstream" in function_list):
        functional=1
    elif ("intergenic" in templine) and ("intergenic" in function_list):
        functional=1
    elif ("downstream" in templine) and ("downstream" in function_list):
        functional=1
    return functional

# returns true if the given line passes ExAC tolerance thresholds (based on parameters specified in YAML)
# NOTE: currently, we check all specified tolerance columns to see if ANY of their z-scores are above the specified z-score cutoff, returning true if this is the case.
def tolerance_pass(line, headlist, yaml_commands):
    lineList = line.rstrip("\n").split("\t")
    tiering_cmds = yaml_commands[yaml_keys.kModules][yaml_keys.kTiering]
    tolerance_zscore_colHeaders = yaml_utils.convertColumns(tiering_cmds[yaml_keys.kTToleranceZScoreCols], yaml_commands)
    tolerance_zscore_cutoff = tiering_cmds[yaml_keys.kTToleranceZScoreCutoff]
    
    for zscore_col_header in tolerance_zscore_colHeaders:
        print lineList, zscore_col_header, headlist
        zscores = lineList[headlist.index(zscore_col_header)].split(yaml_utils.get_dataset_defaults(yaml_commands)[yaml_keys.kDMultimatchDelimiter]) # our delimiter
        for zscore in zscores:
            if(zscore != '' and float(zscore) > float(tolerance_zscore_cutoff)):
                return True
    #else
    return False

#boolean, returns 1 for variant that is conserved according to user-defined criteria
def is_conserved(templine, headlist, yaml_commands):
    logger = logging.getLogger(__name__)
    total = 0
    templinelist = templine.rstrip("\n").split("\t")
    
    tiering_cmds = yaml_commands[yaml_keys.kModules][yaml_keys.kTiering]
    colHeaders = yaml_utils.convertColumns(tiering_cmds[yaml_keys.kTConservationCols], yaml_commands)
    colThresholds = tiering_cmds[yaml_keys.kTConservationThresholds]
    thresh = tiering_cmds[yaml_keys.kTConservationGlobalThreshold]
    
    for idx,colHeader in enumerate(colHeaders):
        colThreshold = colThresholds[idx]
        col = templinelist[headlist.index(colHeader)]
        try:
            if(col != '' and
               (((type(colThreshold) is float or type(colThreshold) is int) and float(col) >= colThreshold)
                or col == colThreshold)
            ):
                total += 1
        except ValueError:
            logger.debug('headlist ' + str(headlist))
            logger.debug('templinelist ' + str(templinelist))
            logger.debug('colHeader: ' + str(colHeader))
            logger.debug('headlist index: ' + str(headlist.index(colHeader)))
            logger.debug('col: ' + str(col))
            raise
        #debug
        #elif(col == ''):
         #   logger.warn('warning: no value for ')

    if total >= thresh:
        return 1 # True
    else:
        return 0 # False


def getClinvarInfoCol(templine, headlist):
    logger = logging.getLogger(__name__)
    templinelist = templine.split("\t")
    clinvarInfo = templinelist[headlist.index(vcfHeaders.kClinvar)] if vcfHeaders.kClinvar in headlist else ''
    if(vcfHeaders.kClinvar not in headlist):
        logger.warn('warning: ' + vcfHeaders.kClinvar + ' not found in annotated header')
    return clinvarInfo

# TODO be consistent -- when doing clinvar star annotation, need clinvar clinrevstatus without "clinvar_" prefix, but when doing tiering with final annotated VCF, need clinvar clinsig with "clinvar_" prefix
def getClinvarClinsigHeader(yaml_cmds):
    return yaml_utils.get_datasets(yaml_cmds)['clinvar'][yaml_keys.kDAnnotation]+'_'+vcfHeaders.kClinvarClinSigHeader
#     return vcfHeaders.kClinvarClinSigHeader
def getClinvarClinRevStatusHeader(yaml_cmds):
#     return yaml_cmds['clinvar']['Annotation']+'_'+vcfHeaders.kClinvarClinRevStatusHeader
    return vcfHeaders.kClinvarClinRevStatusHeader

# return clinsig value(s) as a string (may have multiple reported values in clinvar)
def getClinvarClinsigStr(templine, headlist, yaml_cmds):
    logger = logging.getLogger(__name__)
    clinvarClinsigHeader = getClinvarClinsigHeader(yaml_cmds)
    if(clinvarClinsigHeader in headlist):
        idx = headlist.index(clinvarClinsigHeader)
        tempLineElts = templine.rstrip("\n").split("\t")
        return tempLineElts[idx]
    else:
        logger.error('could not get index of clinvar clinical significance column (' + clinvarClinsigHeader + ') in annotated file')
        raise IndexError('could not get index of clinvar clinical significance column (' + clinvarClinsigHeader + ') in annotated file')
        #print 'error could not get index of clinvar clinical review status column in annotated file'
        return ''

def clinvarClinsigStr2Array(clinsig_str):
    return clinsig_str.split('|')
def getClinvarClinsigArray(templine, headlist):
    return clinvarClinsigStr2Array(getClinvarClinsigStr(templine, headlist))

def isClinvarPathogenicOrLikelyPathogenic(line, headlist, yaml_cmds):
    clinsigStr = getClinvarClinsigStr(line, headlist, yaml_cmds)
    if(str(vcfHeaders.kCLINVAR_CLINSIG_LIKELY_PATHOGENIC_CODE) in clinsigStr or (str(vcfHeaders.kCLINVAR_CLINSIG_PATHOGENIC_CODE) in clinsigStr 
#         and str(vcfHeaders.kCLINVAR_CLINSIG_OTHER_CODE) != clinsigStr
        )):
        return True
    return False

#return clinrevstatus values as a string (similar to clinsig above)
def getClinvarClinReviewStatusStr(templine, headlist, yaml_cmds):
    clinvarClinRevStatusStr = getClinvarClinRevStatusHeader(yaml_cmds)
    if(clinvarClinRevStatusStr in headlist):
        idx = headlist.index(clinvarClinRevStatusStr)
        tempLineElts = templine.rstrip("\n").split("\t")
        return tempLineElts[idx]
    else:
        raise IndexError('could not get index of clinvar clinical review status column (' + clinvarClinRevStatusStr +') in annotated file')
        logger.error('error could not get index of clinvar clinical review status column in annotated file')
        return ''
    
# clinrevstatus convenience methods
def clinvarClinRevStatusStr2Array(clinrevstat_str):
    return clinrevstat_str.split('|')
def getClinvarClinReviewStatusArray(templine, headlist):
    return clinvarClinRevStatusStr2Array(getClinvarClinReviewStatusStr(templine, headlist))


# computes clinvar stars based on clinsig + clin review status (based on revised ClinVar mid 2015 criteria)
# from http://www.ncbi.nlm.nih.gov/clinvar/docs/variation_report/#review_status
# 0 stars: No submitter provided an interpretation with assertion criteria (no assertion criteria provided), or no interpretation was provided (no assertion provided)
# 1 star: At least one submitter provided an interpretation with assertion criteria (criteria provided, single submitter) or multiple submitters provided assertion criteria but there are conflicting interpretations in which case the independent values are enumerated for clinical significance (criteria provided, conflicting interpretations)
# 2 stars: Two or more submitters providing assertion criteria provided the same interpretation (criteria provided, multiple submitters, no conflicts)
# 3 stars: reviewed by expert panel (http://www.ncbi.nlm.nih.gov/clinvar/docs/review_guidelines/)
# 4 stars: practice guideline (http://www.ncbi.nlm.nih.gov/clinvar/docs/review_guidelines/)
def clinvarStars(templine, headlist, yaml_cmds):
    clinRevStatStr = getClinvarClinReviewStatusStr(templine, headlist, yaml_cmds)
    # new star code (for clinvar xml)
    if('practice guideline' in clinRevStatStr):
        return 4
    if('reviewed by expert panel' in clinRevStatStr):
        return 3
    if('criteria provided, multiple submitters, no conflicts' in clinRevStatStr):
        return 2
    if('criteria provided, single submitter' in clinRevStatStr or 'criteria provided, conflicting interpretations' in clinRevStatStr):
        return 1
    if('no assertion criteria provided' in clinRevStatStr or 'no assertion provided' in clinRevStatStr):
        return 0
    else:
        return '' # unknown
    

# checks whether the given VCF line passes the filtering criteria
def filter_pass(yaml_commands, line, vcf_headers):
    filter_pass_values = yaml_commands[yaml_keys.kModules][yaml_keys.kMDefaults][yaml_keys.kMFilterPassValues]
    line_contents = line.split('\t')
    filter_col_text = line_contents[vcf_headers.index("FILTER")]
    for value in filter_pass_values:
        if value in filter_col_text:
            return True
    return False

# whether the given column passes the given criterion (> if a given threshold = number, or = if threshold = string)
def passes_criteria(col, colThresholds):
    for colThreshold in colThresholds:
        if(col != '' and
           (((type(colThreshold) is float or type(colThreshold) is int) and float(col) >= colThreshold)
            or col == colThreshold)
        ):
            return True
    #else
    return False


#boolean, returns 1 for variant that is pathogenic according to user-defined criteria
def is_pathogenic(templine, headlist, yaml_commands):
    templinelist = templine.split("\t")
    pathogenic = 0
    
    tiering_cmds = yaml_commands[yaml_keys.kModules][yaml_keys.kTiering]
    nalg = tiering_cmds[yaml_keys.kTPathogenicityGlobalThreshold]
    
    colHeaders = yaml_utils.convertColumns(tiering_cmds[yaml_keys.kTPathogenicityCols], yaml_commands)
    colsThresholds = tiering_cmds[yaml_keys.kTPathogenicityThresholds]
    
    for idx,colHeader in enumerate(colHeaders):
        if(isinstance(colHeader, list)):
            passed = False
            colThresholdsList = colsThresholds[idx]
            for idx,singleColHeader in enumerate(colHeader):
                colThresholds = yaml_utils.split_multiple_col_thresholds(colThresholdsList[idx], yaml_commands)
                col = templinelist[headlist.index(singleColHeader)]
                if(passes_criteria(col, colThresholds)):
                    passed = True
                    break
            if(passed):
                pathogenic += 1
            continue

        #else
        colThresholds = yaml_utils.split_multiple_col_thresholds(colsThresholds[idx], yaml_commands)
        col = templinelist[headlist.index(colHeader)]
        passed = passes_criteria(col, colThresholds)
        if(passed):
            pathogenic += 1
            
    if pathogenic >= int(nalg):
        #debug
#         print 'is pathogenic'
        return 1
    else:
        #debug
#         print 'not pathogenic'
        return 0


#CURRENTLY UNUSED
#finds variants with regulatory annotations
def is_regulatory(templine, headlist):
    templinelist = templine.split("\t")
    regulatory = 0
    dnase = 0
    adult_enhancer = 0
    fetal_enhancer = 0
    tfbs = 0
    transfac=0
    mirna_coding=0
    mirna_target=0
    if templinelist[headlist.index(vcfHeaders.kEncode_dnaseCluster)] != "":
        if int(templinelist[headlist.index(vcfHeaders.kEncode_dnaseCluster)]) >= 300:
            regulatory=1
            dnase = 1
    if templinelist[headlist.index("encode_tfbsChip_score")] != "": # warning not in annotation
        if int(templinelist[headlist.index("encode_tfbsChip_score")]) >= 300:
            regulatory=1
            tfbs=1
    if templinelist[headlist.index("heart_adult_enhancer")] != "": # warning not in annotation
        regulatory=1
        adult_enhancer=1
    if templinelist[headlist.index("heart_fetal_enhancer")] != "": # warning not in annotation
        regulatory=1
        fetal_enhancer=1
    if templinelist[headlist.index("transFac_score")] != "": # warning not in annotation
        regulatory=1
        transfac=1
    if templinelist[headlist.index("miRNA")] != "": # warning not in annotation
        regulatory=1
        mirna_coding=1
    if templinelist[headlist.index("targetScan")] != "": # warning not in annotation
        regulatory=1
        mirna_target=1
    return regulatory, dnase, tfbs, adult_enhancer, fetal_enhancer, transfac, mirna_coding, mirna_target
    

# CURRENTLY UNUSED
#finds variants in topologically central location in network according to user defined criteria
def is_central_mod(templine, headlist, rank_thresh, phenotype):
    templinelist = templine.split("\t")
    if phenotype == "normal":
        if templinelist[headlist.index("normal_heart_Kin")] != "": # warning not in annotation
            if float(templinelist[headlist.index("normal_heart_Kin")]) >= rank_thresh:
                return 1
            else:
                return 0
        else:
            return 0
    if phenotype == "hypertrophy":
        if templinelist[headlist.index("hypertrophy_heart_Kin")] != "": # warning not in annotation
            if float(templinelist[headlist.index("hypertrophy_heart_Kin")]) >= rank_thresh:
                return 1
            else:
                return 0
        else:
            return 0
    if phenotype == "failure":
        if templinelist[headlist.index("failure_heart_Kin")] != "": # warning not in annotation
            if float(templinelist[headlist.index("failure_heart_Kin")]) >= rank_thresh:
                return 1
            else:
                return 0
        else:
            return 0


# CURRENTLY UNUSED
#finds variants in topologically central location in network according to user defined criteria
def is_central_all(templine, headlist, rank_thresh, phenotype):
    templinelist = templine.split("\t")
    if phenotype == "normal":
        if templinelist[headlist.index("normal_heart_globalK")] != "":
            if float(templinelist[headlist.index("normal_heart_globalK")]) >= rank_thresh:
                return 1
            else:
                return 0
        else:
            return 0
    if phenotype == "hypertrophy":
        if templinelist[headlist.index("hypertrophy_heart_globalK")] != "":
            if float(templinelist[headlist.index("hypertrophy_heart_globalK")]) >= rank_thresh:
                return 1
            else:
                return 0
        else:
            return 0
    if phenotype == "failure":
        if templinelist[headlist.index("failure_heart_globalK")] != "":
            if float(templinelist[headlist.index("failure_heart_globalK")]) >= rank_thresh:
                return 1
            else:
                return 0
        else:
            return 0


# CURRENTLY UNUSED
#finds variants in expressed regions according to user-defined criteria
def is_expressed(templine, headlist, rank_thresh, phenotype):
    templinelist = templine.split("\t")
    if phenotype == "normal":
        if templinelist[headlist.index("normal_heart_exprs_rank")] != "":
            if float(templinelist[headlist.index("normal_heart_exprs_rank")]) >= rank_thresh:
                return 1
            else:
                return 0
        else:
            return 0
    if phenotype == "hypertrophy":
        if templinelist[headlist.index("hypertrophy_heart_exprs_rank")] != "":
            if float(templinelist[headlist.index("hypertrophy_heart_exprs_rank")]) >= rank_thresh:
                return 1
            else:
                return 0
        else:
            return 0
    if phenotype == "failure":
        if templinelist[headlist.index("failure_heart_exprs_rank")] != "":
            if float(templinelist[headlist.index("failure_heart_exprs_rank")]) >= rank_thresh:
                return 1
            else:
                return 0
        else:
            return 0


# CURRENTLY UNUSED
#finds variants in differentially expressed regions according to user-defined criteria
def is_diffexpr(templine, headlist, phenotype, q_thresh):
    templinelist = templine.split("\t")
    if phenotype == "hypertrophy":
        if templinelist[headlist.index("sam_q_normal_hypertrophy")] != "":
            if float(templinelist[headlist.index("sam_q_normal_hypertrophy")]) <= q_thresh:
                return 1
            else:
                return 0
        else:
            return 0
    if phenotype == "failure":
        if templinelist[headlist.index("sam_q_normal_failure")] != "":
            if float(templinelist[headlist.index("sam_q_normal_failure")]) <= q_thresh:
                return 1
            else:
                return 0
        else:
            return 0


#gives max allele frequency in a list of allele frequency annotations
def max_af(templine, headlist, bp_indexlist):
    templinelist = templine.split("\t")
    af = 0.00
    for i in bp_indexlist:
        if templinelist[i] != "":
            if float(templinelist[i]) > af:
                af = float(templinelist[i])
    return af
                
# (originally in vcfUtils)
# #parses info field of vcf file and returns tuple float for mq and mq0
# def parse_info(infofield):
#     infolist = infofield.split(";")
#     for element in infolist:
#         if "MQ=" in element:
#             mapq=float(element.replace("MQ=", ""))
#         elif "MQ0=" in element:
#             mapq0=float(element.replace("MQ0=", ""))
#     return mapq, mapq0

# (from trio vcfUtils)
#parses info field of vcf file and returns tuple float for mq and mq0
def parse_info(infofield):
    infolist = infofield.split(";")
    #print "\nutils before infofield "+str(infofield)
    mapq = 0.0
    mapq0 = 0.0
    #print "\nutils before mapq, mapq0 "+str(mapq)+" "+str(mapq0)
    for element in infolist:
        if "MQ=" in element:
            mapq=float(element.replace("MQ=", ""))
        elif "MQ0=" in element:
            mapq0=float(element.replace("MQ0=", ""))
    #print "\nutils after infofield "+str(infofield)
    #print "\nutils mapq, mapq0 "+str(mapq)+" "+str(mapq0)
    return mapq, mapq0

#merges vcf files split by chromosome, writing head from chromosome 1 to X only for now
def mergeFiles(fin_stem, f_out):
    chrom_arr = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"]

    fout = open(f_out, "w")

    for chrom in chrom_arr:
        fin = open(fin_stem+"_chr"+chrom+"_annotated_hmmfiltered.txt", "r")
        head = fin.readline()
        if chrom == "1":
            fout.write(head)

        for line in fin.readlines():
            fout.write(line)
        fin.close()
        os.system("rm "+fin_stem+"_chr"+chrom+"_annotated_hmmfiltered.txt")
    fout.close()


#splits vcf file by chromosome, to X only for now
def splitFiles(f_in, fout_stem):
    logger = logging.getLogger(__name__)
    chrom_arr = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"]
    chrom_filehandles = []

    try:
        fin = open(f_in, "r")
    except IOError:
        logger.error("Error in vcfUtils.splitFiles: Cannot open input vcf file "+f_in)
        print >> sys.stderr, "Error in vcfUtils.splitFiles: Cannot open input vcf file "+f_in
        exit(1)
                  
    head = get_head(f_in)
    header = head.rstrip('\n').split('\t')
    chrom_idx = header.index('CHROM')
    try:
        for chrom_num in chrom_arr:
            chrom_filehandle = open(fout_stem+'_chr'+str(chrom_num)+'.txt', 'w')
            chrom_filehandle.write(head)
            chrom_filehandles.append(chrom_filehandle)

    except IOError:
        logger.error("Error in vcfUtils.splitFiles: Cannot open input vcf file "+f_in)
        print >> sys.stderr, "Error in vcfUtils.splitFiles: Improper output file specification "+f_in 
        exit(1)
    
    written_chrs = {}
    
    
    while 1:
        temp = fin.readline()
        found_chr = False
        if not temp:
            break
        elif(temp.startswith('#')):
            continue
        else:
            line_contents = temp.rstrip('\n').split('\t')
            chrom_str = str(line_contents[chrom_idx])
            for idx,chrom_num in enumerate(chrom_arr):
                if(chrom_str == str(chrom_num) or chrom_str == 'chr'+str(chrom_num)):
                    chrom_filehandle = chrom_filehandles[idx]
                    chrom_filehandle.write(temp)
                    written_chrs[chrom_num] = 1
                    found_chr = True
                    break
            if not found_chr:
                logger.warn('Warning: could not find file to put variant in chr ' + chrom_str +' into.')
    
    for chrom_num in chrom_arr:
        if chrom_num not in written_chrs:
            logger.warn('Warning: no variants found for chr ' + str(chrom_num))

def getSampleTableName(vcf_file_loc):
    return 'sample_' + getSampleName(vcf_file_loc)

#ALERT WE DONT THINK THIS FUNCTION IS WORKING PROPERLy
def getSampleName(vcf_file_loc):
    return vcf_file_loc[0:len(vcf_file_loc) - 4].replace('.', '_').replace('-', '_').split('/')[-1]

def getFilename(file_path):
    return os.path.split(file_path)[-1]

#gzips a file path.  It is important to perform this after many of the vcf processing steps, as the by default return a vcf, not a vcf.gz
def gzip_file(file_path):
    gzcmd = 'gzip {vcf}'.format(vcf = file_path)
    subprocess.Popen(gzcmd, shell=True).wait()
    return file_path + '.gz'



def stripChrPrefix(vcf_filepath, out_dir, skip_if_exists=False):
    logger = logging.getLogger(__name__)
    logger.info('Stripping chr prefix from VCF CHROM column (if present)')
    vcf_name = getSampleName(vcf_filepath)
    outfilepath = os.path.join(out_dir, vcf_name+'_strippedChr' + '.vcf')
    if(skip_if_exists and os.path.isfile(outfilepath)):
        logger.debug('Reusing existing stripped chr VCF file at ' + outfilepath)
        return outfilepath
    #else
    #strip chr prefix

    vcf = stmp_annotation_util.open_compressed_or_regular(vcf_filepath, 'r')
    out = stmp_annotation_util.open_compressed_or_regular(outfilepath, 'w')
    
    for line in vcf:
        if(line.startswith('#')):
            out.write(line)
            continue
        #else
        lineContents = line.rstrip("\n").split("\t")
        chrom = lineContents[0]
        if(chrom.startswith('chr')):
            chrom = chrom[3:] # strip 'chr' prefix
        lineContents[0] = chrom
        newline = "\t".join(lineContents)+"\n"
        out.write(newline)

    return outfilepath


# splits multiallelic site in VCF into separate rows using bcftools
def splitMultiallelic(vcf_filepath, out_dir, skip_if_exists=False):
    logger = logging.getLogger(__name__)
    logger.info('Converting VCF so there is only 1 allele per line')
    vcf_name = getSampleName(vcf_filepath)
    #THIS HAS BEEN CHANGED BC PRAG DID SPURIOUS SHIT WITH THE .GZ suffix
    outfilepath = os.path.join(out_dir, vcf_name + '_noMultiallelic' + '.vcf')
    if(skip_if_exists and os.path.isfile(outfilepath)):
        logger.debug('Reusing existing converted VCF file at ' + outfilepath)
        return outfilepath
    #else
    # bcftools norm: Left-align and normalize indels, check if REF alleles match the reference, split multiallelic sites into multiple rows; recover multiallelics from multiple rows.
    # -m, --multiallelics: split multiallelic sites into biallelic records (-) or join biallelic sites into multiallelic records (+). 
    # -f, --fasta-ref FILE: reference sequence. 
    # -o, --output FILE
    cmd = "bcftools norm -m - --fasta-ref='{fasta}' '{vcf_filepath}' -o z -o '{out_dir}'".format(vcf_filepath=vcf_filepath, out_dir=outfilepath, fasta='/share/PI/euan/apps/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa')
    logger.debug('bcftools split multiallelic command: ' + str(cmd))
    
    subprocess.Popen(cmd, shell=True).wait()
    outfilepath = gzip_file(outfilepath)
    return outfilepath

#take from the udn pipelines general utils
def rreplace(s, old, new, num_occurrences=1):
    li = s.rsplit(old, num_occurrences)
    return new.join(li)

#reheaders the vcf, so that sample id = udnid
#formerly housed in the UDN pipeline
#Currently this file has an unsolved error 

def reheaderVcf(vcf_filepath, out_dir, id='UDN390934'):
    #bcftools reheader [OPTIONS] file.vcf.gz
    #Modify header of VCF/BCF files, change sample names.
    #-h, --header FILE
    #new VCF header
    #-o, --output FILE
    #see Common Options
    #-s, --samples FILE
    #for some reason echoes the id to create a log file I think?
    reheadered_vcf_path = rreplace(vcf_filepath, '.vcf', '.reheader.vcf', num_occurrences=1)
    #output files come uncompressed so we write them to a properly uncompressed filename
    reheadered_vcf_path = reheadered_vcf_path.strip('.gz')
    cmd = 'echo "{id}"|bcftools reheader -s - -o "{reheadered_vcf_path}" "{non_reheadered_vcf_path}"'.format(id=id, reheadered_vcf_path=reheadered_vcf_path, non_reheadered_vcf_path=vcf_filepath)
    subprocess.Popen(cmd, shell=True).wait()

    outfilepath = gzip_file(reheadered_vcf_path)

    return outfilepath

def stripVCFComments(vcf_filepath, out_dir):
    logger = logging.getLogger(__name__)
    logger.info('Removing extra fields (##) from VCF')
    vcf_name = getSampleName(vcf_filepath)
    outfilepath = os.path.join(out_dir, vcf_name + '_noComments' + '.vcf')
    cmd = 'grep -v "^##" {vcf_filepath} > {vcf_outfile}'.format(vcf_filepath = vcf_filepath, vcf_outfile = outfilepath)
    logger.debug('command to strip vcf comments: ' + str(cmd))
    subprocess.Popen(cmd, shell=True).wait()
    return outfilepath


# helper function to convert a VCF consisting of per-sample allele info to a VCF with just allele freqs using bcftools
# currently waits for completion, but could be easily made async (immediately returning what the output file name would be after completion)
def vcf2allelefreq(vcf_file_loc, output_dir, overwriteIfExists=False):
    logger = logging.getLogger(__name__)
    sample_name = getSampleName(vcf_file_loc)
    outfilename = sample_name + '_alleleFreqs.vcf'
    outfilepath = os.path.join(output_dir, outfilename)
    if((not overwriteIfExists) and os.path.isfile(outfilepath)):
        logger.warn('skipped conversion of ' + sample_name + ' to allele freq VCF because already converted')
    else:
        cmd = "bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AC\t%AN\n' " + vcf_file_loc + " | awk '{OFS=\"\t\";print $0,$6/$7}' > " + outfilepath
        logger.debug('vcf to allele freq bcftools cmd: ' + str(cmd))
        subprocess.Popen(cmd, shell=True).wait()
    
    return outfilename


def getVCFHeader(vcf_file_loc):
    logger = logging.getLogger(__name__)
    vcf = stmp_annotation_util.open_compressed_or_regular(vcf_file_loc, 'r')
    for line in vcf:
        if(line[0] == '#' and line[1] != '#'):
            vcf.close()
            return line[1:].rstrip("\n")
    
    # else
    logger.error('Error no header found in vcf ' + vcf_file_loc, exc_info=True)
    vcf.close()
    return ''

#returns filename without extension
def stripExtension(filename):
    return filename.split('.')[0]

# converts VCF to TSV, extracting each specified tag from the INFO column and listing it as a seperate column in the output TSV
def vcf2tsv(vcf_path, tags, output_dir, output_extension='.tsv'):
    logger = logging.getLogger(__name__)
    vcf_filename = getFilename(vcf_path)
    vcf_filename_noExtension = stripExtension(vcf_filename)
    outpath = os.path.join(output_dir, vcf_filename_noExtension+output_extension)
    vcfHeader = getVCFHeader(vcf_path)
    vcfCols = vcfHeader.split("\t")
    header = '#' + "\t".join(vcfCols[:-1]) + "\t" + "\t".join(tags)
    #write header
    db_utils.writef(outpath, header)
    # then write data
    tagStrs = []
    for tag in tags:
        tagStrs.append('%INFO/{tag}'.format(tag=tag))
    vcfColStrs = []
    for col in vcfCols:
        vcfColStrs.append('%' + col)
    cmd = "/usr/bin/time -v bcftools query -f '{vcfCols}\t{joinedTags}\n' {vcf}".format(vcfCols = "\t".join(vcfColStrs[:-1]), joinedTags="\t".join(tagStrs), vcf=vcf_path)
    logger.debug('vcf to tsv cmd (bcftools extract tags from INFO col): ' + str(cmd))
    subprocess.Popen(cmd, stdout = stmp_annotation_util.open_compressed_or_regular(outpath, 'a'), shell=True).wait()
    return outpath

