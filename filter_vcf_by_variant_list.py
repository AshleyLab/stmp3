#written by noah friedman
#given a vcf, call bcftools to filter a vcf to only include the variants included in the list
import subprocess

def filter_vcf_by_variant_list(variantTextFile, inputVcf, outputFileName):
	#first we need to bgzip the file
	#check if the file is not compressed, if it is compress it
	vcf = inputVcf
	if '.gz' not in vcf:
		cmd = 'bgzip {iVcf}'.format(iVcf = inputVcf)
		vcf = inputVcf + '.gz' 
	cmd = 'bcftools view {iVcf} -R {variantTFile} -o {outputFName}'.format(
			iVcf = vcf,
			variantTFile = variantTextFile,
			outputFName = outputFileName
		)
	print 'bcftools command: ', cmd
	subprocess.Popen(cmd, shell=True).wait()
	print outputFileName
	return outputFileName
