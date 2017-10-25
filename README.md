# stmp3
A third, more modest but less broken version of stmp

**usage**: python analysis_pipeline_master_script.py argument_configs.tsv 
please refer to the file template.tsv for an explanation of potential arguments

**Design philosophies**
"currentWorking[vcf, xls etc]"--> in order to facilitate the 

This script is the 2017 iteration of the pipeline sequence to medical pheontype (stmp).  It has evolved from the goals outlined in stmp2. Whereas stmp aimed to implement all aspects of an annotation pipeline, from intfrastructure to algorithms stmp3 is more modest in its ambition.  It leverages third party tools (vcfanno etc) wherever possible and aims to augment the insights and processes used by the genetic counseling team.
The main script for this repository is analysis_pipeline_master_script.py.  As the name implies, this pipeline aims to perform all aspects of analysis and processing of our sequencing data, from intial processing of bams, fastqs, and vcfs, to endpoint user motivated processing of excel sheets, powerpoint slides, and web based visualization
Here are the stages of processing data in the pipeline goes through, and the associated file suffix that is added (where relevant).  Note that file suffixes can be rendered in plain english by using the script "????.py"

**Part 0**: check coherency of arguments--stub code, doesn't do anything right now

**Part 1**: calling.  Stub for calling rtg.  Not implemented

**Part 2**: preprocess vcfs: for all family vcfs, run preprocessing script.  Specify arguments on the *preprocessing* line of the tsv.  Default if you only specify a proband vcf it only does preprocessing on the proband vcf.  If you have specified family vcfs, after preprocessing all of the family vcfs, it merges them together.  
*options:
  * **smA**--split multiallelics and left normalize
  * **chP**--strip chr prefix
  * **rhP**--reheader vcf
  * **ccP**--concat snp and indel vcfs
  * **rmD**--remove duplicate records

**Part 3**: perform pre annotation filtering. Specify arguments on the *filtering* line of the configuration tsv. The goal of pre annotation filtering is to reduce the size of the vcf with filter steps before the main computationally intensive steps begin. 
----options:
**sgF**--perform segregation filtering to remove all variants that do not pass segregation--currently not implemented
**fbL**--filter by list.  Filters by a new line separated list of chrom\tpos for variants.  The user can provide a specific list of variants to filter in the *variantListToFilterOn* line of the configuration tsv. If not, the user script generates a list of variants to filter on by reading in the chromosomes and positions from the configuration tsv *gcXls* argument

**Part 4**: annotation



all the annotation could be performed with varsomem, but it is too slow to call the API over and over. I would recommend eventually calling varsome with the batch calling option
