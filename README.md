# stmp3
A third, more modest but less broken version of stmp

**usage**: python analysis_pipeline_master_script.py argument_configs.tsv 
please refer to the file template.tsv for an explanation of potential arguments

**Design philosophies**
"currentWorking[vcf, xls etc]"--> in order to facilitate the 

This script is the 2017 iteration of the pipeline sequence to medical pheontype (stmp).  It has evolved from the goals outlined in stmp2. Whereas stmp aimed to implement all aspects of an annotation pipeline, from intfrastructure to algorithms stmp3 is more modest in its ambition.  It leverages third party tools (vcfanno etc) wherever possible and aims to augment the insights and processes used by the genetic counseling team.
The main script for this repository is analysis_pipeline_master_script.py.  As the name implies, this pipeline aims to perform all aspects of analysis and processing of our sequencing data, from intial processing of bams, fastqs, and vcfs, to endpoint user motivated processing of excel sheets, powerpoint slides, and web based visualization
Here are the stages of processing data in the pipeline goes through, and the associated file suffix that is added (where relevant).  Note that file suffixes can be rendered in plain english by using the script "????.py"

**Part 0**: check coherency of arguments--stub code, *Currently not implemented*

**Part 1**: calling.  Stub for calling rtg.  *Currently not implemented*

## VCF Based Processing

**Part 2**: preprocess vcfs: for all family vcfs, run preprocessing script.  Specify arguments on the *preprocessing* line of the tsv.  Default if you only specify a proband vcf it only does preprocessing on the proband vcf.  If you have specified family vcfs, after preprocessing all of the family vcfs, it merges them together.  
**options:**
  * **smA**--split multiallelics and left normalize
  * **chP**--strip chr prefix
  * **rhP**--reheader vcf
  * **ccP**--concat snp and indel vcfs
  * **rmD**--remove duplicate records

**Part 3**: perform pre annotation filtering. Specify arguments on the *filtering* line of the configuration tsv. The goal of pre annotation filtering is to reduce the size of the vcf with filter steps before the main computationally intensive steps begin. 
**options:**
 * **sgF**--perform segregation filtering to remove all variants that do not pass segregation--currently not implemented
 * **fbL**--filter by list.  Filters by a new line separated list of chrom\tpos for variants.  The user can provide a specific list of variants to filter in the *variantListToFilterOn* line of the configuration tsv. If not, the user script generates a list of variants to filter on by reading in the chromosomes and positions from the configuration tsv *gcXls* argument

**Part 4**: annotation.  Annotates using the tool vcfanno.  Paths to the vcfs to annotate from are hardcoded in the code.  Calls the module 'prepare_vcfanno_conf.py' to prepare a conf.toml file that tells vcfanno what to do.  Note that to call vcf anno, we set our current working directory to be the directory where vcf anno is ('/home/noahfrie/noahfrie/devCode/stmp2/vcfanno/'), execute the executable vcfanno program, then reset our current working directory to be where the script itself lives (alert this may be an issue resolve??)
*Note that basically all these things can be accomplished by varsome, but in currently without calling variants in "batch" it is way too slow*
**options:**
 * **exA***--annotate from Exac ('/scratch/PI/euan/common/stmpDatafiles/ExAC.r0.3.1.sites.vep.vcf.gz') with the fields 'KG_AF_POPMAX' and 'ESP_AF_POPMAX'
 * **gnA***--annotate from gNomad ('/scratch/PI/euan/common/gnomad_data/vcf/exomes/gnomad.exomes.r2.0.1.sites.vcf.gz') with the fields 'AF_AFR', 'AF_AMR', 'AF_ASJ', 'AF_EAS', 'AF_FIN', 'AF_NFE', 'AF_OTH', 'AF_SAS', 'AN_AFR', 'AN_AMR', 'AN_ASJ', 'AN_EAS', 'AN_FIN', 'AN_NFE', 'AN_OTH', 'AN_SAS' (allele freq and allele number)
 * **clV***--annotate from clinVar ('/scratch/PI/euan/common/stmpDatafiles/clinvar_20170905.vcf.gz') with the field 'CLNSIG'
 if you would like to add more options, to check out what options you have run 'bcftools view -h vcfName', and add them where needed in the code (where I define lists like 'AF_AFR' etc)
 
 **Part 5**: perform post annotation filtering.  Ie filter out variants that have exac freq annotations that are too high.  *Currently not implemented*
 
 ## XLS based processing
 
 **Part 6**: write annotated vcf to xls. Every info field, every format field, the chromosome, the ref, and the alt becomes a column in an XLS spreadsheet. This spreadsheet will have two sheets, with the first sheet explaining the verbose explanations for every info field included in the xls.
 
 **Part 7**: merge the the ingenuity spreadsheet with the annotated spreadsheet we have generated (using the function 'merge_columns_across_spreadsheets'). We keep all the columns in the ingenuity xls and columns we specifically include from the stmp spreadsheet. Currently we use a hack to match positions from the vcf (stmp) with positions from ingenuity because ingenuity does not left align their variants.  We consider two variants to be the same if their positions are close to each other and use the position and ref alt from stmp instead of ingenuity from here on out (see the function 'find_closest_match_for_pos' for logic) 
 
 **Part 8**: annotate from web search/api based information.  
**options:**
 * **omim***--use microsoft cognitive api to search google for 'omim gene name' and return the first link found that comes from omim (TODO--maybe remove reliance on Ms cognitive api) 
 * **rvis***--annotate with rvis score by parsing the html of rvis website
 * **fathmm***--uses varsome api to get fathmm score ('fathmm_score')
 * **mTaster***--uses varsome api to get mutation taster score ('mutationtaster_converted_rankscore')
 * **sift***--uses varsome api to get sift score ('sift_converted_rankscore')
 * **phylop***--uses varsome api to get phylop score ('phylop100way_vertebrate_rankscore')
 * **mgi***--annotate with mgi information by parsing mgi html information
 
 **Part 9**: improves legibility of xls.  This involves reordering the xls and changing column names to be human readable.  Refer to the variable renameDict in the script merge_and_process_xls.py
 
 **Part 10**: export to powerpoint with the script powerpoint_export.py

### Issues, thoughts and extensions
all the annotation could be performed with varsomem, but it is too slow to call the API over and over. I would recommend eventually calling varsome with the batch calling option
