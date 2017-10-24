# stmp3
A third, more modest but less broken version of stmp

This script is the 2017 iteration of the pipeline sequence to medical pheontype (stmp).  It has evolved from the goals outlined in stmp2. Whereas stmp aimed to implement all aspects of an annotation pipeline, from intfrastructure to algorithms stmp3 is more modest in its ambition.  It leverages third party tools (vcfanno etc) wherever possible and aims to augment the insights and processes used by the genetic counseling team.
The main script for this repository is analysis_pipeline_master_script.py.  As the name implies, this pipeline aims to perform all aspects of analysis and processing of our sequencing data, from intial processing of bams, fastqs, and vcfs, to endpoint user motivated processing of excel sheets, powerpoint slides, and web based visualization
Here are the stages of processing data in the pipeline goes through, and the associated file suffix that is added (where relevant).  Note that file suffixes can be rendered in plain english by using the script "????.py"

Part 0: check coherency of arguments--stub code, doesn't do anything right now
Part 1: calling.  Stub for calling rtg.  Not implemented
Part 2: preprocess vcfs: for all family vcfs, run preprocessing script.  Default if you only specify a proband vcf it only does preprocessing on the 
