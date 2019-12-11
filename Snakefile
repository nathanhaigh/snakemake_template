# A global singularity image to be used for all jobs - need to specify --use-singularity and have singularity available on the command line
#   This image already contains the bioinformatic tools we will be using
singularity:
	"file:///shared/.singularity/workflows-training.simg"
#	"docker://continuumio/miniconda3:4.5.12"

################
# Pseudo-rules #
################
# By convention, the first rule should be called "all" and it's "input" defined as
# the list of ALL the files you want the workflow to create. e.g.:
rule all:
	input:
                expand("references/reference.fasta.gz.{ext}",
                        ext=['amb', 'ann', 'bwt', 'pac', 'sa']
                ),

################
# Rules Proper #
################
rule bwa_index:
	input:
		"{ref}",
	output:
		"{ref}.amb",
		"{ref}.ann",
		"{ref}.bwt",
		"{ref}.pac",
		"{ref}.sa",
	shell:
		"""
		bwa index \
		  -a bwtsw \
		  {input}
		"""
