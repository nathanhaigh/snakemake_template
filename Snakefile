SAMPLES = [
"ACBarrie",
"Alsen",
# "Baxter",
# "Chara",
# "Drysdale",
# "Excalibur",
# "Gladius",
# "H45",
# "Kukri",
# "Pastor",
# "RAC875",
# "Volcanii",
# "Westonia",
# "Wyalkatchem",
# "Xiaoyan",
# "Yitpi",
]

# A global singularity image to be used for all jobs - need to specify --use-singularity and have singularity available on the command line
#   This image already contains the bioinformatic tools we will be using
singularity:
	"file:///shared/.singularity/nextflow-embl-abr-webinar.simg"
#	"docker://rsuchecki/nextflow-embl-abr-webinar"

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
		expand("raw_reads/{SAMPLE}_{read}_fastqc.html",
			SAMPLE=SAMPLES,
			read=['R1', 'R2']
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

rule fastqc:
	input:
		 "raw_reads/{prefix}.fastq.gz",
	output:
		zip = "raw_reads/{prefix}_fastqc.zip",
		html = "raw_reads/{prefix}_fastqc.html",
	shell:
		"""
		fastqc \
		  --threads 1 \
		  {input}
		"""
