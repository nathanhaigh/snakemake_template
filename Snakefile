MAX_THREADS = 32
N_BENCHMARKS = 1
singularity:
	"docker://continuumio/miniconda3:4.5.12"

############################
# Include other Snakefiles #
############################
include:
	"rules/misc.smk"

#######################################
# Convienient rules to define targets #
#######################################
localrules:
	all

rule all:
	input:
#		"reports/raw_reads_multiqc.html",
#		"reports/qc_reads_multiqc.html",
#		expand("mapped/{accession}.bam", accession=ACCESSIONS),


################
# Rules Proper #
################

rule fastqc_raw:
	input:
		"raw_reads/{prefix}.fastq.gz",
	output:
		zip  = "reports/raw_reads/{prefix}_fastqc.zip",
		html = "reports/raw_reads/{prefix}_fastqc.html",
	conda:
		"envs/default.yml"
	threads:
		MAX_THREADS
	benchmark:
		repeat("benchmarks/fastqc_raw/{prefix}.txt", N_BENCHMARKS),
	shell:
		"""
		fastqc --threads {threads} {input}
		"""
