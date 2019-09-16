#!/bin/bash
#SBATCH --job-name analysis
#SBATCH --mem 1G
#SBATCH --ntasks 1 --cpus-per-task 1
#SBATCH --time 00:10:00

# Use a conda software environment for this analysis script
module load \
  miniconda3-4.6.14-gcc-5.4.0-kkzv7zk
conda activate tutorial

SAMPLES=(
  "ACBarrie"
  "Alsen"
  "Baxter"
  "Chara"
  "Drysdale"
  "Excalibur"
  "Gladius"
  "H45"
  "Kukri"
  "Pastor"
  "RAC875"
  "Volcanii"
  "Westonia"
  "Wyalkatchem"
  "Xiaoyan"
  "Yitpi"
)

# Download adapter file for trimmomatic
mkdir -p misc/trimmomatic_adapters
curl "raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-PE.fa" > misc/trimmomatic_adapters/TruSeq3-PE.fa

# Index reference genome file
#####
bwa index -p references/reference.fasta.gz -a bwtsw references/reference.fasta.gz

for SAMPLE in ${SAMPLES[@]}; do
  # FastQC the raw reads
  #####
  fastqc --threads 1 raw_reads/${SAMPLE}_R1.fastq.gz raw_reads/${SAMPLE}_R2.fastq.gz

  # Adapter/quality trim raw reads
  #####
  mkdir -p qc_reads
  trimmomatic PE \
      -threads 1 \
      raw_reads/${SAMPLE}_R1.fastq.gz raw_reads/${SAMPLE}_R2.fastq.gz \
      qc_reads/${SAMPLE}_R1.fastq.gz qc_reads/${SAMPLE}_R1.unpaired.fastq.gz \
      qc_reads/${SAMPLE}_R2.fastq.gz qc_reads/${SAMPLE}_R2.unpaired.fastq.gz \
      ILLUMINACLIP:misc/trimmomatic_adapters/TruSeq3-PE.fa:2:30:10:3:true \
      LEADING:2 \
      TRAILING:2 \
      SLIDINGWINDOW:4:15 \
      MINLEN:36

  # Map QC'd reads to the reference genome
  #####
  mkdir -p mapped
  bwa mem -t 1 \
    references/reference.fasta.gz \
    qc_reads/${SAMPLE}_R1.fastq.gz qc_reads/${SAMPLE}_R2.fastq.gz \
  | samtools view -b \
  > mapped/${SAMPLE}.bam
done
