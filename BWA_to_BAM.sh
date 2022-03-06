#!/bin/bash
#SBATCH --job-name=bwa-test
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --output=obwa.%j

#### SLURM 1 processor BWA test to run for 24 hours should be shorter

# Load the BWA module:
module load BWA/0.7.17-foss-2018b
module load SAMtools/0.1.20-foss-2018b
####
for fastq in /homes/fcombe/stacks/deer_raw/rawdata/pradoutput/*.fq.gz
do
bwa mem -t 12 /homes/fcombe/stacks/fastgbs/refgenome/refgenome.fna.gz $fastq > ${fastq%.*}.sam
################### create bam file from sam
samtools view -S -b -@ 16 ${fastq%.*}.sam > ${fastq%.*}.bam
done




###############################
In beocat then load SAMtools

then
module load SAMtools/0.1.20-foss-2018b
####bamtosortedbam
for file in *.bam; do samtools sort $file > ${file/.bam/.sorted.bam}; done

##index
for file in *.fq.sorted.bam; do samtools index $file; done 

##renameto.bamfiles
for file in *.fq.sorted.bam; do mv -- "$file" "${file%.fq.sorted.bam}.bam"; done
