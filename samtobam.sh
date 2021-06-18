#!/bin/bash
#SBATCH --job-name=bwa-test
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --output=out.%j.samtobam
module load SAMtools/0.1.20-foss-2019b
for file in *.sam; do samtools view -bS $file > ${file/.sam/.bam}; done

for file in *.bam; do sample_name=`echo $i | awk -F "." '{print $1}'`; samtools sort -@ 7 $i >  ${sample_name}.sorted.bam ; done

for file in *.bam; do samtools sort $file > ${file/.bam/sorted.bam}; done
