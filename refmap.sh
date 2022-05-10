#!/bin/bash
#SBATCH --job-name=stacksref
#SBATCH --time=100:00:00
#SBATCH --mem=96G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --output=outputstacks.%j
mkdir ref_map_new_samplenames
#### I set runtime long but will complete quicker
module load Perl
module load Python/3.7.0-iomkl-2018b
#run stacks ref.map.pl script
#Load Stacks module and run on aligned bam files
#the job of assembling the loci is being taken over by BWA and passed to refmappipeline

module load Stacks
ref_map.pl --samples Routetobamfiles/ --popmap Routetopopmap/popmap.txt -T 24 -o ./refmap

#Post-processing with populations pipeline script
