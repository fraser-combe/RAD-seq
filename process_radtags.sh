#!/bin/bash
#
#SBATCH --job-name=bwa-test
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --output=o.%j.pradstacks

module load Perl
module load Python/3.7.0-iomkl-2018b
module load Stacks
#run stacks process radtag with enxyme combination to conduct quality check and make surte samples are ready for stacks pipeline.  
#If working correctly very few reads should be lost as we have already cleaned up reads using GBSpipeline from UMGC
process_radtags -p /homes/fcombe/stacks/gbs.trim  -o /homes/fcombe/stacks/output -c -q -r --renz_1 sbfI --renz_2 taqI 

