#!/bin/bash
#SBATCH --job-name=bwa-index
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --output=obwa.%j
################### _index_genome
module load BWA
bwa index -a bwtsw sorex_cinereus.fa
#END
