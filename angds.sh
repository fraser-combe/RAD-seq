#Install angsd and HTSlib (this is already on beaocat)

#!/bin/bash
#SBATCH --job-name=angsd 
#SBATCH --time=4:00:00 
#SBATCH --mem=24G 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1
#SBATCH --output=outputangsd.%j.txt

module load Perl 
module load Python/3.7.0-iomkl-2018b 
module load GCCcore
module load HTSlib
module load SAMtools
#filter data, output .idx file bam: txt identifying bam files 
#doSaf 1: Calculate the Site allele frequency likelihood based on individual genotype likelihoods assuming HWE 
#anc: refernce genome as the ancestral state
# GL: estimate genotype likelhoods from the mapped reads (1=SAMtools model) 
#out: name of output files 
#minMapQ: Minimum mapQ quality
# minInd: Discard  the sites where we don't have data from -minInd individuals 
#nthreads: number of threads to use

/homes/fcombe/angsd/angsd -bam ~./bamfiles/bamfile.txt(list of bamfiles with location) -doSaf 1 -anc ~/ref_genome/GCA_021461705.1.fna -GL 1 -out outfiles -minMapQ 25 -minInd 5 -nThreads 8


#If you get error Please reindex fasta file
#run
#module load SAMtools
#samtools faidx GCA_021461705.1.fna
