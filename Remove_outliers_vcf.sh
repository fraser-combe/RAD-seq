#!/bin/bash
#SBATCH --job-name=outliers
#SBATCH --time=02:00:00
#SBATCH --mem=24G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --output=outliersformat.%j
#### SLURM 1 processor test to run for 24hours
####
## Final steps to get neutral and adaptive SNP set and correct file formats
module load VCFtools/0.1.16-GCC-9.3.0
# define arguments

inputFile=vcffile.vcf


#### adaptive SNPs ####

## step 1: create a final list of all outlier loci positions

cat pathtobayescanfile/outl_pos_bayesc.txt pathtopcadaptfile/outl_pos_pcadpt.txt > outl_pos_combi.txt

# how many are detected by both methods?
nPosOutDup=`sort outl_pos_combi.txt | uniq -d | wc -l` 
echo ""$nPosOutDup" outlier loci detected with both methods"

# extract total number of unique outlier positions (for file naming)
nPosOut=`sort outl_pos_combi.txt | uniq -u | wc -l` 
echo ""$nPosOut" unique outlier loci in total"

# remove leading whitespaces (from here: https://www.cyberciti.biz/faq/bash-remove-whitespace-from-string/)
shopt -s extglob # turn it on
nPosOut="${nPosOut##*( )}"
shopt -u extglob # turn it off

# extract file with unique positions
sort outl_pos_combi.txt | uniq -u > outl_pos_"$nPosOut".txt
#rm outl2_pos_combi.txt

## step 2: subset filtered vcf file by outlier positions

vcftools --vcf "$inputFile" --positions outl_pos_"$nPosOut".txt --recode --recode-INFO-all --out adaptive_"$nPosOut"
mv adaptive_"$nPosOut".recode.vcf adaptive_"$nPosOut".vcf

#### neutral SNPs ####

## step 1: create a final list of all neutral loci positions

# get list of all positions original vcf file
grep -v "^##" "$inputFile" | cut -f1-2 | sed '1d' > all_pos.txt

# remove previously identified outlier positions to only retain neutral ones
cat all_pos.txt outl_pos_"$nPosOut".txt | sort | uniq -u > ntrl_pos_preHWE.txt

## step 2 : subset filtered vcf file by neutral outlier positions
vcftools --vcf "$inputFile"  --exclude-positions outl_pos_combi.txt --recode --recode-INFO-all --out neutral_preHWE
#vcftools --vcf "$inputFile" --positions ntrl_pos_preHWE.txt --recode --recode-INFO-all --out neutral_preHWE

## step 3 : apply HWE filter
vcftools --vcf neutral_preHWE.recode.vcf --hwe 0.000001 --recode --recode-INFO-all --out neutral
nPosNtrl=`grep  -c -v "^#" neutral.recode.vcf`
mv neutral.recode.vcf neutral_"$nPosNtrl".vcf
mv neutral.log neutral_"$nPosNtrl".log
