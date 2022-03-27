# RAD-seq
Radseq analysis project scripts

# Stacks

1.processradtags.sh
Load Stacjs
Run Process radtags. We already ran trimscript.sh to clean reads and trim to correct length so this script should be quick, we should not loose many reads and makes sure all files are in correct format for refmap.pl Stacks

2.BWA to BAM
Load BWA and SAMtools
Index reference genome and run BWA MEM on all your samples

3.Stacks refmap.pl
Load Stacks
Using reference genome - make sure reference genome is indexed see script (BWA index)
4.Stacks Populations.sh
-Careful of filtering options to make the next step easier

# Filter
Load VCFtools
Filter VCF to set criteria i.e Missing data

Remove outliers from PCAdapt and Bayescan
# Analyses
Structure/ FastStructure and running these programs in Structure_Threader 

Fineradstructure

Dsuite

Treemix

Bayescan

PCAdapt





