#!/bin/bash
#SBATCH --job-name=faststructure
#SBATCH --time=24:00:00
#SBATCH --mem=24G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --output=faststr.%j.txt

module load XZ/5.2.4-GCCcore-8.3.0
module load fastStructure/1.0-foss-2019b-Python-2.7.16
mkdir output
#run faststructure

for k in `seq 1 20`; do for i in `seq 1 10`; do fastStructure -K $k --format=str --input=faststructurecrypneutral --output=output/CPall.run_$i ; done; done


#####This output the data in output folder and runs will be recorded under CPallrun

##In order to choose the appropriate number of model components that explain structure in the dataset we use
chooseK.py --input=CPall


# To plot the best K I save as PNG but you can save as others

distruct.py -K 6 --input=CPall.run_1 --output=bestK.png

###To improve plots use R and pophelper/pophelper shiny app
library(pophelper)
library(pophelperShiny)
runPophelper()

use the shiny app to upload str files or faststructure Q plot files

Otherwise use 
https://lmme.qdio.ac.cn/StructureSelector/
