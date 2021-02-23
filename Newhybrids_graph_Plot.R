

 ## plots the Q-value output from New Hybrids
# make sure file is in CSV format rename headers as below then plot in ggplot2
  
  NH_output<-read.csv(file.choose(),sep=",",header=TRUE)
  ## load required packages
  require(ggplot2)
  require(reshape2)
  ## rename the columns so looks better/easier to interpret 
  colnames(NH_output) <- c("Indv", "Pure1", "Pure2","F1", "F2", "BC1", "BC2")
  NH_output$Indv <-  as.factor(NH_output$Indv) # make the individuals factors to block out the data
  
  NH_melt <- melt(data = NH_output, id.vars = "Indv") ## melt the data to allow the data to be stacked by indivudal
  colnames(NH_melt) <- c("Indv", "PopProb", "CumProb") ## rename so that its prettier
  
  ## lets give the plot some pretty colours
  col.vec <- c("red", "blue", "grey", "green", "black", "yellow", "brown")
  
  ## to be used later if decide to break the data or if there are different numbers of individuals in each population type
  # break.by <- nrow(NH_output)/6
  # break.vec <- c(break.by, (break.by*2), (break.by*3), (break.by*4), (break.by*5), (break.by*6))
  
  ## make a nice pretty little plot
  pretty.plot <- ggplot(NH_melt, aes(x = Indv, y=CumProb, fill=PopProb))
  pretty.plot+geom_bar(stat="identity", position = "stack") + scale_fill_manual(values=col.vec)+ylab("Cumulative Probability")+xlab("Individual")

###*** SET WORKING DIRECTORY HERE ***###

path.hold <- getwd()

##############################################################################

###*** SETTINGS ***###

pofz <- "/final_run/postSim_combined/NH.Results/EAON_combined_Zed.txt_Results/EAON_combined_Zed.txt_PofZ.txt"

pofz <- "C:/Users/HopeLab/Dropbox/My PC (Hope-7820)/Desktop/newhybrids/NEWHYBRIDS/output_run_1/aa-PofZ.txt"
###*** END OF SETTINGS ***###

##############################################################################


NH_plot(paste0(path.hold, pofz))

