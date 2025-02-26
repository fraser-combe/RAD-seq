library(vcfR)
library(hierfstat)
library(adegenet)
library(ggplot2)
library(radiator)
library(pcadapt)
setwd("C:/radiotor")
bayescandeer <- genomic_converter(
  data = "populations.haps.vcf", strata = "strata.filtered.tsv",
  output = c("bayescan","genlight"), filename="deer.geste",parallel.core = 1L)

bayescandeer <- genomic_converter(
  data = "hybrids_v3.vcf", strata = "strata.filtered.tsv",
  output = c("bayescan"), filename="deer.geste",parallel.core = 1L)
##I get the following error so use another method in Hierfstat below
#Error in radiator::write_bayescan(data = input, pop.select = pop.select,  : 
 #                                   object 'pop.select' not found

vcf <- read.vcfR("populations.haps.vcf")
vcf <- read.vcfR("populations.snps.vcf")
vcf <- read.vcfR("deerwriteone_populations.snps.vcf")
vcf <- read.vcfR("reoveoutgroupsnps.vcf") #18k snps with outgroup removed
vcf <- read.vcfR("populations2speciesmd_wt.snps.vcf")#ran populations on just 2 populations 18002k snps
vcf <- read.vcfR("hybrids_v3.vcf")
vcf <- read.vcfR("deer25018removeout.vcf")


pop_map <- read.table("strata.filtered.tsv", header=TRUE, stringsAsFactors = TRUE)
#If not exculding ind run this
genind <- vcfR2genind(vcf)
genind@pop <- pop_map$STRATA
hierfstat <- genind2hierfstat(genind)

#run for exculded individuals
hierfstat <- genind2hierfstat(data2)
####Seperating into just WTD and MD and exclude outgroup
pop_map <- read.table("strata.filtered2popsnoout.tsv", header=TRUE, stringsAsFactors = TRUE)

pop_map <- read.table("strata.filtered3popsnoout.tsv", header=TRUE, stringsAsFactors = TRUE)
removeInd <- c("outg001_R1", "outg002_R1", "outg003_R1")
# remove individuals from *genind object*
# note that in this step there's no longer a comma needed before the closing square bracket
#check ind names
(indNames(data2))
data2 <- genind[!row.names(genind@tab) %in% removeInd]
data2@pop <- pop_map$STRATA

hierfstat <- genind2hierfstat(data2)

#changed to heirfstat2 for the renamed file
write.bayescan(hierfstat,fn="bayescan17481hybswappedpops.bsc",diploid=TRUE)

#Run Bayescan see text file for instructions


### Identifying outliers from Bayescab results
bayescan=read.table("bayescanRunOD10_fst.txt") 
bayescan=read.table("bayescanLongRunOD100_fst.txt") 
bayescan=read.table("outputtest01_fst.txt") 

bayescan=read.table("C:/radiotor/bayescan2popOD10_fst.txt")
bayescan=read.table("outputtest3_fst.txt") 
bayescan=read.table("bayesc_2pop_run2_18k_fst.txt") 

bayescan=read.table("C:/radiotor/bayesc_2pop_run1_18k2_fst.txt") 
bayesc_2pop_run1_18k_fst.txt
#The first column of the bayescan dataframe is the SNP ID. The next three columns (prob, log10(P0), and qval)
#are related to the test of local adaptation considering the logarithm of the posterior odds - log10(PO) - and
#the q-value for the model with selection. The fifth column gives the size of the locus-specific effect (alpha parameter).
#The last one provides the locus-specific FST averaged over all populations.

##run this in beocat on populations.haplotypes.vcf file

grep -v "#" populations.haps.vcf | cut -f 3 > idsnps2.txt



SNPb=read.table("idsnps.txt",header=FALSE)
SNPb=read.table("idsnps17481.txt",header=FALSE) #hybrid_v3.vcf
SNPb=read.table("idsnps17293.txt",header=FALSE)
SNPb=read.table("18ksnps.txt",header=FALSE)

SNPb=read.table("idsnps2pops.txt",header=FALSE)
#Merge the names of the outliers with the results from the bayescan dataframe.

bayescan=cbind(SNPb, bayescan) 
#Rename columns
colnames(bayescan)=c("SNP","PROB","LOG_PO","Q_VALUE","ALPHA","FST") 

write.table(bayescan, "deer-bayescan-results.txt", quote=FALSE, sep="\t", row.names=FALSE) 

#Change the value of the Q_VALUE column: 0 == 0.0001.
attach(bayescan)
class(bayescan$Q_VALUE)

bayescan$Q_VALUE <- as.numeric(bayescan$Q_VALUE) 
bayescan[bayescan$Q_VALUE<=0.001,"Q_VALUE"]=0.001 

#Round the values.

bayescan$LOG_PO <- (round(bayescan$LOG_PO, 4)) 
bayescan$Q_VALUE <- (round(bayescan$Q_VALUE, 4)) 
bayescan$ALPHA <- (round(bayescan$ALPHA, 4)) 
bayescan$FST <- (round(bayescan$FST, 6))

#Add a column for the type of selection grouping based on a Q-VALUE < 0.05. You can also choose a Q-VALUE < 0.01 
#if you want to be more conservative.

bayescan$SELECTION <- ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE<=0.05,"diversifying",ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE>0.05,"neutral","balancing")) 
bayescan$SELECTION<- factor(bayescan$SELECTION)
levels(bayescan$SELECTION) 

#hotho#Save the results of the SNPs potentially under positive (divergent) and balancing selection (qvalue < 0.05).

positive <- bayescan[bayescan$SELECTION=="diversifying",] 
neutral <- bayescan[bayescan$SELECTION=="neutral",] 
balancing <- bayescan[bayescan$SELECTION=="balancing",]

#Check the number of SNPs belonging to each category.

xtabs(data=bayescan, ~SELECTION) 

#Write the results of the SNPs potentially under selection (qvalue < 0.05).

write.table(neutral, "neutral.txt", row.names=F, quote=F)
write.table(balancing, "balancing.txt", row.names=F, quote=F) 
write.table(positive, "positive.txt", row.names=F, quote=F) 

#Transformation Log of the Q value in order to create the ggplot graph.
range(bayescan$Q_VALUE) 
bayescan$LOG10_Q <- -log10(bayescan$Q_VALUE) 

# Use ggplot to create a nice graph
#Create a title for the ggplot graph.

x_title="Log(q-value)" 
y_title="Fst" 
ggplot(bayescan,aes(x=LOG10_Q,y=FST)) +
  geom_point(aes(fill=SELECTION), pch=21, size=2)+ 
  scale_fill_manual(name="Selection",values=c("white","red","orange"))+ 
  labs(x=x_title)+ 
  labs(y=y_title)+   
  theme_classic()

#save pdf
ggsave("bayescan_deer.pdf", dpi=600, width=5, height=5) 
dev.off()

##Or use Rplot function from bayescan



results<-plot_bayescan("C:/radiotor/outputtest3_fst.txt",FDR=0.05)

results<-plot_bayescan("C:/radiotor/bayesc_5pop_run1_25k_fst.txt",FDR=0.05)
results<-plot_bayescan("C:/radiotor/bayesc_2pop_run1_25k_fst.txt",FDR=0.05)
#
results<-plot_bayescan("C:/radiotor/bayesc_2pop_run2_25k_fst.txt",FDR=0.05)

####Run with 18k snps file removed outgroup and renamed to 2 species
results<-plot_bayescan("C:/radiotor/bayesc_2pop_run1_18k_fst.txt",FDR=0.01) ##8651 selection
results<-plot_bayescan("C:/radiotor/bayesc_2pop_run1_18k2_fst.txt",FDR=0.05)
results<-plot_bayescan("C:/radiotor/bayesc_2pop_run2_18k_fst.txt",FDR=0.05)#run 2 8651
results<-plot_bayescan("C:/radiotor/bayesc_2pop_run3_18k_fst.txt",FDR=0.01)#run 3 didnt work


results<-plot_bayescan("C:/radiotor/bayesc_2pop2_run2pr100_18k_fst.txt",FDR=0.05)
results<-plot_bayescan("C:/radiotor/bayesc_2pop_run1_25k_fst.txt",FDR=0.01)


results<-plot_bayescan("C:/radiotor/bayesc_2pop_run3pr1000_25k_fst.txt",FDR=0.01)
results<-plot_bayescan("C:/radiotor/bayesc_2pop_run1_25k_fst.txt",FDR=0.01)

## 
results<-plot_bayescan("C:/Users/frasc/Dropbox/Genomics/Bayescantommy/all13969_bayesc_run1_pr_10_fst.txt",FDR=0.01)


results$PO
results$FDR
results$FNDR
results$p
results$outliers
results$nb_outliers

# detect outliers
outliers_bayescan <-results$outliers # position of outliers
results$nb_outliers #68 outliers for deer 27 


path_to_file<-"C:/Users/frasc/Dropbox/Hope Lab/NM_Shrew_Radseq/referencegenome/VCFfiles/filtering/cryp_075ind_p3_maxmiss25.recode.vcf"

# extract corresponding SNP names and positions
vcfPath<-path_to_file
vcfPath="C:/radiotor/populations2.haps.vcf"

vcf <- read.vcfR(vcfPath, verbose=F)
loci <- as.data.frame(vcf@fix[,1:2])
outlier_loci <- loci[results$nlkjljk,]
outlier_loci <- loci[results$outliers,]

# output positions table and save in your wd
write.table(outlier_loci, file = paste0("outl_pos_bayesc26cry.txt"), sep = "\t", quote = F, row.names = F, col.names = F)
print(paste0("outlier positions table exported to outl_pos_bayesc.txt"), quote = 0)
getwd()
library(boa) ##############################WORKING OUtlier loci in Red
dat<-read.table("C:/radiotor/outputtest3_fst.txt",header = TRUE)
head(dat)
table(dat[,"qval"]<0.1)
outs=which(dat[,"qval"]<0.1)
plot_bayescan("outputtest3_fst.txt",FDR=0.1,add_text=FALSE,size=0.5,highlight=outs)

library(boa)
dat<-read.table("C:/radiotor/bayescanLongRunoD100_fst.txt",header = TRUE)
head(dat)
table(dat[,"qval"]<0.01)
outs=which(dat[,"qval"]<0.1)
plot_bayescan("bayescanLongRunoD100_fst.txt",FDR=0.05,add_text=FALSE,size=0.5,highlight=outs)

library(boa)
dat<-read.table("C:/radiotor/bayesc_2pop_run1_25k_fst.txt",header = TRUE)
head(dat)
table(dat[,"qval"]<0.01)
outs=which(dat[,"qval"]<0.1)
plot_bayescan("bayesc_2pop_run1_25k_fst.txt",FDR=0.05,add_text=FALSE,size=0.5,highlight=outs)



mydata=read.table("bayescanLongRunoD100.sel",colClasses="numeric")
head(mydata)
quartz()
plot(density(mydata[,parameter]), xlab=parameter, main=paste(parameter,"posterior distribution"),xlim=c(0,0.25),ylim=c(0,500))
boa.hpd(mydata[,parameter],0.05)
lines(density(mydata[,"Fst2"]))
lines(density(mydata[,"Fst3"]))
lines(density(mydata[,"Fst4"]))
boa.hpd(mydata[,parameter],0.05)

#VENN DIAGRAM
library(VennDiagram)

#Save bayescan outlier SNP names in a vector.

bayescan_outliers <- positive$SNP
bayescan_outliers <- balancing$SNP
balancing
#Check the number of outliers found by BayeScan.
length(bayescan_outliers)
length(bayescan)
## [1] 8
#Save pcadapt names in a vector.

pcadapt_outliers <- top_1percent$LOCUS

length(pcadapt_outliers)

#Create a Venn diagram.

venn.diagram(
  x = list(bayescan_outliers,pcadapt_outliers),
  category.names = c("Bayescan" , "Pcadapt"),
  filename = 'Venn_diagramm_outliers.png',
  output=TRUE,
  imagetype="png" ,
  height = 400 , 
  width = 400, 
  resolution = 300,
  compression = "lzw",
  cat.cex = 0.6,
  cat.pos = c(-5, 5))

plot_bayescan<-function(res,FDR=0.05,size=1,pos=0.35,highlight=NULL,name_highlighted=F,add_text=T)
{
  if (is.character(res))
    res=read.table(res)
  
  colfstat=5
  colq=colfstat-2
  
  highlight_rows=which(is.element(as.numeric(row.names(res)),highlight))
  non_highlight_rows=setdiff(1:nrow(res),highlight_rows)
  
  outliers=as.integer(row.names(res[res[,colq]<=FDR,]))
  
  ok_outliers=TRUE
  if (sum(res[,colq]<=FDR)==0)
    ok_outliers=FALSE;
  
  res[res[,colq]<=0.0001,colq]=0.0001
  
  # plot
  plot(log10(res[,colq]),res[,colfstat],xlim=rev(range(log10(res[,colq]))),xlab="log10(q value)",ylab=names(res[colfstat]),type="n")
  points(log10(res[non_highlight_rows,colq]),res[non_highlight_rows,colfstat],pch=19,cex=size)
  
  if (name_highlighted) {
    if (length(highlight_rows)>0) {
      text(log10(res[highlight_rows,colq]),res[highlight_rows,colfstat],row.names(res[highlight_rows,]),col="red",cex=size*1.2,font=2)
    }
  }
  else {
    points(log10(res[highlight_rows,colq]),res[highlight_rows,colfstat],col="red",pch=19,cex=size)
    # add names of loci over p and vertical line
    if (ok_outliers & add_text) {
      text(log10(res[res[,colq]<=FDR,][,colq])+pos*(round(runif(nrow(res[res[,colq]<=FDR,]),1,2))*2-3),res[res[,colq]<=FDR,][,colfstat],row.names(res[res[,colq]<=FDR,]),cex=size)
    }
  }
  lines(c(log10(FDR),log10(FDR)),c(-1,1),lwd=1, lty=3)
  
  return(list("outliers"=outliers,"nb_outliers"=length(outliers)))
}



##PCAadapt from Seaconnect script used it for bayes scan
#----- PCAdapt -----#
# vignette: https://cran.r-project.org/web/packages/pcadapt/vignettes/pcadapt.html

# load libraries
library(pcadapt)
library(vcfR)
vcf<- "C:/radiotor/populations.snps.vcf"
vcf<-read.vcfR("C:/radiotor/hybrids_v3.vcf")
# define arguments
path_to_file <- "C:/radiotor/populations2.haps.vcf"
path_to_file <- "C:/radiotor/populations.snps.vcf"
path_to_file <- "C:/radiotor/deerwriteone_populations.snps.vcf"
path_to_file <- "C:/radiotor/populations2speciesmd_wt.snps.vcf"
path_to_file <- "C:/radiotor/hybrids_v3.vcf"


#cryptotis pcadapt on 50 and 75% missing data
path_to_file <- ("C:/Users/frasc/Dropbox/Hope Lab/NM_Shrew_Radseq/referencegenome/VCFfiles/filtering/50perc_maxmiss.recode.vcf")
path_to_file <- ("C:/Users/frasc/Dropbox/Hope Lab/NM_Shrew_Radseq/referencegenome/VCFfiles/filtering/75perc_maxmiss.recode.vcf")



# import vcf file to pcadapt format
dat <- read.pcadapt(path_to_file, type = "vcf")

# choose number K of principal components
x <- pcadapt(input = dat, K = 20)
# set plot output file name
pdf(file=paste0("RPlots.pdf"))
# scree plot of principal components
plot(x, option = "screeplot", K = 10)
# score plot of principal components
plot(x, option = "scores")
plot(x, option = "scores", i = 3, j = 4) # check PC axis 3 and 4
dev.off()
# compute outliers detection with K = 4
x <- pcadapt(dat, K = 3)
# Manhattan plot
#jpeg(filename = paste0(args[2],"_manhattan_pcadpt",".jpeg"))
plot(x, option = "manhattan")
#dev.off()
# QQ plot
plot(x , option = "qqplot", threshold = 0.05)
# histogram of p-values
#jpeg(filename = paste0(args[2],"_histP_pcadpt",".jpeg"))
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
#dev.off()
# histogram of test statistic
plot(x, option = "stat.distribution")

# choose cutoff for outlier detection
# Bonferroni correction
padj <- p.adjust(x$pvalues, method="bonferroni")
alpha <- 0.01
outliers <- which(padj < alpha)
print(paste(length(outliers), "outliers detected"),quote=0)

# extract corresponding SNP names and positions
vcf <- read.vcfR(path_to_file, verbose=F)
loci <- as.data.frame(vcf@fix[,1:2])
outlier_loci <- loci[outliers,]
getwd()
# output positions table
write.table(outlier_loci, file = paste0("247175perperc_pos_pcadpt.txt"), sep = "\t", quote = F, row.names = F, col.names = F)
print(paste0("outlier positions table exported to outl_pos_pcadpt_"), quote = 0)



plot(x, option ="screeplot")

require(qvalue)
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.1
outliers <- which(qval<alpha)
print(paste(length(outliers), "outliers detected"),quote=0)

remove outliers run again
dat <- read.pcadapt(dat[, -outliers], type = "lfmm")
xr<-pcadapt(dat,K=2)

plot(xr, option = "scores")
plot(x, option = "scores") #compare to original
dev.off()


vcf2 <- "deer25018removeout.vcf"
`require(pcadapt)
data <- read.pcadapt("deer25018removeout.vcf",type="vcf",type.out="matrix")

#Because individuals are in lines and SNPs in columns use type="lfmm"
genotypes <- read.pcadapt(data,type="lfmm")
run.pcadapt()
res <- pcadapt(genotypes,K=20)
plot(res, option ="screeplot")

###Compare pop. structure with and without outliers

poplist.names <- c(rep("Cha03",49), rep("Cha15",42))

#Without removal of outliers
x <- pcadapt(data,K=2)

#With removal of outliers
require(qvalue)
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.1
outliers <- which(qval<alpha)
x1<-read.pcadapt(data[,-outliers],type="lfmm")
xr<-pcadapt(x1,K=2)

plot(xr, option = "scores")
plot(x, option = "scores") #compare to original
