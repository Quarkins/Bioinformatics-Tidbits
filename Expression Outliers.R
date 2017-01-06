#First get the libraries required for script
library(edgeR)
library(limma)
library(ggplot2)
library(statmod)
library(RColorBrewer)
library(tidyr)

###############################################
##### Get Pauls Data ##########################
###############################################

#First read in the samples files (these are my local files, so you will need to download your own)
sample1 = read.csv(file="./PaulsData/20150911/counts.txt", sep="\t",stringsAsFactors=FALSE,skip=1,row.names=1)
sample2 = read.csv(file="./PaulsData/20150603/counts.txt", sep="\t",stringsAsFactors=FALSE,skip=1,row.names=1)
sample3 = read.csv(file="./PaulsData/20150701/counts.txt", sep="\t",stringsAsFactors=FALSE,skip=1,row.names=1)
sample4 = read.csv(file="./PaulsData/20160115/counts.txt", sep="\t",stringsAsFactors=FALSE,skip=1,row.names=1)

#Extract the gene lengths
glengths = sample1[,"Length",drop=FALSE]

#Just keep the columns with the counts
s1 = sample1[,6:13]
s2 = sample2[,6:17]
s3 = sample3[,6:12]
s4 = sample4[,6:29]

#From sample 3 remove columns MLM11 and MLM14 which are repeats of two other samples in s2
s3 = s3[,c(1,3,4,6,7)]
#Divide every element by 2 to account for the fact that we are using double counted paired end reads (see README in PaulsDara/20150911 for why)
s1 = s1/2
s2 = s2/2
s3 = s3/2
s4 = s4/2

#Combine samples into one big dataset
paul.counts = cbind(s1,s2,s3,s4)

#Tidy up sample names
colnames(paul.counts) = gsub("_.*","",colnames(paul.counts))

#Produce DGEList
y <- DGEList(counts=paul.counts, genes=rownames(paul.counts))

#### remove genes with very low counts (i.e only keep those rows[genes] with greater than 1 count per million in at least 5 samples)
keep <- rowSums(cpm(y)>=1) >= 5
y.keep <-y[keep,]


dat <- y.keep
dat = calcNormFactors(dat) #Normalise the raw counts
dat <- as.matrix(dat) #Transform to a matrix

#Log it
datlog <- cpm(dat,log=TRUE)

#Make a little function to plot expression outliers for a give gene
expressOut <- function(gene,dataset){
  y <- dataset[gene,] #Select the counts for a given gene
  x <- colnames(dataset) #Get the sample names
  med <- median(y)
  iqr <- IQR(y)
  
  #Get outliers
  outs <- (y > med+1.5*iqr | y < med-1.5*iqr)

  #Plot
  plot(y,xlab="Samples",ylab="Counts",pch=19)
  abline(h=med+1.5*iqr,col=2,lty=2)
  abline(h=med-1.5*iqr,col=2,lty=2)
  abline(h=med,col="gray",lty=1)
  text(which(outs),y[outs],labels=names(y[outs]),cex=0.7,col=2)
  legend("topleft",legend=c("Median","1.5 IQR"),
         col=c("gray",2),lty=c(1,2),cex=0.5)
  
} #Note: we dont even want to use the above fucntion

#Making fancy box plots for all the genes overview

#Read in list of genes which are of interest (this is a csv file i made based on Paul and Ians suggestions)
gene_list <- read.csv(file="../kinase.csv",sep=",")

#Subset on genes
dd <- data.frame(datlog)
dd$GeneName = row.names(dd)
m <-  dd$GeneName %in% gene_list$Gene_name
kin <- dd[m,]

#Now use tidyr to make dataset in the from required by ggplots
tt <- kin %>% gather(Sample,Counts,1:49)


#Check if a point is an outlier
check_outlier <- function(v,coef=1.5){
  iqr <- IQR(v)
  med <- median(v)
  res <- v < (med-coef*iqr) | v > (med+coef*iqr)
  return(res)
}

#check if over expressed
out_up <- function(v,coef=1.5){
  iqr <- IQR(v)
  med <- median(v)
  res <- v > (med+coef*iqr)
  return(res)
}  

#Function to check if gene under expressed  
out_down <- function(v,coef=1.5){
  iqr <- IQR(v)
  med <- median(v)
  res <- v < (med-coef*iqr)
  return(res)
}  

#Apply outlier function to dataset to find outliers
outs <- apply(kin[,1:49],1,check_outlier)
outs <- as.data.frame(outs)

#Find those which are outliers (both over and under expressed)
out_up = as.data.frame(apply(kin[,1:49],1,out_up))
out_down = as.data.frame(apply(kin[,1:49],1,out_down))


#Make a dataset for each sample list which genes it is an outlier for
outlie = data.frame(Over_Outlier=rep(0,ncol(outs)),Under_Outlier=rep(0,ncol(outs)))
for (i in 1:ncol(outs)){ 
  print(i)
  outlie[i,]= c(paste(as.character(row.names(outs)[out_up[,i]]),sep="",collapse=","),
                paste(as.character(row.names(outs)[out_down[,i]]),sep="",collapse=","))
  }
row.names(outlie) <- colnames(outs)
write.csv(file="Outliers.csv",outlie)


#Make Boxplot for each gene
pdf("Boxplot.pdf")
p <- ggplot(tt,aes(factor(GeneName),Counts,fill=factor(GeneName))) #So here we colour and factor by gene (so one boxplot per gene)
q <- p + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +xlab("Gene") + ylab("LogCPM") #Make sure to label the plot
dev.off()

###############################
# To make an interactive plot we use plotly 
###############################

#plotly 
library(plotly)
ggplotly(q)

#Make jitter plot for every gene
library(plotly)
jit <- ggplot(tt,aes(factor(GeneName),Counts,fill=factor(GeneName)))
jit = jit + geom_boxplot() + geom_point(size=1,aes(color=factor(GeneName),alpha=0.5,text=Sample)) 
ggplotly(jit)
