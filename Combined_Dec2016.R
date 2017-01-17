#First get the libraries required for script
library(edgeR)
library(limma)
library(ggplot2)
library(statmod)
library(RColorBrewer)
library(tidyr)

###############################################
##### Get Pauls Data Round 1 ##################
###############################################

#First read in the file
sample1 = read.csv(file="./PaulsData/20150911/counts.txt", sep="\t",stringsAsFactors=FALSE,skip=1,row.names=1)
sample2 = read.csv(file="./PaulsData/20150603/counts.txt", sep="\t",stringsAsFactors=FALSE,skip=1,row.names=1)
sample3 = read.csv(file="./PaulsData/20150701/counts.txt", sep="\t",stringsAsFactors=FALSE,skip=1,row.names=1)
sample4 = read.csv(file="./PaulsData/20160115/counts.txt", sep="\t",stringsAsFactors=FALSE,skip=1,row.names=1)

#Extract the gene lengths
glengths = sample1[,"Length",drop=FALSE]

#Just keep counts
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

#Tidy column names
names(s1) = gsub("_ML.*","",names(s1))
names(s1) = gsub("PE15R*.","",names(s1))
names(s2) = gsub("_S.*","",names(s2))
names(s3) = gsub("_ML.*","",names(s3))
names(s4) = gsub("_ML.*","",names(s4))
names(s4) = gsub("mapped*.","",names(s4))

###############################################
##### Get Pauls Data for round 2 ##############
###############################################

#First read in the file
sample5 = read.csv(file="./PaulsData/20160525/counts.txt", sep="\t",stringsAsFactors=FALSE,row.names=2)
sample6 = read.csv(file="./PaulsData/20160610_20160722/counts.txt", sep="\t",stringsAsFactors=FALSE,skip=1,row.names=1)

#Extract the gene lengths
glengths = sample1[,"Length",drop=FALSE]

#Just keep counts
s5 = sample5[,7:18]
s6 = sample6[,6:42]


#Tidy column names
names(s5) = gsub("_ML.*","",names(s5))
names(s5) = gsub("_merged.Aligned.out.bam","",names(s5))
names(s5) = gsub("mapped*.","",names(s5))

names(s6) = gsub("_ML.*","",names(s6))
names(s6) = gsub("_merged.Aligned.out.bam","",names(s6))
names(s6) = gsub("mapped*.","",names(s6))


####################################################
## Pauls Data Round 3 ##############################
####################################################

#First read in the file
sample7 = read.csv(file="./PaulsData/20160922/batch1/counts.txt", sep="\t",stringsAsFactors=FALSE,skip=1,row.names=1)
sample8 = read.csv(file="./PaulsData/20160922/batch2/counts.txt", sep="\t",stringsAsFactors=FALSE,skip=1,row.names=1)
sample9 = read.csv(file="./PaulsData/20161007/counts.txt", sep="\t",stringsAsFactors=FALSE,skip=1,row.names=1)

#Extract the gene lengths
glengths = sample1[,"Length",drop=FALSE]

#Just keep counts
s7 = sample7[,6:9]
s8 = sample8[,6:9]
s9 = sample9[,6:9]


#Tidy column names
names(s7) = gsub("_ML.*","",names(s7))
names(s7) = gsub("_merged.Aligned.out.bam","",names(s7))
names(s7) = gsub("mapped*.","",names(s7))

names(s8) = gsub("_ML.*","",names(s8))
names(s8) = gsub("_merged.Aligned.out.bam","",names(s8))
names(s8) = gsub("mapped*.","",names(s8))
names(s8) = gsub(".merged","",names(s8))

names(s9) = gsub("_ML.*","",names(s9))
names(s9) = gsub("_merged.Aligned.out.bam","",names(s9))
names(s9) = gsub("mapped*.","",names(s9))
names(s9) = gsub(".merged","",names(s9))               
 
####################################################
## Pauls Data Round 4 ##############################
####################################################

sample10 = read.csv(file="./PaulsData/20161207/counts.txt", sep="\t",stringsAsFactors=FALSE,skip=1,row.names=1)
s10 = sample10[,6:15]
names(s10) = gsub("_ML.*","",names(s10))
names(s10) = gsub("_merged.Aligned.out.bam","",names(s10))
names(s10) = gsub("mapped*.","",names(s10))
names(s10) = gsub(".merged","",names(s10))                
                 
#Combine samples
paul.counts = cbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10)
names(paul.counts) = gsub("EKL\\.","EKL",names(paul.counts))
names(paul.counts) = gsub("\\.","-",names(paul.counts))
names(paul.counts) = gsub("_","-",names(paul.counts))

#Produce DGEList
y <- DGEList(counts=paul.counts, genes=rownames(paul.counts))


#Save batch information
sample_group = c(rep(1,ncol(s1)),rep(2,ncol(s2)),rep(3,ncol(s3)),
                 rep(4,ncol(s4)),rep(5,ncol(s5)),rep(6,ncol(s6)),
                 rep(7,ncol(s7)),rep(8,ncol(s8)),rep(9,ncol(s9)),rep(10,ncol(s10)))
targets = data.frame(Sample=colnames(y$counts),batch = sample_group)

#Some of the files don't have a - in the name, fix this
targets$Sample = gsub('MLM1','MLM-1',targets$Sample)


#### remove genes with very low counts (i.e only keep those rows[genes] with greater than 1 count per million in at least 5 samples)
#keep <- rowSums(cpm(y)>=1) >= 5
#y.keep <-y[keep,]
y = calcNormFactors(y) #Normalise the raw counts

#Transform to logFPKM
FPKM = rpkm(y,log=TRUE,gene.length = glengths$Length)


#Get the sample information from Pauls strange excel sheet (which is incomplete....)
stargets = read.table(file="Sample_type_only.txt", sep="\t",stringsAsFactors=FALSE,header=1)

#Keep order of samples from targets file, which is the same as the FPKM
korder <- targets$Sample
tt = merge(targets,stargets,by="Sample",all.x=TRUE)
tt = tt[match(korder,tt$Sample),]

#Replace NA
tt[is.na(tt)] <- "Unknown"

tt$Cell_Type = ifelse(grepl('T-ALL',tt$Type),'T-ALL', 
                      ifelse(grepl('other',tt$Type) | grepl('other ',tt$Type) | grepl('Other ',tt$Type),'Other',
                           ifelse(grepl('Unknown',tt$Type),'Unknown',
                                  ifelse(grepl('Ph',tt$Type),'Clinical Phlike',tt$Type))))

##########################################################
# Classify these samples with AllSorts ###################
##########################################################

library(AllSorts)
sfpkm = streamline(y$counts,glengths$Length)
classed = classify(sfpkm,c(0.25,0.25,0.75,0.5))

#Combine with targets
ttt = cbind(tt,classed)


saveRDS(FPKM,file='../ShinyDemo/Express_Yourself/Data/MCRI.Rds')
saveRDS(ttt[,c(1,3,4,6,11)],file='../ShinyDemo/Express_Yourself/Data/MCRI_target.Rds')
