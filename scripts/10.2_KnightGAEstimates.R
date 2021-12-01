###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: Gestational age predictor (Knight method)

# input: processed beta values generated using the script:10preprocessingKnight.R
#         applied to IVF samples only

# output: gestational age estimates (.csv)

# code source: https://github.com/akknight/PredictGestationalAge
#              also source of: suppl22.csv, suppl21.csv, suppl24.txt, cgprobesGApredictor.csv
# method replicated from: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1068-z 

###########################################################################################################################

# load packages

suppressPackageStartupMessages(library(WGCNA))
suppressPackageStartupMessages(library(impute))
suppressPackageStartupMessages(library(dynamicTreeCut))
suppressPackageStartupMessages(library(RPMM))
suppressPackageStartupMessages(library(sqldf))
suppressPackageStartupMessages(library(flashClust))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))

# load the data

#load in supplemental file 22 (file path adjusted - full location supplied)
probeAnnotation21kdatMethUsed=read.csv("suppl22.csv")

#load supplemental file 21

probeAnnotation27k=read.csv("suppl21.csv")

#Input the set of CpG coefficients here.

datClock=read.csv("cgprobesGApredictor.csv",as.is=T)

#Input your test dataset

dat0=fread("betas_GA.csv") %>%
  rename("CpGName" = "ID") %>% select(-Chromosome, -Start, -End, -Strand)

#removes probes from the annotation file that were QC'ed out of your merged dataset
#The CpG identifiers in dat0 should be named CpGName

probeAnnotation21kdatMethUsed <- probeAnnotation21kdatMethUsed %>% filter(Name %in% dat0$CpGName)

#insert supplemental file 24

source("suppl24.txt")

# get dimensions of sample data
nSamples=dim(dat0)[[2]]-1
nProbes= dim(dat0)[[1]]

#Create a log file which will be output into your directory

match1=probeAnnotation21kdatMethUsed$Name %in% dat0$CpGName

dat1= dat0 %>% filter(CpGName %in% probeAnnotation21kdatMethUsed$Name)
asnumeric1=function(x) {as.numeric(as.character(x))}
dat1[,-1]=apply(as.matrix(dat1[,-1]),2,asnumeric1)
set.seed(1)

#Normalize the Test Dataset
normalizeData=T

#load in "Normalization.R"
source("10NormalizeAndPredictGA.R")

#This Output file will contain the age prediction
write.table(datout,"Outputfile_GAEstimation.csv", row.names=F, sep="," )
dat0UsedNormalized= datMethUsedNormalized %>% column_to_rownames("samples") %>% t() %>%
  as.data.frame() %>% rownames_to_column("CpGName")
write.table(dat0UsedNormalized,file="dat0UsedNormalized_GAEstimation.csv",sep=",",row.names=F)




