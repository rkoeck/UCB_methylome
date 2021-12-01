###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: Estimate cellular composition of umbilical cord blood samples

# tissue: umbilical cord blood
# method: Houseman algorithm
# reference: FlowSorted.CordBloodCombined.450K

# input: .idat files from each cohort / array type processed separately: IVF data (GSE189531), FLEHS 450k (GSE110128), 
#                                                                 ENVIRONAGE 450K (GSE151042), ENVIRONAGE EPIC (unpublished)

# output: cell composition estimnates for the specified cell types (.csv) per sample
#         the estimates are combined with the sample annotation file for downstream processing

###########################################################################################################################

# load packages

suppressPackageStartupMessages(library(FlowSorted.Blood.EPIC))
suppressPackageStartupMessages(library(minfi))
suppressPackageStartupMessages(library(FlowSorted.CordBloodCombined.450k))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ExperimentHub))
suppressPackageStartupMessages(library(FlowSorted.Blood.EPIC))

# set up the environment for the cell type estimation

# load the list of IDOL optimised CpGs for cord blood deconvolution
data("IDOLOptimizedCpGsCordBlood")

# load the reference data (flow sorted umbilical cord blood samples)
hub <- ExperimentHub()

myfiles <- query(hub, "FlowSorted.CordBloodCombined.450k")

FlowSorted.CordBloodCombined.450k <- myfiles[[1]]

# extract from it the names of the probes (the EPIC data will need to be filtered to contain just these probes)

flow.manifest <- getManifest(FlowSorted.CordBloodCombined.450k)

probes <- getManifestInfo(flow.manifest, type = "locusNames")

# probes contains all the IDOL optimised CpGs for cell type deconvolution

message("loading sample data")

# set data input and output directories

base.dir <- "CultMed"

data.dir <- file.path(base.dir, "idat")

sample.annotation <- "annotation.txt"

report.dir <- "cellComposition"  

#the Houseman algorithm is used via the minfi package

# load in metadata and creare Basename column for desired filenames

targets <- fread(sample.annotation) %>% 
  mutate(Basename = paste(Sentrix_ID, Sentrix_Position , sep = "_"))

# read in the red and green files

RGSet <- read.metharray.exp(base = data.dir, targets = targets)

# Add a "Sample_Name"column to the RGSet to match the reference data

colData(RGSet)$Sample_Name <- colData(RGSet)$Basename

# Add sample names to the RGSet object

sampleNames(RGSet) <- colData(RGSet)$Sample_Name

# extract the sample names from reference data and culture media data

message("correcting object dimensions")

referenceSamples <- sampleNames(FlowSorted.CordBloodCombined.450k)

mediaSamples <- sampleNames(RGSet)

# combine the 2 datasets so that the probes are matched up

test <- combineArrays(FlowSorted.CordBloodCombined.450k, RGSet, outType = "IlluminaHumanMethylation450k")

# split the test object back into the reference dataset and the samples of interest

RGSet1 <- test[ , mediaSamples]

FlowSorted.CordBloodCombined.450k <- test[ , referenceSamples]

# calculate the cell type composition

# method =  IDOL 

# estimate the cell counts

message("starting IDOL")

countsIDOL <- estimateCellCounts2(RGSet1,
                                  compositeCellType = "CordBloodCombined",
                                  processMethod = "preprocessNoob",
                                  probeSelect = "IDOL",
                                  cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC"),
                                  referencePlatform = "IlluminaHumanMethylation450k",
				                          IDOLOptimizedCpGs = IDOLOptimizedCpGsCordBlood,
                                  refereceset = "FlowSorted.CordBloodCombined.450k")

# write the resultant cell compositions to a .csv file

message("writing csv")

write.csv(x = countsIDOL, file = file.path(report.dir, "cellCompositionIDOL.csv"))