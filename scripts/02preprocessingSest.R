# use minfi package to extract information from .idat files needed for sex prediction using SEst package

# input: data from each cohort / array type processed separately: IVF data (GSE189531), FLEHS 450k (GSE110128), 
#                                                                 ENVIRONAGE 450K (GSE151042), ENVIRONAGE EPIC (unpublished)

# output:  detection p-values of all sites & samples (.csv) 
#          raw (unnormalised/unfiltered) beta values for all sites & samples (.csv)


# load the required packages

suppressPackageStartupMessages(library(minfi))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))

# set the data directories

base.dir <- "CultMed"

data.dir <- file.path(base.dir, "idat")

sample.annotation <- "annotation.txt"

report.dir <- "cohort2/sEST"

# reaad the raw data files using minfi

targets <- fread(sample.annotation) %>%
  mutate(Basename = paste(Sentrix_ID, Sentrix_Position , sep = "_"))

RGSet <- read.metharray.exp(base = data.dir, targets = targets)

# proces the RGSet to an MSet  which has beta values

MSet <- preprocessRaw(RGSet)

# calculate beta values by converting this to a ratioSet object

RSet <- ratioConvert(MSet, what = "beta")

# extract beta values and detection p values for all sites

detPval <- detectionP(RGSet)

betas <- getBeta(RSet)

# write the resultant files as csv files for further processing

write.csv(detPval, file = file.path(report.dir,  "detPval.csv"), row.names = T)

write.csv(betas, file = file.path(report.dir, "betasRaw.csv"), row.names = T)
