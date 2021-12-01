###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: applcation of a mixed effects linear model using the variancePartition package

# input: preprocessed beta values generated using 01preprocessingSWAN.R
#        conducted to compare G5 and HTF within the IVF cohort (including and excluding pregnancy complications)
#       and once to compare all IVF to all naturally conceived individuals

# targeted analysis: for the targeted analysis, sites were restricted to those contained in the following documents:
#                     imprinted sites: manifestImprinted.csv (source: Ginjala    V. Gene imprinting gateway. Genome Biology 2001.)
#                     birthweight associated sites: manifestBirthweight.csv (source: DOI: 10.1038/s41467-019-09671-3)

# output: (multiple testing corrected) statistics per CpG site (.csv)
#         model design matrix (.csv)

###########################################################################################################################

# load packages

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(variancePartition))

# set up the parallel environment

param = SnowParam(workers = 30, type = "SOCK", progressbar = TRUE)
register(param)

# set data directories

message("loading files")

betas.file <- "betas.csv"

annotation.file <- "annotation.csv"

# load the data files

annotation <- read.csv(annotation.file)

betas <- fread(betas.file) %>% 
	select(ID, as.character(annotation$Sample_ID)) %>% 
	column_to_rownames("ID")

# adjust the sample plate variable in the annotation file

annotation$samplePlate <- if_else(annotation$Sample_Plate == "WG5839007-BCD", "a", "b")

annotation$complication <- if_else(annotation$pregnancy_complic == "no", "no", "yes")

# for consistency of naming rownames of annotation file should correspond to colnames of betas file

annotation <- annotation %>% column_to_rownames("Sample_ID")

# create the complex model using the dream package (variancePartition)

# specify the model

 model <- ~ MEDIUM + counts.CD8T + counts.CD4T + counts.NK + counts.Bcell + 
   counts.Mono + counts.Gran + counts.nRBC + gender + gestational_age + 
   maternal_age + complication + (1 | samplePlate) + (1 | center)

# calculate the weights (this is only possible on sites containing no NA values)

weights <- voomWithDreamWeights(na.omit(betas), model, annotation)

fit <- dream(weights, model, annotation)

design <- fit$design

table <- topTable(fit, coef = "MEDIUMvg5", number = Inf)

# save the resultant objects as .csv files

write.csv(x = design, file = "dmpsDesignDream.csv")

write.csv(x = table, file = "dmpsFitDream.csv")