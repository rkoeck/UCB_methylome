# outlier site identification using threshold-based technique
###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: threshold: lower quartile - 3x IQR or upper quartile + 3x IQR

# input: preprocessed beta values generated using 01preprocessingSWAN.R
#        conducted to compare G5 and HTF within the IVF cohort (including and excluding pregnancy complications)

# output: outlier status per site per sample (.csv)
#         summary of outliers (per type) per sample (.csv)

###########################################################################################################################

# load packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(parallel))

# load the data

annotation = fread("annotation.csv") %>%
  mutate(Sample_ID = as.character(Sample_ID))

betas = fread("betas.csv") %>%
  select(ID, all_of(annotation$Sample_ID)) %>%
  column_to_rownames("ID") %>%
  as.matrix()

# calculate the IQR, lower quartile and upper quartile at each probe
iqr = rowIQRs(betas, na.rm = T)

lowerQ = rowQuantiles(betas, probs = 0.25, na.rm = T)

upperQ = rowQuantiles(betas, probs = 0.75, na.rm = T)

# add the iqr values to the betas object

betas2 = as.data.frame(betas) %>%
  rownames_to_column("ID")

betas2$iqr = iqr

betas2$lowerQ = lowerQ

betas2$upperQ = upperQ

# make a new variable for 3x IQR and upper and lower threhsolds (25th percentile - 3x IQR and 75th percentile + 3x IQR)

betas2 = betas2 %>% mutate(tripleIQR = 3 * iqr,
                         lower = lowerQ - tripleIQR,
                         upper = upperQ + tripleIQR)

# determine whether each beta value lies within the normal range or represents a hypo/hyper methylation outlier

betas3 = betas2 %>% mutate(across(.cols = all_of(annotation$Sample_ID), .fns = ~if_else(.x < lower, "hypo", 
                                                                                        if_else(.x > upper, "hyper", "FALSE"))))
# write the raw data as a file in case this is required for further processing

write.csv(betas3, file = "rawThresholdData.csv", row.names = FALSE)

# for each sample calculate how many outliers/non-outliers/NAs there were

hypo = betas3 %>% select(all_of(annotation$Sample_ID)) %>%
  summarise_all(~sum(. == "hypo", na.rm = T))

row.names(hypo) = "hypo"

hyper = betas3 %>% select(all_of(annotation$Sample_ID)) %>%
  summarise_all(~sum(. == "hyper", na.rm = T))

row.names(hyper) = "hyper"

outlier = betas3 %>% select(all_of(annotation$Sample_ID)) %>%
  summarise_all(~sum(. %in% c("hyper", "hypo"), na.rm = T))

row.names(outlier) = "outlier"

false= betas3 %>% select(all_of(annotation$Sample_ID)) %>%
  summarise_all(~sum(. == "FALSE", na.rm = T))

row.names(false) = "false"

samples = rbind(false, hyper, hypo, outlier)

# write this summary out

write.csv(samples, "sampleThresholdOutliers.csv", row.names = T)
