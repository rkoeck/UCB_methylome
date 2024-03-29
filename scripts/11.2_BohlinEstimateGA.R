###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: application of the Bohlin GA estimator

# input: beta values generated by 11preprocessingBohlin.R
#         epigenetic gestational age was only calculated for samples without pregnancy complications

# output: epigenetic gestational age estimates per sample added to the annotation sheet (.csv)

# method replicated from: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1063-4

###########################################################################################################################

# load packages

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(predictGA))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(ggplot2))

# set the data directories

base.dir <- "output/"

annotation.file <- file.path(base.dir, "annotation.csv")

betas.file <- file.path(base.dir, "betasBMIQ.csv")

  
# load the annotation file and add the additional necessary columns

annotation <- fread(annotation.file) %>% as.data.frame() %>%
  mutate(Sample_ID = as.character(Sample_ID)) %>%
  dplyr::select(-V1)

annotation$samplePlate <- if_else(annotation$Sample_Plate == "WG5839007-BCD", 1, 2)

annotation$complication <- if_else(annotation$pregnancy_complic == "no", "no", "yes")

# load the betas file - this has been BMIQ (within sample) normalised using RnBeads (as suggested by the Bohlin paper)
# only select the samples that are included in the annotation file

betas <- fread(betas.file) %>% select(-Chromosome, -Start, -End, - Strand) %>%
  column_to_rownames("ID") %>% select(all_of(annotation$Sample_ID))

# Because of the observed plate effects Bohlin recommend to carry out between sample normalisation as well

# method: ComBat
# package: sva

betas <- ComBat(betas, batch = annotation$samplePlate)

sites <- extractSites(type = "se")

sites.na <- sites[!sites %in% rownames(betas)]

# complete all missing sites with 0s - effectively removes them from the model

for(site in sites.na){
  
  betas[site ,] <- 0
  
}

mypred <- predictGA(betas)

# predictGA gives you a GA estimated in days, create also a variable in weeks
mypred$predictedGA <- mypred$GA / 7

# merge the GA predictions with the annotation file

annotation <- mypred %>% rownames_to_column("Sample_ID") %>%
  full_join(annotation, ., by = "Sample_ID")


# plot the predicted vs. actual gestational age

annotation %>% ggplot(aes(x = gestational_age, y = predictedGA)) +
  geom_point(aes(colour = MEDIUM)) +
  theme_classic() +
  geom_smooth(method = "lm", colour = "black", se = FALSE, fullrange = T) +
  scale_x_continuous(limits = c(36, 43)) +
  scale_y_continuous(limits = c(36,43)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey")

cor(annotation$gestational_age, annotation$predictedGA, method = "pearson")

# save the output

write.csv(x = annotation, file = file.path(base.dir, "estimatedGABohlin.csv"))
