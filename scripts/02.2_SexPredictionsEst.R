###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: sex prediction the sEst package & visualisation of the results

# input: detection p-values (.csv) and raw beta values (.csv) generated in the script 02preprocessingSest.R
#         each cohort is run separately as described in 02preprocessingSest.R

# output: updated annotation file containing pedicted sex columns (.csv)
#         scatter plot of the results

###########################################################################################################################

# load packages

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(sest))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(ggplot2))

# load the data

betas = fread("rawBetaValues.csv")

pvals = fread("detPval.csv")

ann = fread("annotation.csv")

# correctly format the betas and pvals for the sest tool

betas = betas %>% column_to_rownames("V1") %>% as.matrix()

pvals = pvals %>% column_to_rownames("V1") %>% as.matrix()

# conduct sex estimation

sex = estimateSex(beta.value = betas, detecP = pvals)

# plot the results

prediction = sex$test %>% as.data.frame() %>% rownames_to_column("Basename")

ann = full_join(ann, prediction, by = "Basename")

plot = ann %>% ggplot(aes(x = X.PC1, y = Y.PC1, colour = predicted, shape = gender)) +
  geom_point() +
  theme_classic() +
  scale_shape_manual(values = c("F" = 16, "M" = 17), labels = c("Female", "Male")) +
  scale_color_manual(values = c("F" = "plum2", "M" = "skyblue2", "N" = "grey"), labels = c("Female", "Male", "Not specified")) +
  labs(shape = "Recorded sex", colour = "Predicted Sex")

# compare the annotated gender to the predicted gender

ann = sex$test %>% rownames_to_column("Basename") %>% 
  select(Basename, predicted.X, predicted.Y, predicted) %>%
  full_join(ann, ., by = "Basename")

mismatch = ann %>% filter(gender != predicted)

# samples with a mismatched sex prediction and recorded sex are removed from downstream analyses



