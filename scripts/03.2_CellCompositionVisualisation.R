###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: visualisation of the cell count esimates generated using the 03CellCompositionEstimation.R script

# input: cell composition file generated as output from 03CellCompositionEstimation.R
#        visualisation conducted separately for IVF samples and naturally conceived samples (all 3 cohorts together)

# output: violin plot overlaid with box plot & individual data points

###########################################################################################################################

#load packages
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))

#set directories

base.dir <- "output"

annotation.file <- "annotation.csv"

idol.file <- file.path(base.dir, "cellCompositionIDOL.csv")

# load the files

annotation <- fread(annotation.file) %>% select(-contains("counts"))

idol <- fread(idol.file)

# make a corresponding column in the annotation file with matching sample ID name

annotation$V1 <- annotation$barcode

# merge the annotation data with each of the cell estimation datasets

idol <- idol %>% filter(V1 %in% annotation$V1) %>%
  full_join(annotation, ., by = "V1")

# reformat the data for plotting

idol.plot <- idol %>% select(barcode, MEDIUM, contains("counts")) %>%
  gather(., cellType, counts, counts.CD8T:counts.nRBC, factor_key = T)

# draw a violin plot showing the distribution of the different cell types

idol.plot$Medium = str_replace(string = idol.plot$Medium, pattern = "vg5", replacement = "G5")
idol.plot$Medium = str_replace(string = idol.plot$Medium, pattern = "htf", replacement = "HTF")
idol.plot$cellType = str_remove(idol.plot$cellType, "counts.")


plot = idol.plot %>% ggplot(aes(x = cellType, y = counts, group = interaction(Medium, cellType))) + 
  geom_violin(aes(fill = Medium), position = position_dodge(width = 1)) +
  geom_boxplot(position = position_dodge(width = 1), width = 0.2) +
  geom_point(position = position_dodge(width = 1), size = 0.2) +
  theme_classic() +
  coord_flip() +
  scale_fill_manual(values = c(G5 = "darkgoldenrod3", HTF = "deepskyblue2")) +
  xlab("cell type") +
  ylab("proportional composition")

