###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: Visualisation of the output of 05EWAS.R and 06RegionStatistics.R

# input: preprocessed beta values corresponding to the conducted statistical testing (output from 01preprocessingSWAN.R)
#        statistical testing results generated using 05EWAS.R or 06RegionStatistics.R
#        manifest file containing sites of interest (imprinted sites/genes & birth-weight associated sites & genes)

# output: volcano plot(s)

###########################################################################################################################

#load packages
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(qqman))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(VennDiagram))

# set data directories

base.dir <- "output"

annotation.file <- file.path(base.dir, "annotation.csv")

dmps.file <- file.path(base.dir, "dmpsFitDream.csv")

betas.file <- file.path(base.dir, "betas.csv")

# load the data

annotation <- fread(annotation.file) %>% mutate(Sample_ID = as.character(Sample_ID)) %>% select(-V1)

dmps <- fread(dmps.file) %>% rename("ID" = "V1")

betas.annotation <- fread(betas.file) %>% 
  select(ID, Start)

# work out the group means and differences in beta values

vg5 <- annotation %>% filter(MEDIUM == "vg5") %>% .$Sample_ID
htf <- annotation %>% filter(MEDIUM == "htf") %>% .$Sample_ID

group_betas <- betas

group_betas$htf_mean <- group_betas %>% select(all_of(htf)) %>% rowMeans(na.rm = TRUE)
group_betas$vg5_mean <- group_betas %>% select(all_of(vg5)) %>% rowMeans(na.rm = TRUE)

group_betas$htf_SD <- group_betas %>% select(all_of(htf)) %>% as.matrix() %>% rowSds(na.rm = TRUE)
group_betas$vg5_SD <- group_betas %>% select(all_of(vg5)) %>% as.matrix() %>% rowSds(na.rm = TRUE)

group_betas <- group_betas %>% select(-all_of(htf), -all_of(vg5)) %>%
  rownames_to_column("ID")

#use vg5 group as reference group
group_betas <- group_betas %>% mutate(mean_diff = vg5_mean - htf_mean)

# join the data:

dmps <- full_join(dmps, group_betas, by = "ID")

# draw a volcano plot of the results (use p.value threshold of 0.1)

dmps.volcano <- dmps %>% mutate(significant = adj.P.Val < 0.1)

plot <- dmps.volcano %>% 
  ggplot(aes(x = mean_diff, y = -(log10(P.Value)), colour = significant)) +
  geom_point() +
  xlab("mean difference (mean G5 - mean HTF)") +
  ylab("-log10 p-value") +
  labs(colour = "FDR significance") +
  scale_color_manual(labels = c("Not significant (p > 0.1)", "p < 0.1"), values = c("grey", "red")) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom") +
  scale_x_continuous(breaks = c( -0.2, -0.1, 0, 0.1, 0.2 ), 
                     limits = c(-0.25, 0.2))

####################### load the data for imprinted sites and birthweight sites to get position IDs ##############

bwt= fread("manifestBirthweight.csv") 

bwt = bwt$ID

imp = fread("manifestImprinted.csv") 

imp = imp$ID

# add columns to the original data to signify which positions/genes are associated with birthweight/imprinting genes

dmps.volcano = dmps.volcano %>% mutate(bwt = ID %in% bwt,
                                                       imp = ID %in% imp)

# re-draw the volcano plot highlighting the sites / genes of interest (specify whether to include imp/bwt sites/genes)

plot.imprinted = dmps.volcano %>% ggplot(aes(x = mean_diff, y = -(log10(P.Value)))) +
  geom_point(colour = "grey", size = 0.4) +
  geom_point(data = function(x){filter(x, imp)}, colour = "deeppink4", size = 0.4) +
  xlab("mean difference \n (mean G5 - mean HTF)") +
  ylab("-log10 p-value") +
    theme_classic(base_size = 12) +
   scale_x_continuous(breaks = c( -0.2, -0.1, 0, 0.1, 0.2 ), 
                      limits = c(-0.25, 0.2))
 