# apply the iEVORA algorithm to preprocessed beta values

# input: preprocessed beta values generated using 01preprocessingSWAN.R
#        conducted to compare G5 and HTF within the IVF cohort (including and excluding pregnancy complications)

# load packages
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(matrixTests))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(ggplot2))

# load the annotation & betas files 

ann = fread("annotation.csv") %>%
  mutate(Sample_ID = as.character(Sample_ID)) %>%
  select(-V1)

betas = fread("betas.csv") %>%
  select(ID, all_of(ann$Sample_ID)) %>%
  column_to_rownames("ID") %>%
  .[ , order(colnames(.))] %>%
  as.matrix()

# ensure that the annotation has the same sample order & create grouping variable

ann = ann %>% mutate(groups = if_else(MEDIUM == "vg5", 0, 1)) %>%
  arrange(Sample_ID)

# group 0 = vg5
# group 1 = htf

# apply row_iEVORA

result = row_ievora(betas, ann$groups, cutT = 0.05, cutBfdr = 0.001)

result = result %>% 
  rownames_to_column("ID") %>%
  arrange(rank)

result = result %>% rename("obs.vg5" = obs.0, "obs.htf" = obs.1, 
                           "mean.vg5" = mean.0, "mean.htf" = mean.1, 
                           "var.vg5" = var.0, "var.htf" = var.1)

result = result %>% mutate(var_diff = var.vg5 - var.htf)

# of the significant diffVar sites see how many are more variable with which culture medium

sig = result %>% filter(significant == TRUE)

table(sig$var.vg5 > sig$var.htf)

# annotate the significant sites with their gene names

manifest = fread("manifest_EPIC_selected_combined.csv") %>%
  filter(ID %in% sig$ID) %>%
  select(ID, custom_Name)

sig = full_join(sig, manifest, by = "ID")

# generate a list of gene names from the sites with significantly different variance

gene_names = sig %>% filter(nchar(custom_Name) > 0) %>% select(custom_Name)

write.csv(file = "diffVarGenes.csv", gene_names, row.names = F)

# check if there is any overlap with imprinted genes

manifest = fread("manifest_EPIC_imprinted.csv") %>%
  select(-V1)

table(sig$custom_Name %in% manifest$custom_Name)

sig %>% filter(custom_Name %in% manifest$custom_Name)

# check if there is any overlap with birthweight associated CpG sites (meta-EWAS)

birthweight = fread("dmpsComplexBirthwtDreamEWAS.csv") %>%
  rename("ID" = V1) %>% select(ID)

sig %>% filter(ID %in% birthweight$ID) %>% nrow()