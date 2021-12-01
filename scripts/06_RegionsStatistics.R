###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: statistics testing of CpG sites aggregated into regions (genes, promoters and CpG islands)

# input: preprocessed beta values generated using 01preprocessingSWAN.R to calculate average methylation values per gene
#        aggregated beta values per promoter and CpG island extracted frmo the output of 01preprocessingSWAN.R
#        conducted to compare G5 and HTF within the IVF cohort (including and excluding pregnancy complications)
#        and once to compare all IVF to all naturally conceived individuals for each region type (genes, promoters, CpG islands)

# targeted analysis: for the targeted analysis, genes were restricted to those contained in the following documents:
#                     imprinted sites: manifestImprinted.csv (source: Ginjala    V. Gene imprinting gateway. Genome Biology 2001.)

# output: (multiple testing corrected) statistics per region - gene/promoter/CGI (.csv)

###########################################################################################################################

# load packages

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(variancePartition))
suppressPackageStartupMessages(library(reshape2))

# set up the parallel environment

param = SnowParam(workers = 30, type = "SOCK", progressbar = TRUE)
register(param)

# set data directories

betas.file <- "betas.csv"

annotation.file <- "annotation.csv"

manifest.file <- "manifest_EPIC_selected_combined.csv"

# load the data files

annotation <- read.csv(annotation.file)

betas <- fread(betas.file) %>% 
	select(ID, as.character(annotation$Sample_ID)) %>% 
	column_to_rownames("ID")

manifest <- fread(manifest.file) %>% select(ID, CHR, custom_Name, custom_Accession, custom_Group)

# adjust the sample plate variable in the annotation file

annotation$samplePlate <- if_else(annotation$Sample_Plate == "WG5839007-BCD", "a", "b")

annotation$complication <- if_else(annotation$pregnancy_complic == "no", "no", "yes")

# for consistency of naming rownames of annotation file should correspond to colnames of betas file

annotation <- annotation %>% column_to_rownames("Sample_ID")

###############################################################################################################
###################### group the sites by gene & calculate the average methylation level per gene ##############

# filter the manifest file to contain only probes passing QC

manifest <- manifest %>% filter(ID %in% rownames(betas))

# annotate the betas file with the corresponding information from the manifest file

betas.annotated <- betas %>%
        rownames_to_column("ID") %>%
        full_join(manifest, . , by = "ID")

# select valid gene names (genes uniquely mapping to 1 chromosome and containing at least 3 CpG sites)

sites_per_gene <- manifest %>% 
	filter(nchar(custom_Name) > 0) %>% 
	group_by(custom_Name, CHR) %>% 
	summarise(sites = n())

chroms_per_gene <- sites_per_gene %>% 
	group_by(custom_Name) %>% 
	summarise(chrs = n())

#remove genes that aren't annotated to a unique chromosome
chroms_per_gene <- chroms_per_gene %>% filter(chrs == 1)

# remove genes that aren't represented by 3 or more individual CpG sites
sites_per_gene <- sites_per_gene %>% filter(custom_Name %in% chroms_per_gene$custom_Name) %>% filter(sites >=3)

# determine the number of sites per region (within each gene name
# remove unannotated sites and genes

sites.per.region <- betas.annotated %>% 
	filter(nchar(custom_Name) > 0, nchar(custom_Group) > 0) %>% 
	filter(custom_Name %in% sites_per_gene$custom_Name) %>% 
	group_by(custom_Name, custom_Group) %>% 
	tally()

# only keep sites with 3 or more CpG sites

sites.per.region <- sites.per.region %>% 
	filter(n >= 3)

# make filtering variables

sites.per.region[ , "filtering_name"] <- paste0(sites.per.region$custom_Name, "_._", sites.per.region$custom_Group)

betas.annotated[ , "filtering_name"] <- paste0(betas.annotated$custom_Name, "_._", betas.annotated$custom_Group)

# filter the betas file to contain only the relevant sites & calculate the grouped means

betas.regions <- betas.annotated %>% 
	filter(filtering_name %in% sites.per.region$filtering_name) %>% 
	select(-ID, -CHR, -custom_Name, -custom_Accession, - custom_Group) %>% 
	group_by(filtering_name) %>% 
	summarise_all(., list(~mean(.,na.rm = T)))

betas.regions <- colsplit(betas.regions$filtering_name, pattern = "_._", 
				names = c("custom_Name", "custom_Group")) %>%
	cbind(betas.regions)

# create a vector of unique regions types

regions <- unique(betas.regions$custom_Group)


###############################################################################################################
################################ conduct statistical testing #################################################

# specify the model

model <- ~ MEDIUM + counts.CD8T + counts.CD4T + counts.NK + counts.Bcell + counts.Mono + 
  counts.Gran + counts.nRBC + gender + gestational_age + maternal_age + 
  complication + (1 | samplePlate) + (1 | center)

complex.regions <- list()

# to apply the model to all region types the regional methylation values should be combined into one
# data frame with an extra column "regions" denoting the region type

for(region in regions){
	
	data = betas.regions %>% filter(custom_Group == region) %>%
		select(-filtering_name, -custom_Group) %>%
		column_to_rownames("custom_Name")

	weights = voomWithDreamWeights(na.omit(data), model, annotation)	

	fit = dream(weights, model, annotation)

	regions[[region]] <- topTable(fit.complex, coef = "MEDIUMvg5", number = Inf)

	regions[[region]]$region <- region

}

# create 1 data frame for each output list

regions.complete <- rbind(regions[[1]], regions[[2]], regions[[3]],
				 regions[[4]], regions[[5]], regions[[6]],
				 regions[[7]])

write.csv(x =  regions.complete, file = "dmrsRegionsFitDream.csv")
