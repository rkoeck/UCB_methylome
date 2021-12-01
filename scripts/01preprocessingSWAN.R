# preprocessing of raw .idat files using RnBeads

# input: data from each cohort / array type processed separately: IVF data (GSE189531), FLEHS 450k (GSE110128), 
#                                                                 ENVIRONAGE 450K (GSE151042), ENVIRONAGE EPIC (unpublished)

# output:  methylation beta value per site per sample (format: .csv), 
#         aggregated methylation values for promoters and CpG islands per sample (format: .csv)


rm(list=ls(all=T))

# load packages
suppressPackageStartupMessages(library(RnBeads))

#set up parallel processing
num.cores <- 16

parallel.setup(num.cores)

#set file directories

data.dir <- "Projects/CultMed"

idat.dir <- file.path(data.dir, "idat/")

sample.annotation <- file.path("annotation.txt")

analysis.dir <- "cohort2"

report.dir <- file.path(analysis.dir, "preprocessing")

#name the Cohort

Cohort = "CultureMedia_Cohort2"

#set the options for RnBeads

rnb.options(analysis.name = "CultureMedia_Cohort2_preprocessing", #name the analysis
            logging = TRUE, #creates a log in the automatic run of the pipeline
            assembly = "hg19", #assembly genome
            analyze.sites = TRUE, #analyse per site/probe - always done for preprocessing steps
            identifiers.column = "Sample_ID", # column name in table of phenotype information to use as sample identifiers otherwise rownames
            gz.large.files = TRUE, #large files should be compressed
            import = TRUE, #carry out import module, only false if using previous RnBSet
            import.sex.prediction = TRUE, #does sex prediction when data is imported
            qc = TRUE, # QC module carried out
            preprocessing = TRUE, #preprocessing steps are completed
            normalization = TRUE, #normalisation is never carried out on sequencing data
            normalization.method = "swan", #select normalisation method
            normalization.background.method = "none", #no background normalisation is carried out for EPIC arrays
            filtering.context.removal = c("CC", "CAG", "CAH", "CTG", "CTG"), #only retain CpG probes
            filtering.snp = "any", #remove all probes that overlap with SNPs
            filtering.greedycut = TRUE, #run greedycut filtering to remove low quality probes and samples
            filtering.sex.chromosomes.removal = TRUE, #all probes on sex chromosomes are removed
            filtering.missing.value.quantile = 0.05, # proportion of samples that must have value for probe for it to be included
            imputation.method = "none", #no imputation to be carried out,
            inference = FALSE, #no covariate inference to be done
            exploratory = FALSE, #carry out some steps of the exploratory module
            differential = FALSE, #differential module not to be carried out,
            export.to.bed = FALSE, #export data to bed file
            export.to.trackhub = NULL, #disable export to trackhub (disc full??)
            export.to.csv = TRUE #methylation values are exported to csv files
)

options(fftempdir= analysis.dir)
options(ffcaching="ffeachflush")

rnb.run.analysis(dir.reports=report.dir, sample.sheet=sample.annotation,data.dir=idat.dir, data.type="infinium.idat.dir")

# samples and sites deemed poor quality by this preprocessing are removed from subsequent analyses