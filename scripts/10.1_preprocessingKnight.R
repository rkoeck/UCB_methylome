###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: preprocessing with RnBeads for gestational age prediction (Knight)

# input: .idat files from IVF samples

# output: beta values (.csv)

# method replicated from: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1068-z 

###########################################################################################################################

rm(list=ls(all=T))

#attach libraries
suppressPackageStartupMessages(library(RnBeads))

#set up parallel processing
num.cores <- 16

parallel.setup(num.cores)

#set file directories

data.dir <- "CultMed"

idat.dir <- file.path(data.dir, "idat/")

sample.annotation <- file.path("annotation.txt")

analysis.dir <- "cohort2"

report.dir <- file.path(analysis.dir, "preprocessingGA")

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
            import.sex.prediction = FALSE, #does sex prediction when data is imported
            qc = FALSE, # QC module carried out
            preprocessing = TRUE, #preprocessing steps are completed
            normalization = NULL, #normalisation is never carried out on sequencing data
	    normalization.method = "none", #no normalisisation method to be applied
	    normalization.background.method = "none", #no background correction applied
            filtering.greedycut = TRUE, #run greedycut filtering to remove low quality probes and samples
            filtering.greedycut.pvalue.threshold = 0.05, #detection p=value threhsold for greedycut
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

