# UCB_methylome
In this project we have compared the umbilical cord blood (UCB) methylome of neonates that were born after in vitro fertilisation (IVF) procedures using different IVF culture media, namely G5 (Vitrolife) and HTF (Lonza). Data were generated using Illumina’s Infinium MethylationEPIC bead chip that profiles around 850,000 CpG sites across the human genome. Furthermore, data from birth cohorts (ENVIRONAGE and FLEHS) were included to compare the methylomes of IVF and naturally conceived neonates.

The manifest_filter directory contains the manifest files used for the targeted analyses and the scripts directory contains the scripts used for data processing, analysis and visualization. A brief description of the contents of the scripts is given below:

### Singularity container
The build files can be used to create singularity images containing all the necessary software and packages for the following purposes:
1. [Data processing](scripts/RnBeadsBlood_build.txt)
2. [statistical testing](scripts/dream_build.txt) with mixed effects linear models as implemented in the dream package

### Data preprocessing:

1. [Data were first pre-processed using the RnBeads package which applied the following functions:](scripts/01_preprocessingSWAN.R)
    + Normalisation using the SWAN method  
    + greedycut (p-value threshold 0.05) removal of poor quality sites and samples
    + removal of probes on the sex chromosomes
    + removal of probes containing SNPs
    + removal of probes not in the CpG context
    + removal of probes with missing values in 5% or more of the samples
2.	Sex estimation using the method from the sEst package:
    + [extraction of raw beta values (unnormalized) and of detection p-values for all samples](scripts/02.1_preprocessingSest.R)
    + [sex prediction & visualisation](scripts/02.2_SexPredictionsEst.R)
3. Cell composition estimates are calculated using the IDOL optimized probes from the package FlowSorted.Blood.Epic:
    + [calculation of the cell compositions](scripts/03.1_CellCompositionEstimation.R)
    + [visualization of the results](scripts/03.2_CellCompositionVisualisation.R)

### Analysis & visualization:

4.	[Principal component analysis & determination of associations between PCs and sample characteristics](scripts/04_PCACovariates.R)
5.	[Differential methylation analysis (sites) using eBayes moderated mixed effects linear models as implemented in the VariancePartition package](scripts/05_EWAS.R)
6.	[Differential methylation analysis on a regional level using eBayes moderated mixed effects linear models as implemented in the VariancePartition package. Groups of interest: genes, promoters, CpG islands](scripts/06_RegionsStatistics.R)
7.	[Visualisation of the differential methylation analyses as volcano plots](scripts/07_VolcanoPlot.R)
8.	Determination of outliers using a threshold-based approach
    + [Calculation of the outliers](scripts/08.1_Outliers.R)
    + [visualization of the outliers per sample](scripts/08.2_OutliersVisualisation.R)
9.	[Application of the iEVORA algorithm to identify differentially variable CpG sites](scripts/09_iEVORA.R)
10.	Application of the Knight epigenetic gestational age predictor
    + [Preparation of the data according to the Knight method](scripts/10.1_preprocessingKnight.R)
    + [A script to conduct the epigenetic gestational age prediction for each sample which calls](scripts/10.2_KnightGAEstimates.R)
    + [the script for normalization and estimation](scripts/10.3_NormalizeAndPredictGA.R)
11.	Calculation of epigenetic gestational age based on the Bohlin method
    + [preprocessing of the data according to the Bohlin method](scripts/11.1_preprocessingBohlin.R)
    + [estimating epigenetic gestational age according to the Bohlin method](scripts/11.2_BohlinEstimateGA.R)
    + [calculating the residuals and visualizing the results](scripts/11.3_BohlinResiduals.R)
12.	[Combination of the resultant plots from the scripts above for publication](scripts/12_combinePlots.R)

