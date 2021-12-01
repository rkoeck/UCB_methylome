#NormalizeAndPredictGA.R
#This script is slightly modified from https://github.com/akknight/PredictGestationalAge

# namely: adjustments are made to allow the prediction to work with EPIC data (namely data that is missing some of the required CpGs)

# method replicated from: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1068-z 


fastImputation=FALSE

#STEP 1: DEFINE QUALITY METRICS

meanMethBySample =as.numeric(apply(as.matrix(dat1[,-1]),2,mean,na.rm=TRUE))
minMethBySample   =as.numeric(apply(as.matrix(dat1[,-1]),2,min,na.rm=TRUE))
maxMethBySample  =as.numeric(apply(as.matrix(dat1[,-1]),2,max,na.rm=TRUE))


datMethUsed= t(dat1[,-1])
colnames(datMethUsed)=as.character(dat1[,1])


noMissingPerSample=apply(as.matrix(is.na(datMethUsed)),1,sum)
table(noMissingPerSample)

#STEP 2: Imputing 
if (! fastImputation & nSamples>1 & max(noMissingPerSample,na.rm=TRUE)<3000 ){
  
  # run the following code if there is at least one missing
  if ( max(noMissingPerSample,na.rm=TRUE)>0 ){
    dimnames1=dimnames(datMethUsed)
    datMethUsed= data.frame(t(impute.knn(t(datMethUsed))$data))
    dimnames(datMethUsed)=dimnames1
  } # end of if
} # end of if (! fastImputation )

if ( max(noMissingPerSample,na.rm=TRUE)>=3000 ) fastImputation=TRUE


if ( fastImputation | nSamples==1 ){
  noMissingPerSample=apply(as.matrix(is.na(datMethUsed)),1,sum)
  table(noMissingPerSample)
  if ( max(noMissingPerSample,na.rm=TRUE)>0 & max(noMissingPerSample,na.rm=TRUE) >= 3000 ) {normalizeData=FALSE}
  
  # run the following code if there is at least one missing
  if ( max(noMissingPerSample,na.rm=TRUE)>0 & max(noMissingPerSample,na.rm=TRUE) < 3000 ){
    dimnames1=dimnames(datMethUsed)
    for (i in which(noMissingPerSample>0) ){
      selectMissing1=is.na(datMethUsed[i,])
      datMethUsed[i,selectMissing1] = as.numeric(probeAnnotation21kdatMethUsed$goldstandard2[selectMissing1])
    } # end of for loop
    dimnames(datMethUsed)=dimnames1
  } # end of if
} # end of if (! fastImputation )

# STEP 3: Data normalization (each sample requires about 8 seconds). It would be straightforward to parallelize this operation.

if (normalizeData ){
  datMethUsedNormalized=BMIQcalibration(datM=datMethUsed,goldstandard.beta= probeAnnotation21kdatMethUsed$goldstandard2,plots=FALSE)
}
if (!normalizeData ){ datMethUsedNormalized=datMethUsed }
rm(datMethUsed); gc()

#STEP 4: Predict age and create a data frame for the output (referred to as datout)
selectCpGsClock = datClock %>% filter(datClock$CpGmarker %in% dimnames(datMethUsedNormalized)[[2]]) %>% .[ ,1]

# adjust the normalised object to contain only CpG sites that are part of datClock

datMethUsedNormalized <- rownames_to_column(datMethUsedNormalized, "samples")
datMethClock0 = datMethUsedNormalized %>% select(samples, all_of(selectCpGsClock))

# object must contain all of the probes from datClock otherwise the prediction model will fail. 7 of these probes have been removed from the EPIC array
# as per advice for how to adjust the method to the EPIC array - the probes are included and filled with NA

removed = datClock %>% filter(!CpGmarker %in% selectCpGsClock) %>% .[-1,1]

datMethClock0 = datMethClock0 %>% mutate(cg02941816 = NA, cg19564877 = NA, cg20394284 = NA,
                                         cg23858360 = NA, cg25306927 = NA, cg25374854 = NA, cg26656135 = NA)

datMethClock0 = column_to_rownames(datMethClock0, "samples")

datMethClock= data.frame(datMethClock0[ as.character(datClock$CpGmarker[-1])])
dim(datMethClock)

#create a vector of coefficients that contains coefficients * nSamples
coefs <- rep(datClock$CoefficientTraining[-1], times = nSamples)

# turn these into a matrix where each row contains the coefficient corresponding to the relevant probe column in the new data file
coefMatrix <- t(matrix(coefs, ncol = nSamples))

#multiply the coefmatrix and the sample matrix (using the * to multiply does an element-wise multiplication ie 1,1 * 1,1, 1,2*1,2 etc)

predictedAgeMatrix <- as.matrix(datMethClock) * coefMatrix

predictedAge = data.frame(sample = rownames(predictedAgeMatrix), 
                          predictedGA = as.numeric(datClock$CoefficientTraining[1]) + rowSums(predictedAgeMatrix, na.rm = T))

datout= predictedAge %>% mutate(NAs = noMissingPerSample, meanMethSample = meanMethBySample, minMethSample = minMethBySample, maxMethSample = maxMethBySample)
