################
## Run Trendy ##
################


###############
## Libraries ##
###############
suppressMessages(library(Matrix))
suppressMessages(library(Trendy))


##########
## Data ##
##########
setwd("~/Documents/Research/Chris_Mixing_2019/Final_MixData_March25/")
load("PreprocessingAndNormalization/NormDat.RData")


##################
## Filter Genes ##
##################
## Remove genes with 80th percentile expression across experiments < 20 ##
hGenes <- sort(rownames(dat.h100$allDay)[(pmax(
    apply(dat.h100$allDay,1,quantile,probs=0.8),
    apply(dat.h85$allDay,1,quantile,probs=0.8),
    apply(dat.h10$allDay,1,quantile,probs=0.8)
  ))>=20])
mGenes <- sort(rownames(dat.m100$allDay)[(pmax(
    apply(dat.m100$allDay,1,quantile,probs=0.8),
    apply(dat.m90$allDay,1,quantile,probs=0.8),
    apply(dat.m15$allDay,1,quantile,probs=0.8)
  ))>=20])


################
## Run Trendy ##
################
MK <- 4 # maximum number of breakpoints
mN <- 5 # minimum number of points in segment
nIter <- 1000 # number of seeds to try
pValThresh <- 0.1 # p-value at which to call slope
nCores <- parallel::detectCores()

h100.trendy <- results(trendy(dat.h100$allDay[hGenes,],dat.h100$days_jitter,
                              meanCut=0,maxK=MK,minNumInSeg=mN,
                              pvalCut=pValThresh,numTry=nIter,NCores=nCores))
h85.trendy <- results(trendy(dat.h85$allDay[hGenes,],dat.h85$days_jitter,
                             meanCut=0,maxK=MK,minNumInSeg=mN,
                             pvalCut=pValThresh,numTry=nIter,NCores=nCores))
h10.trendy <- results(trendy(dat.h10$allDay[hGenes,],dat.h10$days_jitter,
                             meanCut=0,maxK=MK,minNumInSeg=mN,
                             pvalCut=pValThresh,numTry=nIter,NCores=nCores))

m100.trendy <- results(trendy(dat.m100$allDay[mGenes,],dat.m100$days_jitter,
                              meanCut=0,maxK=MK,minNumInSeg=mN,
                              pvalCut=pValThresh,numTry=nIter,NCores=nCores))
m90.trendy <- results(trendy(dat.m90$allDay[mGenes,],dat.m90$days_jitter,
                             meanCut=0,maxK=MK,minNumInSeg=mN,
                             pvalCut=pValThresh,numTry=nIter,NCores=nCores))
m15.trendy <- results(trendy(dat.m15$allDay[mGenes,],dat.m15$days_jitter,
                             meanCut=0,maxK=MK,minNumInSeg=mN,
                             pvalCut=pValThresh,numTry=nIter,NCores=nCores))


##################
## Save Results ##
##################
save(h100.trendy,h85.trendy,h10.trendy,hGenes,
     m100.trendy,m90.trendy,m15.trendy,mGenes,
     file="Trendy/trendyFit.RData",
     compress="xz")





