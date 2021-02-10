###############
## Run Trend ##
###############
WorkDir <- "~/Documents/"
setwd(WorkDir)


###############
## Libraries ##
###############
suppressMessages(library(Trendy))


##########
## Data ##
##########
load("FACS_Data/NormDat.RData")


##################
## Filter Genes ##
##################
## Remove genes with 80th percentile expression across experiments < 20 ##
hGenes <- sort(rownames(h100.facs33$normDat)[(pmax(
  apply(h100.facs33$normDat, 1, quantile, probs = 0.8),
  apply(hMix.facs33$normDat, 1, quantile, probs = 0.8)
)) >= 20])
mGenes <- sort(rownames(m100.facs33$normDat)[(pmax(
  apply(m100.facs33$normDat, 1, quantile, probs = 0.8),
  apply(mMix.facs33$normDat, 1, quantile, probs = 0.8)
)) >= 20])


################
## Run Trendy ##
################
MK <- 4 # maximum number of breakpoints
mN <- 2 # minimum number of points in segment
nIter <- 20 # number of seeds to try
pValThresh <- 0.1 # p-value at which to call slope
nCores <- parallel::detectCores()

h100FACS33.trendy <- results(trendy(h100.facs33$normDat[hGenes, ],
                                    h100.facs33$Day,
                                    meanCut = 0, maxK = MK, minNumInSeg = mN,
                                    pvalCut = pValThresh, numTry = nIter,
                                    NCores = nCores))
hMixFACS33.trendy <- results(trendy(hMix.facs33$normDat[hGenes, ],
                                    hMix.facs33$Day,
                                    meanCut = 0, maxK = MK, minNumInSeg = mN,
                                    pvalCut = pValThresh, numTry = nIter,
                                    NCores = nCores))

m100FACS33.trendy <- results(trendy(m100.facs33$normDat[mGenes, ],
                                    m100.facs33$Day,
                                    meanCut = 0, maxK = MK, minNumInSeg = mN,
                                    pvalCut = pValThresh, numTry = nIter,
                                    NCores = nCores))
mMixFACS33.trendy <- results(trendy(mMix.facs33$normDat[mGenes, ],
                                    mMix.facs33$Day,
                                    meanCut = 0, maxK = MK, minNumInSeg = mN,
                                    pvalCut = pValThresh, numTry = nIter,
                                    NCores = nCores))


##################
## Save Results ##
##################
save(h100FACS33.trendy, hMixFACS33.trendy, hGenes,
     m100FACS33.trendy, mMixFACS33.trendy, mGenes,
     file = "FACS_Data/trendyFit.RData",
     compress = "xz")




