#################################################
## Pre filtering by total sequenced expression ##
#################################################
WorkDir <- "~/Documents/"

# Normalization to allow testing for poorly aligned 
# or otherwise outlier samples

## Libraries ##
suppressMessages(library(Matrix))
suppressMessages(library(scran))
suppressMessages(library(splines))
suppressMessages(library(lmtest))
suppressMessages(library(BiocParallel))


## Data ##
setwd(paste0(WorkDir, "/PreprocessingAndNormalization/"))
load("FormatedRawData.RData")


## Format for analysis ##
hExpr.list <- list(
  H100 = cbind(H100.1, H100.2, H100.3),
  H85 = cbind(H85.1, H85.2, H85.3),
  H10 = cbind(H10.1, H10.2, H10.3),
  H0 = cbind(H0.1, H0.2, H0.3)
)
hTime.list <- list(
  H100 = c(c(0:8, 2 * (5:21)), c(0:8, 2 * (5:15), 2 * (17:21)), c(0:8, 2 * (5:15), 2 * (17:21))),
  H85 = c(c(1:8, 2 * (5:21)), c(1:8, 2 * (5:21)), c(1:8, 2 * (5:21))),
  H10 = c(c(1:8, 2 * (5:21)), c(1:8, 2 * (5:21)), c(1:8, 2 * (5:21))),
  H0 = c(c(0:8, 2 * (5:21)), c(0:8, 2 * (5:9), 2 * (11:21)), c(0:8, 2 * (5:9), 2 * (11:21)))
)

mExpr.list <- list(
  M0 = cbind(M0.1, M0.2, M0.3),
  M15 = cbind(M15.1, M15.2, M15.3),
  M90 = cbind(M90.1, M90.2, M90.3),
  M100 = cbind(M100.1, M100.2, M100.3)
)
mTime.list <- list(
  M0 = c(c(0:8, 2 * (5:21)), c(0:8, 2 * (5:15), 2 * (17:21)), c(0:8, 2 * (5:15), 2 * (17:21))),
  M15 = c(c(1:8, 2 * (5:21)), c(1:8, 2 * (5:21)), c(1:8, 2 * (5:21))),
  M90 = c(c(1:8, 2 * (5:21)), c(1:8, 2 * (5:21)), c(1:8, 2 * (5:21))),
  M100 = c(c(0:8, 2 * (5:21)), c(0:8, 2 * (5:9), 2 * (11:21)), c(0:8, 2 * (5:9), 2 * (11:21)))
)


## Compute filtering thresholds ##
filterThresh.func <- function(hDat, mDat) {
  l10Expr <- log10(colSums(rbind(hDat, mDat)))
  quantile(l10Expr, 0.5) - 1.5 * diff(quantile(l10Expr, c(0.25, 0.75)))
}

filterObs.func <- function(hDat, mDat, tVec) {
  filterThresh <- filterThresh.func(hDat, mDat)
  remObs <- which(log10(colSums(rbind(hDat, mDat))) < filterThresh)
  if(length(remObs) > 0) {
    list(
      tVec = tVec[-remObs],
      hDat = hDat[, -remObs],
      mDat = mDat[, -remObs]
    )
  } else {
    list(
      tVec = tVec,
      hDat = hDat,
      mDat = mDat
    )
  }
}

## Filter for low expression ##
h100.out <- filterObs.func(hExpr.list$H100, mExpr.list$M0, hTime.list$H100)
hExpr.list$H100 <- h100.out$hDat
mExpr.list$M0 <- h100.out$mDat
hTime.list$H100 <- h100.out$tVec
mTime.list$M0 <- h100.out$tVec

h85.out <- filterObs.func(hExpr.list$H85, mExpr.list$M15, hTime.list$H85)
hExpr.list$H85 <- h85.out$hDat
mExpr.list$M15 <- h85.out$mDat
hTime.list$H85 <- h85.out$tVec
mTime.list$M15 <- h85.out$tVec

h10.out <- filterObs.func(hExpr.list$H10, mExpr.list$M90, hTime.list$H10)
hExpr.list$H10 <- h10.out$hDat
mExpr.list$M90 <- h10.out$mDat
hTime.list$H10 <- h10.out$tVec
mTime.list$M90 <- h10.out$tVec

h0.out <- filterObs.func(hExpr.list$H0, mExpr.list$M100, hTime.list$H0)
hExpr.list$H0 <- h0.out$hDat
mExpr.list$M100 <- h0.out$mDat
hTime.list$H0 <- h0.out$tVec
mTime.list$M100 <- h0.out$tVec

## Add FPKM expression ##
hFpkm.list <- list(
  H100 = cbind(H100.1.fpkm, H100.2.fpkm, H100.3.fpkm)[, colnames(hExpr.list$H100)],
  H85 = cbind(H85.1.fpkm, H85.2.fpkm, H85.3.fpkm)[, colnames(hExpr.list$H85)],
  H10 = cbind(H10.1.fpkm, H10.2.fpkm, H10.3.fpkm)[, colnames(hExpr.list$H10)],
  H0 = cbind(H0.1.fpkm, H0.2.fpkm, H0.3.fpkm)[, colnames(hExpr.list$H0)]
)

mFpkm.list <- list(
  M0 = cbind(M0.1.fpkm, M0.2.fpkm, M0.3.fpkm)[, colnames(mExpr.list$M0)],
  M15 = cbind(M15.1.fpkm, M15.2.fpkm, M15.3.fpkm)[, colnames(mExpr.list$M15)],
  M90 = cbind(M90.1.fpkm, M90.2.fpkm, M90.3.fpkm)[, colnames(mExpr.list$M90)],
  M100 = cbind(M100.1.fpkm, M100.2.fpkm, M100.3.fpkm)[, colnames(mExpr.list$M100)]
)


#################################
## Pre-filtering normalization ##
#################################

## Compute scale-factors ##
scaleFactor.func <- function(exprList, tList) {
  timePointSF.list <- mapply(function(subExpr, subT) {
    SF <- rep(1, ncol(subExpr))
    for(i in sort(unique(subT))) {
      tInd <- subT == i
      if(sum(tInd) == 1) {next}
      SF[tInd] <- 1 / (mean(colSums(subExpr[, tInd])) / colSums(subExpr[, tInd]))
    }
    return(SF)
  }, exprList, tList, SIMPLIFY = FALSE)
  
  pooledSamples <- mapply(function(subExpr, subT) {
    poolExpr <- matrix(0, nrow(subExpr), length(unique(subT)))
    tVec <- sort(unique(subT))
    for(i in 1:length(tVec)) {
      poolExpr[, i] <- rowMeans(subExpr[, subT == tVec[i], drop = FALSE])
    }
    return(poolExpr)
  }, exprList, tList, SIMPLIFY = FALSE)
  poolMat <- do.call(cbind, pooledSamples)
  pooledSF <- calculateSumFactors(poolMat, 
                                  sizes = 3:min(20, ncol(poolMat) - 1))
  pooledSF.list <- list()
  startInd <- 1
  for(i in 1:length(pooledSamples)) {
    endInd <- (startInd + ncol(pooledSamples[[i]]) - 1)
    pooledSF.list[[i]] <- pooledSF[startInd:endInd]
    startInd <- endInd + 1
    names(pooledSF.list)[i] <- names(pooledSamples)[i]
  }
  
  sampleSF.list <- mapply(function(subTimePointSF, subPooledSF, subT) {
    tVec <- sort(unique(subT))
    for(i in 1:length(tVec)) {
      tInd <- subT == tVec[i]
      subTimePointSF[tInd] <- subTimePointSF[tInd] * subPooledSF[i]
    }
    return(subTimePointSF)
  }, timePointSF.list, pooledSF.list, tList, SIMPLIFY = FALSE)
  sfMu <- mean(unlist(sampleSF.list))
  sampleSF.list <- lapply(sampleSF.list, function(x) {x / sfMu})
  
  return(sampleSF.list)
}

h.SF <- scaleFactor.func(hExpr.list[c("H100", "H85", "H10")], 
                         hTime.list[c("H100", "H85", "H10")])

m.SF <- scaleFactor.func(mExpr.list[c("M100", "M90", "M15")], 
                         mTime.list[c("M100", "M90", "M15")])


## Normalize ##
hNorm.list <- mapply(function(subExpr, subSF) {
  t(t(subExpr) / subSF)
}, hExpr.list[c("H100", "H85", "H10")], h.SF)

mNorm.list <- mapply(function(subExpr, subSF) {
  t(t(subExpr) / subSF)
}, mExpr.list[c("M100", "M90", "M15")], m.SF)

hNorm.list$H0 <- t(t(hExpr.list$H0) / m.SF$M100)
mNorm.list$M0 <- t(t(mExpr.list$M0) / h.SF$H100)

save(hNorm.list, mNorm.list, hExpr.list, hTime.list, mExpr.list, mTime.list,
     file = "PreFilterNorm.RData", 
     compress = "xz")


########################################
## Remove missaligned/outlier samples ##
########################################


## Libraries ##
suppressMessages(library(splines))


## Testing Functions ##
calcPC <- function(dat, nPC) {
  geneVarOrd <- order(apply(dat, 1, var), decreasing = TRUE)
  prcomp(t(as.matrix(dat[geneVarOrd[1:1000], ])), 
         center = TRUE, scale. = TRUE)$x[, 1:nPC]
}

calcSplines <- function(PC, tVec, doFDR = TRUE) {
  stdResid <- matrix(nrow = ncol(PC), ncol = length(tVec))
  for(i in 1:ncol(PC)) {
    x <- PC[, i]
    mod <- lm(x ~ bs(tVec, degree = 4))
    stdResid[i, ] <- rstudent(mod)
  }
  stdDist <- sqrt(colSums(stdResid^2))
  pVal <- apply(stdResid, 2, function(x) {
    lp <- ncol(PC) * pt(max(abs(x)), mod$df.residual - 1, lower.tail = TRUE, log.p = TRUE)
    if(-lp < 0.01) {
      2 * (-expm1(lp))
    } else {
      2 * (1 - exp(lp))
    }
  })
  if(doFDR) {pVal <- p.adjust(pVal, method = "BH")}
  return(list(
    p = pVal,
    err = stdDist
  ))
}

runSelection <- function(dat, tVec, nPC, pThresh) {
  subDat <- dat
  subT <- tVec
  remVec <- c()
  pVec <- rep(0, length(subT))
  while(any(pVec <= pThresh)) {
    splineOut <- calcSplines(calcPC(subDat, nPC), subT, doFDR = FALSE)
    pVec <- splineOut$p
    if(any(splineOut$p <= pThresh)) {
      dupTInd <- which(subT %in% unique(subT[duplicated(subT)]))
      if(!any(splineOut$p[dupTInd] <= pThresh)) {pVec <- 1; break}
      remInd <- dupTInd[which(splineOut$p[dupTInd] == min(splineOut$p[dupTInd]))]
      if(length(remInd) > 1) {
        remInd <- remInd[which.max(splineOut$err[remInd])]
      }
      if(length(remVec) > 0) {
        remVec <- c(seq_len(ncol(dat))[-remVec][remInd], remVec)
      } else {
        remVec <- remInd
      }
      subDat <- subDat[, -remInd]
      subT <- subT[-remInd]
      pVec <- 0
    }
  }
  while(TRUE) {
    if(length(remVec) == 0) {break}
    addDat <- FALSE
    subDat <- cbind(0, subDat)
    subT <- c(0, subT)
    for(i in 1:length(remVec)) {
      subDat[, 1] <- dat[, remVec[i]]
      subT[i] <- tVec[remVec[i]]
      splineOut <- calcSplines(calcPC(subDat, nPC), subT, doFDR = TRUE)
      pVec <- splineOut$p
      if(pVec[1] >= pThresh) {
        addDat <- TRUE
        subDat <- cbind(subDat, dat[, remVec[i]])
        subT <- c(subT, tVec[remVec[i]])
        remVec <- remVec[-i]
        break
      }
    }
    subDat <- subDat[, -1]
    subT <- subT[-1]
    if((!addDat) | length(remVec) == 0) {break}
  }
  if(length(remVec) == 0) {remVec <- c()}
  return(remVec)
}

subSetDat <- function(datList, dat.raw, tVec, nPC, pThresh, minMu = 1, minExp = 10) {
  dat <- do.call(rbind, datList)
  dat.raw <- do.call(rbind, dat.raw)
  geneInd <- (rowMeans(dat.raw) >= minMu) | (apply(dat.raw, 1, max) >= minExp)
  remInd <- runSelection(dat[geneInd, ], tVec, nPC, pThresh)
  return(remInd)
}

## Filter Data ##
nPC = 10; sigThresh = 1e-5
H100.filter <- subSetDat(list(hNorm.list$H100, mNorm.list$M0), 
                         list(hExpr.list$H100, mExpr.list$M0), 
                         hTime.list$H100, nPC, sigThresh)
H85.filter <- subSetDat(list(hNorm.list$H85, mNorm.list$M15), 
                        list(hExpr.list$H85, mExpr.list$M15), 
                        hTime.list$H85, nPC, sigThresh)
H10.filter <- subSetDat(list(hNorm.list$H10, mNorm.list$M90), 
                        list(hExpr.list$H10, mExpr.list$M90), 
                        hTime.list$H10, nPC, sigThresh)
H0.filter <- subSetDat(list(hNorm.list$H0, mNorm.list$M100), 
                       list(hExpr.list$H0, mExpr.list$M100), 
                       hTime.list$H0, nPC, sigThresh)

save(H100.filter, H85.filter, H10.filter, H0.filter,
     file = "FilterSamples.RData",
     compress = "xz")




##################################
## Post-filtering normalization ##
##################################

## Enumerate filtered time points ##
t.mix100 <- hTime.list$H100
if(length(H100.filter) > 0) {t.mix100 <- t.mix100[-H100.filter]}
t.mix85 <- hTime.list$H85
if(length(H85.filter) > 0) {t.mix85 <- t.mix85[-H85.filter]}
t.mix10 <- hTime.list$H10
if(length(H10.filter) > 0) {t.mix10 <- t.mix10[-H10.filter]}
t.mix0 <- hTime.list$H0
if(length(H0.filter) > 0) {t.mix0 <- t.mix0[-H0.filter]}

## Structure filtered data ##
dat.h100 <- list(allDay = hExpr.list$H100)
if(length(H100.filter) > 0) {dat.h100$allDay <- dat.h100$allDay[, -H100.filter]}
dat.h85 <- list(allDay = hExpr.list$H85)
if(length(H85.filter) > 0) {dat.h85$allDay <- dat.h85$allDay[, -H85.filter]}
dat.h10 <- list(allDay = hExpr.list$H10)
if(length(H10.filter) > 0) {dat.h10$allDay <- dat.h10$allDay[, -H10.filter]}

dat.m100 <- list(allDay = mExpr.list$M100)
if(length(H0.filter) > 0) {dat.m100$allDay <- dat.m100$allDay[, -H0.filter]}
dat.m90 <- list(allDay = mExpr.list$M90)
if(length(H10.filter) > 0) {dat.m90$allDay <- dat.m90$allDay[, -H10.filter]}
dat.m15 <- list(allDay = mExpr.list$M15)
if(length(H85.filter) > 0) {dat.m15$allDay <- dat.m15$allDay[, -H85.filter]}

## Add timepoints to data ##
dat.h100$days <- t.mix100
dat.h85$days <- t.mix85
dat.h10$days <- t.mix10
dat.m100$days <- t.mix0
dat.m90$days <- t.mix10
dat.m15$days <- t.mix85

dat.h100$days_jitter <- t.mix100 + rnorm(length(t.mix100), sd = 0.01)
dat.h85$days_jitter <- t.mix85 + rnorm(length(t.mix85), sd = 0.01)
dat.h10$days_jitter <- t.mix10 + rnorm(length(t.mix10), sd = 0.01)
dat.m100$days_jitter <- t.mix0 + rnorm(length(t.mix0), sd = 0.01)
dat.m90$days_jitter <- t.mix10 + rnorm(length(t.mix10), sd = 0.01)
dat.m15$days_jitter <- t.mix85 + rnorm(length(t.mix85), sd = 0.01)

## Add raw EC data ##
dat.h100$EC <- dat.h100$allDay
dat.h85$EC <- cbind(dat.h100$EC[, dat.h100$days == 0], dat.h85$allDay)
dat.h10$EC <- cbind(dat.h100$EC[, dat.h100$days == 0], dat.h10$allDay)

dat.m100$EC <- dat.m100$allDay
dat.m90$EC <- cbind(dat.m100$EC[, dat.m100$days == 0], dat.m90$allDay)
dat.m15$EC <- cbind(dat.m100$EC[, dat.m100$days == 0], dat.m15$allDay)

dat.h100$csDepth <- log(Matrix::colSums(dat.h100$allDay))
dat.h85$csDepth <- log(Matrix::colSums(cbind(dat.h100$allDay[, dat.h100$days==0],dat.h85$allDay)))
dat.h10$csDepth <- log(Matrix::colSums(cbind(dat.h100$allDay[, dat.h100$days==0],dat.h10$allDay)))

dat.m100$csDepth <- log(Matrix::colSums(dat.m100$allDay))
dat.m90$csDepth <- log(Matrix::colSums(cbind(dat.m100$allDay[, dat.m100$days==0],dat.m90$allDay)))
dat.m15$csDepth <- log(Matrix::colSums(cbind(dat.m100$allDay[, dat.m100$days==0],dat.m15$allDay)))

## Add fpkm data ##
dat.h100$fpkm <- hFpkm.list$H100
if(length(H100.filter) > 0) {dat.h100$fpkm <- dat.h100$fpkm[, -H100.filter]}
dat.h85$fpkm <- cbind(dat.h100$fpkm[, dat.h100$days == 0], hFpkm.list$H85)
if(length(H85.filter) > 0) {dat.h85$fpkm <- dat.h85$fpkm[, -H85.filter]}
dat.h10$fpkm <- cbind(dat.h100$fpkm[, dat.h100$days == 0], hFpkm.list$H10)
if(length(H10.filter) > 0) {dat.h10$fpkm <- dat.h10$fpkm[, -H10.filter]}

dat.m100$fpkm <- mFpkm.list$M100
if(length(H0.filter) > 0) {dat.m100$fpkm <- dat.m100$fpkm[, -H0.filter]}
dat.m90$fpkm <- cbind(dat.m100$fpkm[, dat.m100$days == 0], mFpkm.list$M90)
if(length(H10.filter) > 0) {dat.m90$fpkm <- dat.m90$fpkm[, -H10.filter]}
dat.m15$fpkm <- cbind(dat.m100$fpkm[, dat.m100$days == 0], mFpkm.list$M15)
if(length(H85.filter) > 0) {dat.m15$fpkm <- dat.m15$fpkm[, -H85.filter]}

## Normalize counts ##
hExpr.list <- list(H100 = dat.h100$allDay, H85 = dat.h85$allDay, H10 = dat.h10$allDay)
hTime.list <- list(H100 = dat.h100$days, H85 = dat.h85$days, H10 = dat.h10$days)
h.SF <- scaleFactor.func(hExpr.list, hTime.list)

mExpr.list <- list(M100 = dat.m100$allDay, M90 = dat.m90$allDay, M15 = dat.m15$allDay)
mTime.list <- list(M100 = dat.m100$days, M90 = dat.m90$days, M15 = dat.m15$days)
m.SF <- scaleFactor.func(mExpr.list, mTime.list)

dat.h100$allDay <- t(t(dat.h100$allDay) / h.SF$H100)
dat.h85$allDay <- t(t(dat.h85$allDay) / h.SF$H85)
dat.h10$allDay <- t(t(dat.h10$allDay) / h.SF$H10)

dat.m100$allDay <- t(t(dat.m100$allDay) / m.SF$M100)
dat.m90$allDay <- t(t(dat.m90$allDay) / m.SF$M90)
dat.m15$allDay <- t(t(dat.m15$allDay) / m.SF$M15)

dat.h100$SF <- h.SF$H100
dat.h85$SF <- c(h.SF$H100[dat.h100$days == 0], h.SF$H85)
dat.h10$SF <- c(h.SF$H100[dat.h100$days == 0], h.SF$H10)

dat.m100$SF <- m.SF$M100
dat.m90$SF <- c(m.SF$M100[dat.m100$days == 0], m.SF$M90)
dat.m15$SF <- c(m.SF$M100[dat.m100$days == 0], m.SF$M15)

dat.h85$allDay <- cbind(dat.h100$allDay[, dat.h100$days == 0], dat.h85$allDay)
dat.h85$days <- c(dat.h100$days[dat.h100$days == 0], dat.h85$days)
dat.h85$days_jitter <- c(dat.h100$days_jitter[dat.h100$days == 0], dat.h85$days_jitter)
dat.h85$csDepth <- c(dat.h100$csDepth[dat.h100$days == 0], dat.h85$csDepth)

dat.h10$allDay <- cbind(dat.h100$allDay[, dat.h100$days == 0], dat.h10$allDay)
dat.h10$days <- c(dat.h100$days[dat.h100$days == 0], dat.h10$days)
dat.h10$days_jitter <- c(dat.h100$days_jitter[dat.h100$days == 0], dat.h10$days_jitter)
dat.h10$csDepth <- c(dat.h100$csDepth[dat.h100$days == 0], dat.h10$csDepth)

dat.m90$allDay <- cbind(dat.m100$allDay[, dat.m100$days == 0], dat.m90$allDay)
dat.m90$days <- c(dat.m100$days[dat.m100$days == 0], dat.m90$days)
dat.m90$days_jitter <- c(dat.m100$days_jitter[dat.m100$days == 0], dat.m90$days_jitter)
dat.m90$csDepth <- c(dat.m100$csDepth[dat.m100$days == 0], dat.m90$csDepth)

dat.m15$allDay <- cbind(dat.m100$allDay[, dat.m100$days == 0], dat.m15$allDay)
dat.m15$days <- c(dat.m100$days[dat.m100$days == 0], dat.m15$days)
dat.m15$days_jitter <- c(dat.m100$days_jitter[dat.m100$days == 0], dat.m15$days_jitter)
dat.m15$csDepth <- c(dat.m100$csDepth[dat.m100$days == 0], dat.m15$csDepth)

dat.h100$sumDay <- sapply(sort(unique(dat.h100$days)), function(d) {rowMeans(dat.h100$allDay[, dat.h100$days == d, drop = FALSE])})
dat.h85$sumDay <- sapply(sort(unique(dat.h85$days)), function(d) {rowMeans(dat.h85$allDay[, dat.h85$days == d, drop = FALSE])})
dat.h10$sumDay <- sapply(sort(unique(dat.h10$days)), function(d) {rowMeans(dat.h10$allDay[, dat.h10$days == d, drop = FALSE])})

dat.m100$sumDay <- sapply(sort(unique(dat.m100$days)), function(d) {rowMeans(dat.m100$allDay[, dat.m100$days == d, drop = FALSE])})
dat.m90$sumDay <- sapply(sort(unique(dat.m90$days)), function(d) {rowMeans(dat.m90$allDay[, dat.m90$days == d, drop = FALSE])})
dat.m15$sumDay <- sapply(sort(unique(dat.m15$days)), function(d) {rowMeans(dat.m15$allDay[, dat.m15$days == d, drop = FALSE])})

dat.h100$sumDayEC <- sapply(sort(unique(dat.h100$days)), function(d) {rowSums(dat.h100$EC[, dat.h100$days == d, drop = FALSE])})
dat.h85$sumDayEC <- sapply(sort(unique(dat.h85$days)), function(d) {rowSums(dat.h85$EC[, dat.h85$days == d, drop = FALSE])})
dat.h10$sumDayEC <- sapply(sort(unique(dat.h10$days)), function(d) {rowSums(dat.h10$EC[, dat.h10$days == d, drop = FALSE])})

dat.m100$sumDayEC <- sapply(sort(unique(dat.m100$days)), function(d) {rowSums(dat.m100$EC[, dat.m100$days == d, drop = FALSE])})
dat.m90$sumDayEC <- sapply(sort(unique(dat.m90$days)), function(d) {rowSums(dat.m90$EC[, dat.m90$days == d, drop = FALSE])})
dat.m15$sumDayEC <- sapply(sort(unique(dat.m15$days)), function(d) {rowSums(dat.m15$EC[, dat.m15$days == d, drop = FALSE])})

dat.h100$sumDayFpkm <- sapply(sort(unique(dat.h100$days)), function(d) {rowSums(dat.h100$fpkm[, dat.h100$days == d, drop = FALSE])})
dat.h85$sumDayFpkm <- sapply(sort(unique(dat.h85$days)), function(d) {rowSums(dat.h85$fpkm[, dat.h85$days == d, drop = FALSE])})
dat.h10$sumDayFpkm <- sapply(sort(unique(dat.h10$days)), function(d) {rowSums(dat.h10$fpkm[, dat.h10$days == d, drop = FALSE])})

dat.m100$sumDayFpkm <- sapply(sort(unique(dat.m100$days)), function(d) {rowSums(dat.m100$fpkm[, dat.m100$days == d, drop = FALSE])})
dat.m90$sumDayFpkm <- sapply(sort(unique(dat.m90$days)), function(d) {rowSums(dat.m90$fpkm[, dat.m90$days == d, drop = FALSE])})
dat.m15$sumDayFpkm <- sapply(sort(unique(dat.m15$days)), function(d) {rowSums(dat.m15$fpkm[, dat.m15$days == d, drop = FALSE])})

## Save results ##
save(dat.h100,dat.h85,dat.h10,
     dat.m100,dat.m90,dat.m15,
     file="NormDat.RData",
     compress="xz")


