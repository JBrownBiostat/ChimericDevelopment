###############################
## Alignment quality control ##
###############################


## Libraries ##
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))


## Data ##
setwd("~/Documents/Research/Chris_Mixing_2019/Final_MixData_March25/")
load("PreprocessingAndNormalization/FormatedRawData.RData")
load("PreprocessingAndNormalization/NormDat.RData")


#####################################
## Enumerate filtered observations ##
#####################################
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

H100.FT <- filterThresh.func(hExpr.list$H100, mExpr.list$M0)
H85.FT <- filterThresh.func(hExpr.list$H85, mExpr.list$M15)
H10.FT <- filterThresh.func(hExpr.list$H10, mExpr.list$M90)
H0.FT <- filterThresh.func(hExpr.list$H0, mExpr.list$M100)


## Removed observations ##
H100.remInd <- which(!(colnames(hExpr.list$H100) %in% colnames(dat.h100$allDay)))
H85.remInd <- which(!(colnames(hExpr.list$H85) %in% colnames(dat.h85$allDay)))
H10.remInd <- which(!(colnames(hExpr.list$H10) %in% colnames(dat.h10$allDay)))
H0.remInd <- which(!(colnames(hExpr.list$H0) %in% colnames(dat.m100$allDay)))


#####################################
## Objective alignment error rates ##
#####################################
## Functions ##
ErrRates.func <- function(pureDat, errorDat) {
  useGenes.pure <- toupper(rownames(pureDat)) %in% toupper(rownames(errorDat))
  useGenes.pure <- rownames(pureDat)[useGenes.pure]
  useGenes.pure <- useGenes.pure[order(toupper(useGenes.pure))]
  useGenes.pure <- useGenes.pure[!duplicated(toupper(useGenes.pure))]
  useGenes.error <- rownames(errorDat)[toupper(rownames(errorDat)) %in% toupper(useGenes.pure)]
  useGenes.error <- useGenes.error[order(toupper(useGenes.error))]
  useGenes.error <- useGenes.error[!duplicated(toupper(useGenes.error))]
  
  pureSub <- pureDat[useGenes.pure, ]
  pureSub <- pureSub[!duplicated(rownames(pureSub)), ]
  errorSub <- errorDat[useGenes.error, ]
  errorSub <- errorSub[!duplicated(rownames(errorSub)), ]
  errorSub <- errorSub[rowSums(pureSub > 0) >= 5, ]
  pureSub <- pureSub[rowSums(pureSub > 0) >= 5, ]
  
  errRates <- colSums(errorSub) / (colSums(errorSub) + colSums(pureSub))
  return(errRates)
  
  # errorRates <- as.matrix(errorSub / (errorSub + pureSub))
  # if(any(is.nan(errorRates))) {
  #   errorRates[is.nan(errorRates)] <- 0
  # }
  # 
  # purePC <- prcomp(t(pureSub), scale. = TRUE, rank. = 10)
  # pcErr <- abs(t(purePC$rotation)) %*% errorRates / colSums(abs(purePC$rotation))
  # return(pcErr)
}

plotErrRates.func <- function(errRates, tVec, remInd, title) {
  plotDat <- data.frame(
    errRates = errRates,
    t = tVec,
    filterVec = factor("Retained", levels = c("Retained", "Removed"))
  )
  plotDat$filterVec[remInd] <- "Removed"
  
  p <- ggplot(plotDat, aes(x = t, y = errRates, color = filterVec)) + 
    theme_classic() + 
    geom_jitter(size = 1, width = 0.2, height = 0) +
    labs(x = "day", y = "total error", title = title, color = "Filtering") + 
    ylim(c(0, max(ErrRates.func(hExpr.list$H100, mExpr.list$M0), 
                  ErrRates.func(mExpr.list$M100, hExpr.list$H0))))
  return(p)
}

## Results ##
H100.errPlot <- plotErrRates.func(ErrRates.func(hExpr.list$H100, mExpr.list$M0),
                                  hTime.list$H100, H100.remInd, "H100 empirical misalignment")
H0.errPlot <- plotErrRates.func(ErrRates.func(mExpr.list$M100, hExpr.list$H0),
                                hTime.list$H0, H0.remInd, "M100 empirical misalignment")


############################
## Sequencing depth plots ##
############################
## Functions ##
plotSeqDepth.func <- function(hDat, mDat, tVec, filterThresh, remInd, title) {
  plotDat <- data.frame(
    lExpr = log10(colSums(hDat) + colSums(mDat)),
    t = tVec,
    filterVec = factor("Retained", levels = c("Retained", "Removed"))
  )
  plotDat$filterVec[remInd] <- "Removed"
  
  p <- ggplot(plotDat, aes(x = t, y = lExpr, color = filterVec)) + 
    theme_classic() + 
    geom_jitter(size = 1, width = 0.2, height = 0) +
    labs(x = "day", y = "Total sequenced\nexpression (log 10)", title = title, color = "Filtering") + 
    geom_hline(yintercept = filterThresh, linetype = "dashed")
  return(p)
}

## Results ##
H100.expPlot <- plotSeqDepth.func(hExpr.list$H100, mExpr.list$M0, hTime.list$H100,
                                  H100.FT, H100.remInd, "H100 total expression")

H85.expPlot <- plotSeqDepth.func(hExpr.list$H85, mExpr.list$M15, hTime.list$H85,
                                 H85.FT, H85.remInd, "H85 total expression")

H10.expPlot <- plotSeqDepth.func(hExpr.list$H10, mExpr.list$M90, hTime.list$H10,
                                 H10.FT, H10.remInd, "H10 total expression")

H0.expPlot <- plotSeqDepth.func(hExpr.list$H0, mExpr.list$M100, hTime.list$H0,
                                H0.FT, H0.remInd, "M100 total expression")


################
## Final plot ##
################
jpeg(filename = "AlignmentQC/qcPlot.jpeg", quality = 0.9, height = 10.5, width = 8, res = 300, units = "in")
grid.arrange(H100.errPlot, H0.errPlot, H100.expPlot, H0.expPlot, H10.expPlot, H85.expPlot,
             ncol = 2)
dev.off()





