#################################
## Basic autocorrelation plots ##
#################################
WorkDir <- "~/Documents/"


###############
## Libraries ##
###############
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(scales))
suppressMessages(library(reshape2))
suppressMessages(library(viridis))
suppressMessages(library(topGO))
suppressMessages(library(BiocParallel))


##########
## Data ##
##########
setwd(WorkDir)
load("PreprocessingAndNormalization/NormDat.RData")
load("Trendy/GeneClasses.RData")
load("Trendy/trendyFit.RData")
load("Trendy/EnrichmentResults.RData")

if (.Platform$OS.type == "windows") {
  prll <- SnowParam(workers = 8, progressbar = F)
  register(BPPARAM = prll, default = TRUE)
} else {
  prll <- MulticoreParam(workers = 8, progressbar = F)
  register(BPPARAM = prll, default = TRUE)
}

themeObj <- theme(
  title = element_text(face = "bold", size = 8),
  axis.title = element_text(size = 7),
  axis.text = element_text(size = 6), 
  legend.title = element_text(size = 6),
  legend.text = element_text(size = 6)
)


#########################
## Find variable genes ##
#########################
H100.cv <- sapply(h100.trendy, function(subTrend) {
  sd(subTrend$Fitted.Values) / mean(subTrend$Fitted.Values)
})
H85.cv <- sapply(h85.trendy, function(subTrend) {
  sd(subTrend$Fitted.Values) / mean(subTrend$Fitted.Values)
})
H10.cv <- sapply(h10.trendy, function(subTrend) {
  sd(subTrend$Fitted.Values) / mean(subTrend$Fitted.Values)
})
maxCV <- apply(cbind(H100.cv, H85.cv, H10.cv), 1, max)

topGenes <- names(h100.trendy)[order(maxCV, decreasing = T)[1:2000]]


#######################
## Correlation plots ##
#######################
lossF <- function(dat, w, lCoef) {
  d <- abs(dat[, 1] * lCoef[1] + dat[, 2] * lCoef[2] + lCoef[3]) / 
    sqrt(sum(lCoef[1:2] ^ 2))
  projX <- (lCoef[2] * (lCoef[2] * dat[, 1] - lCoef[1] * dat[, 2]) - lCoef[1] * lCoef[3]) / 
    sum(lCoef[1:2] ^ 2)
  projY <- (lCoef[1] * (-lCoef[2] * dat[, 1] + lCoef[1] * dat[, 2]) - lCoef[2] * lCoef[3]) / 
    sum(lCoef[1:2] ^ 2)
  return(cbind(d, projX, projY))
}

f.func1 <- function(slopeVec, dat, w) {
  sum(lossF(dat, w, c(slopeVec, 0)))
}

f.func2 <- function(coefVec, BP, dat, w) {
  coef1 <- c(-coefVec[1], 1, 0)
  coef2 <- c(-coefVec[2], 1, -(coefVec[1] * BP - coefVec[2] * BP))
  loss1 <- lossF(dat, w, coef1)
  loss2 <- lossF(dat, w, coef2)
  loss3 <- sqrt((BP - dat[, 1]) ^ 2 + 
                      (coefVec[1] * BP - dat[, 2]) ^ 2)
  
  BPY <- BP * coefVec[1]
  
  d <- ifelse((loss1[, 2] <= BP & loss1[, 3] <= BPY) & (loss2[, 2] >= BP & loss2[, 3] >= BPY),
              pmin(loss1[, 1], loss2[, 1]), 
              ifelse((loss1[, 2] <= BP & loss1[, 3] <= BPY), 
                     loss1[, 1], 
                     ifelse((loss2[, 2] >= BP & loss2[, 3] >= BPY), 
                            loss2[, 1], 
                            loss3)))
  
  sum((w ^ 2 * d ^ 2))
}

bootOptim.func <- function(optimDat, lRef, lTest, nBoot) {
  w <- buildCorVec.func(lRef, lTest)$value
  optimOut <- optim(c(1, 1), f.func2, BP = 16, dat = optimDat, w = w,
                    lower = 1 / 10, upper = 10, method = "L-BFGS-B")
  optimPar <- optimOut$par
  bootSlope <- bplapply(1:nBoot, function(i, optimDat, lRef, lTest) {
    w <- buildCorVec.func(lRef, lTest, boot = T)
    optimOut <- optim(c(1, 1), f.func2, 
                      BP = 16, dat = optimDat, w = w,
                      lower = 1 / 10, upper = 10, method = "L-BFGS-B")
    return(optimOut$par)
  }, optimDat = optimDat, lRef = lRef, lTest = lTest, BPPARAM = prll)
  bootSlope <- do.call(rbind, bootSlope)
  return(list(
    slope1 = optimPar[1],
    slope1SE = sd(bootSlope[, 1]),
    slope2 = optimPar[2],
    slope2SE = sd(bootSlope[, 2]),
    acc1Rate = 1 / optimPar[1],
    acc1SE = sd(1 / bootSlope[, 1]),
    acc2Rate = 1 / optimPar[2],
    acc2SE = sd(1 / bootSlope[, 2])
  ))
}

buildCorVec.func <- function(lRef, lTest, boot = F) {
  if(boot) {
    rowInd <- sample.int(nrow(lRef), nrow(lRef), replace = T)
    corMat <- cor(lRef[rowInd, ], lTest[rowInd, ], method = "s")
  } else {
    corMat <- cor(lRef, lTest, method = "s")
  }
  
  if(any(corMat < 0)) {
    corMat[corMat < 0] <- 0
  }
  
  if(boot) {
    return(c(corMat))
  } else {
    return(melt(corMat))
  }
}

corPlot.func <- function(testDat, refDat, saveFile, plotTitle, 
                         xtitle, ytitle, corLim = NULL, nBoot){
  lRef <- log1p(refDat); lTest <- log1p(testDat)
  w <- buildCorVec.func(lRef, lTest)
  plotDat <- cbind(w,
                   x.cent=rep(c(0:7,8.25,2*(5:21)),26),
                   x.width=rep(c(rep(1,8),1.5,rep(2,17)),26),
                   y.cent=rep(c(0:7,8.25,2*(5:21)),each=26),
                   y.width=rep(c(rep(1,8),1.5,rep(2,17)),each=26))
  optimDat <- plotDat[, c("x.cent", "y.cent")]
  bootOut <- bootOptim.func(optimDat, lRef, lTest, nBoot)
  
  plotBoot <- data.frame(
    x = c(0, seq(16, 42, length = 1e2)),
    y = c(0, seq(bootOut$slope1 * 16, 
                 bootOut$slope2 * 42 - bootOut$slope2 * 16 + bootOut$slope1 * 16,
                 length = 1e2))
  )
  colnames(plotDat)[1:3] <- c("testInd","refInd","sCor")
  if(is.null(corLim)){corLim=quantile(plotDat$sCor,probs=c(0,0.5,1),na.rm=TRUE)}
  plotDat$plotCor <- pmax(corLim[1],pmin(corLim[3],plotDat$sCor))
  p <- ggplot(plotDat,aes(x=x.cent,y=y.cent,width=x.width,height=y.width,fill=sCor))+
    theme_classic() + 
    geom_tile()+
    scale_fill_gradientn(colors=c("darkblue","white","darkred")) +
    geom_abline(slope = 1, intercept = 0, col = "white") +
    geom_line(data = plotBoot, aes(x = x, y = y), col = "black", inherit.aes = F, lwd = 1) +
    labs(x=xtitle, y=ytitle, fill="S-cor",
         title = plotTitle,
         subtitle = paste0("Acc1 = ", round(bootOut$acc1Rate, 3), " (",
                           round(bootOut$acc1SE, 3), ")\nAcc2 = ",
                           round(bootOut$acc2Rate, 3), " (", 
                           round(bootOut$acc2SE, 3), ")"))+
    scale_y_continuous(breaks = 4 * (0:10), limits = c(0, 41)) +
    scale_x_continuous(breaks = 4 * (0:10), limits = c(0, 41)) +
    themeObj
  jpeg(saveFile, width = 4, height = 4, units = "in", res = 350, quality = 0.95)
  plot(p)
  dev.off()
  return(list(plotDat = plotDat, plot = p))
}

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

##########
nBoot <- 1e4

h100.null <- corPlot.func(dat.h100$sumDay[topGenes, ], dat.h100$sumDay[topGenes, ], 
                          "AutoCorrelation/BasePlots/h100null.jpeg",
                          "H100 reference", "H100", "H100", c(0.0, 0.845, 0.99), nBoot)
h85.null <- corPlot.func(dat.h85$sumDay[topGenes, ], dat.h85$sumDay[topGenes, ], 
                         "AutoCorrelation/BasePlots/h85null.jpeg",
                         "H85 reference", "H85","H85", c(0.0, 0.845, 0.99), nBoot)
h10.null <- corPlot.func(dat.h10$sumDay[topGenes, ], dat.h10$sumDay[topGenes, ], 
                         "AutoCorrelation/BasePlots/h10null.jpeg",
                         "H10 reference", "H10", "H10", c(0.0, 0.845, 0.99), nBoot)

h100_10 <- corPlot.func(dat.h10$sumDay[topGenes, ], dat.h100$sumDay[topGenes, ], 
                        "AutoCorrelation/BasePlots/h100_10.jpeg",
                        "H10 by H100", "H100", "H10", c(0.0, 0.845, 0.99), nBoot)
h100_85 <- corPlot.func(dat.h85$sumDay[topGenes, ], dat.h100$sumDay[topGenes, ],
                        "AutoCorrelation/BasePlots/h100_85.jpeg",
                        "H85 by H100","H100","H85",c(0.0,0.845,0.99), nBoot)
h85_10 <- corPlot.func(dat.h10$sumDay[topGenes, ], dat.h85$sumDay[topGenes, ], 
                       "AutoCorrelation/BasePlots/h85_10.jpeg",
                       "H10 by H85", "H85", "H10", c(0.0, 0.845, 0.99), nBoot)

HGenes <- rownames(dat.h100$allDay)[rownames(dat.h100$allDay) %in% toupper(rownames(dat.m100$allDay))]
HGenes <- intersect(HGenes, topGenes)
MGenes <- rownames(dat.m100$allDay)[toupper(rownames(dat.m100$allDay)) %in% topGenes]
MGenes <- MGenes[!duplicated(toupper(MGenes))]

h100_0 <- corPlot.func(dat.m100$sumDay[MGenes, ], dat.h100$sumDay[HGenes, ], 
                       "AutoCorrelation/BasePlots/h100_0.jpeg",
                       "M100 by H100","H100","M100",c(0.0,0.845,0.99), nBoot)

LP <- get_legend(h100_10$plot)

save(h100_10, h100_85, h85_10, h100_0, h100.null, h85.null, h10.null,
     file = "AutoCorrelation/BasePlots/Plots.RData",
     compress = "xz")
jpeg("AutoCorrelation/BasePlots/H10Triplicate.jpeg", width = 12, height = 4, 
     units = "in", res = 350, quality = 0.95)
grid.arrange(h100_10$plot + theme(legend.position = "none"), 
             h100_85$plot + theme(legend.position = "none"), 
             h10.null$plot + theme(legend.position = "none"), 
             LP,
             ncol = 4, widths = c(1, 1, 1, 0.2))
dev.off()

jpeg("AutoCorrelation/BasePlots/Fig1Triplicate.jpeg", width = 12, height = 4, 
     units = "in", res = 350, quality = 0.95)
grid.arrange(h100_0$plot + theme(legend.position = "none"), 
             h100_10$plot + theme(legend.position = "none"), 
             h10.null$plot + theme(legend.position = "none"), 
             LP,
             ncol = 4, widths = c(1, 1, 1, 0.2))
dev.off()




