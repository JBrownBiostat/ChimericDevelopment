########################################
## Trendy Results: Pathway Enrichment ##
########################################
WorkDir <- "~/Documents/"


###############
## Libraries ##
###############
suppressMessages(library(piano))
suppressMessages(library(snowfall))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(gridExtra))
suppressMessages(library(grid))


##########
## Data ##
##########
setwd(WorkDir)
# load("PreprocessingAndNormalization/NormDat.RData")
load("Trendy/trendyFit.RData")
load("Trendy/GeneClasses.RData")

## Download gmt references for use ##
pathGSC <- loadGSC(file = paste0("Trendy/gmtReferences/c2.cp.v7.2.symbols.gmt"), type = "gmt")
tfGSC <- loadGSC(file = paste0("Trendy/gmtReferences/c3.all.v7.2.symbols.gmt"), type = "gmt")

themeObj <- theme(
  title = element_text(face = "bold", size = 8),
  axis.title = element_text(size = 7),
  axis.text = element_text(size = 6), 
  legend.title = element_text(size = 6),
  legend.text = element_text(size = 6)
)


###############
## Functions ##
###############
checkPeak <- function(trendy.ref, trendy.test){
  ret <- list(gene = c(), refBreak = c(), testBreak = c())
  for(i in names(trendy.ref)){
    ref <- trendy.ref[[i]]
    test <- trendy.test[[i]]
    if(!is.na(ref$Breakpoints[1]) & !is.na(test$Breakpoints[1])){
      breakRef <- 0; breakTest <- 0
      for(j in 1:length(ref$Breakpoints)) {
        if(ref$Segment.Trends[j] == 1 & ref$Segment.Trends[j + 1] %in% c(0, -1)){
          breakRef <- ref$Breakpoints[j]
          break
        }
        if(j == length(ref$Breakpoints)) {next}
      }
      for(j in 1:length(test$Breakpoints)) {
        if(test$Segment.Trends[j] == 1 & test$Segment.Trends[j + 1] %in% c(0, -1)){
          breakTest <- test$Breakpoints[j]
          break
        }
        if(j == length(test$Breakpoints)) {next}
      }
      if(breakRef != 0 & breakTest != 0) {
        l <- length(ret$gene)
        ret$gene[l+1] <- i
        ret$refBreak[l+1] <- breakRef
        ret$testBreak[l+1] <- breakTest
      }
    }
  }
  data.frame(
    gene = ret$gene,
    refBreak = ret$refBreak,
    testBreak = ret$testBreak
  )
}

checkUpTrend <- function(trendy.ref, trendy.test){
  ret <- list(gene = c(), refBreak = c(), testBreak = c(), 
              refSlope = c(), testSlope = c(),
              refSlopeAtTest = c(), testSlopeAtRef = c())
  for(i in names(trendy.ref)){
    ref <- trendy.ref[[i]]
    test <- trendy.test[[i]]
    if(any(ref$Segment.Trends == 1) & any(test$Segment.Trends == 1)){
      if(ref$Segment.Trends[1] != 0){
        upSeg <- which.max(ref$Segment.Trends == 1)
        if(any(ref$Segment.Trends[1:(upSeg-1)] == -1)){next}
      }
      if(test$Segment.Trends[1] != 0){
        upSeg <- which.max(test$Segment.Trends == 1)
        if(any(test$Segment.Trends[1:(upSeg-1)] == -1)){next}
      }
      l <- length(ret$gene)
      ret$gene[l+1] <- i
      ret$refBreak[l+1] <- c(0, ref$Breakpoint)[which.max(ref$Segment.Trends == 1)]
      ret$testBreak[l+1] <- c(0, test$Breakpoint)[which.max(test$Segment.Trends == 1)]
      ret$refSlope[l+1] <- ref$Segment.Slopes[which.max(ref$Segment.Trends == 1)]
      ret$testSlope[l+1] <- test$Segment.Slopes[which.max(test$Segment.Trends == 1)]
      ret$refSlopeAtTest[l+1] <- ifelse(is.na(ref$Breakpoints[1]),
                                        ref$Segment.Slopes[1],
                                        ref$Segment.Slopes[which.max(c(ref$Breakpoints, 42) >= ret$testBreak[l + 1])])
      ret$testSlopeAtRef[l+1] <- ifelse(is.na(test$Breakpoints[1]),
                                        test$Segment.Slopes[1],
                                        test$Segment.Slopes[which.max(c(test$Breakpoints, 42) >= ret$testBreak[l + 1])])
    }
  }
  data.frame(
    gene = ret$gene,
    refBreak = ret$refBreak,
    testBreak = ret$testBreak,
    refSlope = ret$refSlope,
    testSlope = ret$testSlope,
    refSlopeAtTest = ret$refSlopeAtTest,
    testSlopeAtRef = ret$testSlopeAtRef
  )
}

gsaTable <- function(gsaRes) {
  out <- data.frame(
    geneSet = names(gsaRes$gsc),
    dirStat = gsaRes$statDistinctDir,
    pAdjAcc = c(gsaRes$pAdjDistinctDirUp),
    pAdjDec = c(gsaRes$pAdjDistinctDirDn)
  )
  out <- out[order(out$pAdjAcc), ]
  return(out)
}

buildLevelStats <- function(trendy.ref, trendy.test) {
  peakOut <- checkPeak(trendy.ref, trendy.test)
  upOut <- checkUpTrend(trendy.ref, trendy.test)
  lsMat <- data.frame(
    Gene = names(trendy.test),
    dDay = 0
  )
  matchVec <- match(upOut$gene, lsMat$Gene)
  lsMat$dDay[matchVec] <- upOut$refBreak - upOut$testBreak
  matchVec <- match(peakOut$gene, lsMat$Gene)
  lsMat$dDay[matchVec] <- peakOut$refBreak - peakOut$testBreak
  
  lsMat <- matrix(c(lsMat$dDay), ncol = 1)
  rownames(lsMat) <- names(trendy.test)
  
  return(lsMat)
}

calculateEnrichment <- function(LSVec, pathGSC, tfGSC, nPerm, ncpus, minSize = 1, maxSize = 250) {
  pathOut <- runGSA(LSVec, 
                    gsc = pathGSC, 
                    nPerm = ceiling(nPerm / ncpus) * ncpus,
                    ncpus = ncpus, 
                    gsSizeLim = c(minSize, maxSize))
  pathTable <- gsaTable(pathOut)
  
  tfOut <- runGSA(LSVec, 
                  gsc = tfGSC, 
                  nPerm = ceiling(nPerm / ncpus) * ncpus,
                  ncpus = ncpus, 
                  gsSizeLim = c(minSize, maxSize))
  tfTable <- gsaTable(tfOut)
  
  return(list(
    pathEnrich = pathTable,
    tfEnrich = tfTable
  ))
}

plotGO <- function(GOdat, plotTitle = "Top terms"){
  GOdat <- cleanGOdat(GOdat)
  GOdat$pAdjAcc <- log10(GOdat$pAdjAcc)
  if(any(is.na(GOdat$pAdjAcc))) {
    GOdat$pAdjAcc[is.na(GOdat$pAdjAcc)] <- min(GOdat$pAdjAcc[!is.na(GOdat$pAdjAcc)])
  }
  xRange <- max(GOdat$pAdjAcc) - min(GOdat$pAdjAcc)
  p <- ggplot(GOdat, aes(x = pAdjAcc, y = forcats::fct_rev(geneSet))) +
    theme_linedraw() +
    geom_point() +
    xlim(min(GOdat$pAdjAcc) - 0.1 * xRange,
         max(GOdat$pAdjAcc) + 0.1 * xRange) + 
    themeObj + 
    theme(axis.ticks.y = element_blank(),
          strip.text.y = element_text(size=6,face="bold",angle=0),
          plot.background = element_rect(fill="white", colour="white"),
          axis.text=element_text(size=6))+
    labs(x="FDR (log 10)",y="", title = plotTitle)
  return(p)
}

cleanGOdat <- function(GOdat) {
  GOdat <- GOdat[1:10, ]
  GOdat$geneSet <- gsub("...", "", GOdat$geneSet, fixed = TRUE)
  GOdat$geneSet <- unlist(lapply(strsplit(GOdat$geneSet, "_"), function(x) {
    n <- cumsum(nchar(x))
    n <- n + seq_len(length(n)) - 1
    if(length(x) == 1){return(x)}
    ret <- x[1]
    m <- 0
    for(i in 2:length(x)) {
      if(n[i - 1] >= m + 16) {
        m <- n[i - 1]
        ret <- c(ret, "\n", x[i])
      } else {
        ret <- c(ret, " ", x[i])
      }
    }
    return(paste(ret, collapse = ""))
  }))
  GOdat$geneSet <- factor(as.character(GOdat$geneSet), 
                       levels = as.character(GOdat$geneSet), 
                       ordered = TRUE)
  return(GOdat)
}


##############
## Analysis ##
##############
H10.LS <- buildLevelStats(h100.trendy, h10.trendy)
H10.enrich <- calculateEnrichment(H10.LS, pathGSC, tfGSC, 4e6, 8)
head(H10.enrich$pathEnrich, 10)
head(H10.enrich$tfEnrich, 10)

H85.LS <- buildLevelStats(h100.trendy, h85.trendy)
H85.enrich <- calculateEnrichment(H85.LS, pathGSC, tfGSC, 4e6, 8)
head(H85.enrich$pathEnrich, 10)
head(H85.enrich$tfEnrich, 10)

hGenes <- names(h100.trendy)
mGenes <- names(m100.trendy)
mGenes <- mGenes[toupper(mGenes) %in% hGenes]
mGenes <- mGenes[!duplicated(toupper(mGenes))]
hGenes <- hGenes[hGenes %in% toupper(mGenes)]
all(hGenes == toupper(mGenes))
hTest.trendy <- lapply(hGenes, function(hGene) {h100.trendy[[hGene]]})
mTest.trendy <- lapply(mGenes, function(mGene) {m100.trendy[[mGene]]})
names(hTest.trendy) <- hGenes
names(mTest.trendy) <- hGenes
MTest.LS <- buildLevelStats(hTest.trendy, mTest.trendy)
MTest.enrich <- calculateEnrichment(MTest.LS, pathGSC, tfGSC, 4e6, 8)
head(MTest.enrich$pathEnrich, 10)
head(MTest.enrich$tfEnrich, 10)

head(MTest.enrich$tfEnrich[MTest.enrich$tfEnrich$geneSet %in% 
                             H10.enrich$tfEnrich$geneSet[H10.enrich$tfEnrich$pAdjAcc <= 1e-2], ], 10)
head(H10.enrich$tfEnrich, 10)

save(H10.enrich, H85.enrich, MTest.enrich,
     file = "Trendy/pathEnrich.RData",
     compress = "xz")

################
## Write .csv ##
################
system("mkdir Trendy/pathEnrichDocs")
write.csv(H10.enrich$pathEnrich, "Trendy/pathEnrichDocs/H10_PathEnrich.csv")
write.csv(H10.enrich$tfEnrich, "Trendy/pathEnrichDocs/H10_TFEnrich.csv")

write.csv(H85.enrich$pathEnrich, "Trendy/pathEnrichDocs/H85_PathEnrich.csv")
write.csv(H85.enrich$tfEnrich, "Trendy/pathEnrichDocs/H85_TFEnrich.csv")

write.csv(MTest.enrich$pathEnrich, "Trendy/pathEnrichDocs/M100_PathEnrich.csv")
write.csv(MTest.enrich$tfEnrich, "Trendy/pathEnrichDocs/M100_TFEnrich.csv")

write.csv(MTest.enrich$pathEnrich[
  MTest.enrich$pathEnrich$geneSet %in% 
    H10.enrich$pathEnrich$geneSet[H10.enrich$pathEnrich$pAdjAcc <= 1e-2], 
], "Trendy/pathEnrichDocs/M100_Subset_PathEnrich.csv")
write.csv(MTest.enrich$tfEnrich[
  MTest.enrich$tfEnrich$geneSet %in% 
    H10.enrich$tfEnrich$geneSet[H10.enrich$tfEnrich$pAdjAcc <= 1e-2], 
], "Trendy/pathEnrichDocs/M100_Subset_TFEnrich.csv")


##########
## Plot ##
##########
p.h10Path <- plotGO(H10.enrich$pathEnrich, "H10 Pathways")
p.h10TF <- plotGO(H10.enrich$tfEnrich, "H10 TF/miRNA")
p.h10 <- grid.arrange(p.h10Path, p.h10TF, nrow = 1)
p.h10 <- annotate_figure(p.h10, fig.lab = "A", fig.lab.size = 14, fig.lab.face = "bold")

MTest.enrich$tfEnrich <- MTest.enrich$tfEnrich[MTest.enrich$tfEnrich$geneSet %in% 
                                                 H10.enrich$tfEnrich$geneSet[H10.enrich$tfEnrich$pAdjAcc <= 1e-2], ]
MTest.enrich$pathEnrich <- MTest.enrich$pathEnrich[MTest.enrich$pathEnrich$geneSet %in% 
                                                     H10.enrich$pathEnrich$geneSet[H10.enrich$pathEnrich$pAdjAcc <= 1e-2], ]
p.hSubPath <- plotGO(MTest.enrich$pathEnrich, "H10/M100 Subset Pathways")
p.hSubTF <- plotGO(MTest.enrich$tfEnrich, "H10/M100 Subset TF/miRNA")
p.hSub <- grid.arrange(p.hSubPath, p.hSubTF, nrow = 1)
p.hSub <- annotate_figure(p.hSub, fig.lab = "B", fig.lab.size = 14, fig.lab.face = "bold")

jpeg("Trendy/TrendyFigs/PathEnrich.jpeg", width = 8, height = 8, 
     units = "in", res = 300, quality = 0.95)
grid.arrange(p.h10, p.hSub, ncol = 1)
dev.off()




