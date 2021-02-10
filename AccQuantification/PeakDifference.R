############################################
## Acceleration from differences in peaks ##
############################################
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
suppressMessages(library(splines))


##########
## Data ##
##########
load("PreprocessingAndNormalization/NormDat.RData")
load("Trendy/GeneClasses.RData")
load("Trendy/trendyFit.RData")
load("Trendy/EnrichmentResults.RData")

exptCol <- c(H10 = "deeppink", H85 = "darkorchid", H100 = "cornflowerblue", M100 = "darkolivegreen2")
themeObj <- theme(
  title = element_text(face = "bold", size = 10),
  axis.title = element_text(size = 9),
  axis.text = element_text(size = 8)
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
  return(ret)
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
  return(ret)
}


##############
## Analysis ##
##############
peak85Dat <- checkPeak(h100.trendy, h85.trendy)
peak10Dat <- checkPeak(h100.trendy, h10.trendy)
genes <- intersect(peak10Dat$gene, peak85Dat$gene)
peakDat <- data.frame(
  gene = genes,
  H100Peak = peak85Dat$refBreak[match(genes, peak85Dat$gene)],
  H85Peak = peak85Dat$testBreak[match(genes, peak85Dat$gene)],
  H10Peak = peak10Dat$testBreak[match(genes, peak10Dat$gene)]
)

up85Dat <- checkUpTrend(h100.trendy, h85.trendy)
up10Dat <- checkUpTrend(h100.trendy, h10.trendy)
genes <- intersect(up10Dat$gene, up85Dat$gene)
upDat <- data.frame(
  gene = genes,
  H100Peak = up85Dat$refBreak[match(genes, up85Dat$gene)],
  H85Peak = up85Dat$testBreak[match(genes, up85Dat$gene)],
  H10Peak = up10Dat$testBreak[match(genes, up10Dat$gene)]
)
# upDat <- upDat[upDat$H10Peak > 0 & upDat$H85Peak > 0 & upDat$H100Peak > 0, ]

plotDat <- rbind(peakDat, upDat)

M100_Genes <- toupper(names(m100.trendy))
subGenes <- intersect(M100_Genes, names(h100.trendy))
subM100.trendy <- m100.trendy
names(subM100.trendy) <- M100_Genes
subM100.trendy <- subM100.trendy[match(subGenes, M100_Genes)]
subH100.trendy <- h100.trendy[match(subGenes, names(h100.trendy))]

peakM100Dat <- checkPeak(subH100.trendy, subM100.trendy)
mPeakDat <- data.frame(
  gene = peakM100Dat$gene,
  H100Peak = peakM100Dat$refBreak,
  M100Peak = peakM100Dat$testBreak, 
  M100FC = peakM100Dat$refBreak / peakM100Dat$testBreak
)

upM100Dat <- checkUpTrend(subH100.trendy, subM100.trendy)
mUpDat <- data.frame(
  gene = upM100Dat$gene,
  H100Peak = upM100Dat$refBreak,
  M100Peak = upM100Dat$testBreak, 
  M100FC = upM100Dat$refBreak / upM100Dat$testBreak
)
mUpDat <- mUpDat[mUpDat$H100Peak > 0 & mUpDat$M100Peak > 0, ]

mDat <- rbind(mPeakDat, mUpDat)


##########
## Plot ##
##########
p85 <- ggplot(plotDat) + 
  theme_classic() + 
  geom_abline(slope = 1, intercept = 0, color = "black") + 
  geom_point(aes(x = H85Peak, y = H100Peak), color= "darkorchid1", size = 0.25, alpha = 0.5) + 
  geom_smooth(aes(x = H85Peak, y = H100Peak), color = exptCol[2], 
              method = "lm", formula = y ~ bs(x), lwd = 1) + 
  labs(x = "H85", y = "H100", title = "H85 Acceleration\nrelative to H100") + 
  themeObj
p85

p10 <- ggplot(plotDat) + 
  theme_classic() + 
  geom_abline(slope = 1, intercept = 0, color = "black") + 
  geom_point(aes(x = H10Peak, y = H100Peak), color= "deeppink1", size = 0.25, alpha = 0.5) + 
  geom_smooth(aes(x = H10Peak, y = H100Peak), color = exptCol[1], 
              method = "lm", formula = y ~ bs(x), lwd = 1) + 
  labs(x = "H10", y = "H100", title = "H10 Acceleration\nrelative to H100") + 
  themeObj
p10

p85 <- ggplot(plotDat) + 
  theme_classic() + 
  geom_hline(yintercept = 1, color = "black") + 
  geom_point(aes(x = H85Peak, y = H85FC), color= "darkorchid1", size = 0.25, alpha = 0.5) + 
  geom_smooth(aes(x = H85Peak, y = H85FC), color = exptCol[2], method = "lm", formula = y ~ bs(x), lwd = 1) + 
  labs(x = "H85", y = "Acc. factor", title = "H85 Acceleration\nrelative to H100") + 
  ylim(0, 7.5) + xlim(c(2, 40)) + 
  themeObj

p10 <- ggplot(plotDat) + 
  theme_classic() + 
  geom_hline(yintercept = 1, color = "black") + 
  geom_point(aes(x = H10Peak, y = H10FC), color= "deeppink1", size = 0.25, alpha = 0.5) + 
  geom_smooth(aes(x = H10Peak, y = H10FC), color = exptCol[1], method = "lm", formula = y ~ bs(x), lwd = 1) + 
  labs(x = "H10", y = "Acc. factor", title = "H10 Acceleration\nrelative to H100") + 
  ylim(0, 7.5) + xlim(c(2, 40)) + 
  themeObj

p1085 <- ggplot(plotDat) + 
  theme_classic() + 
  geom_hline(yintercept = 1, color = "black") + 
  geom_point(aes(x = H10Peak, y = H10H85FC), color= "chocolate1", size = 0.25, alpha = 0.5) + 
  geom_smooth(aes(x = H10Peak, y = H10H85FC), color = "chocolate3", method = "lm", formula = y ~ bs(x), lwd = 1) + 
  labs(x = "H10", y = "Acc. factor", title = "H10 Acceleration\nrelative to H85") + 
  ylim(0, 7.5) + xlim(c(2, 40)) + 
  themeObj

pM <- ggplot(mDat) + 
  theme_classic() + 
  geom_hline(yintercept = 1, color = "black") + 
  geom_point(aes(x = M100Peak, y = M100FC), color= "darkolivegreen1", size = 0.25, alpha = 0.5) + 
  geom_smooth(aes(x = M100Peak, y = M100FC), color = "darkolivegreen2", method = "lm", formula = y ~ bs(x), lwd = 1) + 
  labs(x = "H10", y = "Acc. factor", title = "M100 Acceleration\nrelative to H100") + 
  ylim(0, 7.5) + xlim(c(2, 40)) + 
  themeObj

jpeg("AccQuantification/RelPeakUp.jpeg", 
     height = 3, width = 8, units = "in", res = 350, quality = 0.95)
grid.arrange(p10, p85, p1085, pM, nrow = 1)
dev.off()




