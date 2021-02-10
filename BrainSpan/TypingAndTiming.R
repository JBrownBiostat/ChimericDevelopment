###############################
## In-Vivo typing and timing ##
###############################
WorkDir <- "~/Documents/"
setwd(paste0(WorkDir, "/BrainSpan/"))


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
suppressMessages(library(irlba))
suppressMessages(library(scran))
suppressMessages(library(MuSiC))
suppressMessages(library(Matrix))
suppressMessages(library(Biobase))
suppressMessages(library(xbioc))
suppressMessages(library(DirichletReg))
suppressMessages(library(splines))
suppressMessages(library(glmpca))


##########
## Data ##
##########
load("../PreprocessingAndNormalization/NormDat.RData")
load("../Trendy/GeneClasses.RData")
load("../Trendy/trendyFit.RData")
load("../Trendy/EnrichmentResults.RData")

themeObj <- theme(
  title = element_text(face = "bold", size = 8),
  axis.title = element_text(size = 7),
  axis.text = element_text(size = 6), 
  legend.title = element_text(size = 6),
  legend.text = element_text(size = 6)
)

exptCol <- c(H10 = "deeppink", H85 = "darkorchid", H100 = "cornflowerblue", M100 = "darkolivegreen2")

# BS: All Genes #
AG <- read.csv("RawData/AllGenes/Expression.csv",header=FALSE,row.names=1)
AG.col <- read.csv("RawData/AllGenes/Columns.csv")
AG.row <- read.csv("RawData/AllGenes/Rows.csv")
AG.age <- strsplit(as.character(AG.col$age)," ")
AG.age <- data.frame(weeks=sapply(AG.age,function(x)x[1]),typ=sapply(AG.age,function(x)x[2]))
AG.dat <- matrix(NA,nrow(AG),length(unique(AG.col$structure_acronym))*6); colnames(AG.dat) <- paste0("X",1:ncol(AG.dat))
age <- unique(AG.col$age)
str <- unique(AG.col$structure_acronym); for(i in 1:6){str <- intersect(str,AG.col$structure_acronym[which(AG.col$age==age[i])])}
str <- sort(str)
index <- 1
for(i in 1:6){
  for(j in 1:length(str)){
    newDat <- AG[,which(AG.col$age==age[i]&AG.col$structure_acronym==str[j])]
    if(!is.null(nrow(newDat))){newDat <- rowMeans(newDat)}
    AG.dat[,index] <- newDat
    colnames(AG.dat)[index] <- paste0(str[j],".",strsplit(as.character(age[i])," ")[[1]][1])
    index <- index+1
  }
}
AG.dat <- AG.dat[,c(5*(1:3)-4,5*(1:3)-3,5*(1:3)-2,5*(1:3)-1,5*(1:3))]
rownames(AG.dat) <- AG.row$gene_symbol
AG.dat <- AG.dat[which(rownames(AG.dat)%in%rownames(dat.h100$sumDay)),]
AG.dat <- AG.dat[-which(duplicated(rownames(AG.dat))),]
AG.dat <- AG.dat[order(rownames(AG.dat)),]
week9Col <- c("AMY.9", "DFC.9", "HIP.9", "MFC.9", "OFC.9")


#########################
## Dimension reduction ##
#########################
# Input data #
exptDat <- cbind(dat.h100$sumDayFpkm, dat.h85$sumDayFpkm, dat.h10$sumDayFpkm)
IVDat <- AG.dat[, week9Col]
geneVec <- intersect(rownames(exptDat), rownames(IVDat))
geneVec <- geneVec[rowSums(cbind(exptDat[geneVec, ], IVDat[geneVec, ])) >= 5]
ivSF <- calculateSumFactors(cbind(exptDat[geneVec, ], IVDat[geneVec, ]))
h100NC <- ncol(dat.h100$sumDayEC)
h85NC <- ncol(dat.h85$sumDayEC)
h10NC <- ncol(dat.h10$sumDayEC)

# Variable genes #
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

IV.cv <- apply(IVDat, 1, function(x) {sd(x) / mean(x)})

geneVec <- unique(c(
  names(maxCV)[order(maxCV, decreasing = T)[1:1.5e3]],
  names(IV.cv)[order(IV.cv, decreasing = T)[1:1.5e3]]
))
geneVec <- intersect(geneVec, names(maxCV))
geneVec <- intersect(geneVec, names(IV.cv))

# Dimension reduction #
pcaOut <- glmpca(Y = cbind(exptDat[geneVec, ], IVDat[geneVec, ]), 
                 L = 6, 
                 fam = "nb", 
                 optimizer = "fisher",
                 ctl = list(penalty = 10, minIter = 400),
                 sz = ivSF)


####################
## Plot distances ##
####################
## Calculate distances ##
dist.func <- function(hRed, IVRed, hDat, IVDat, expt) {
  do.call(rbind, lapply(1:ncol(IVDat), function(i) {
    data.frame(
      expt = expt,
      IVref = colnames(IVDat)[i],
      eDist = sqrt(colSums((t(hRed) - t(IVRed)[, i]) ^ 2)),
      sampT = sort(unique(hDat$days))
    )
  }))
}

H100.red <- pcaOut$factors[1:h100NC, ]
H85.red <- pcaOut$factors[(h100NC + 1):(h100NC + h85NC), ]
H10.red <- pcaOut$factors[(h100NC + h85NC + 1):(h100NC + h85NC + h10NC), ]
IV.red <- pcaOut$factors[(h100NC + h85NC + h10NC + 1):nrow(pcaOut$factors), ]

h100.dist <- dist.func(H100.red, IV.red, dat.h100, IVDat, "H100")
h85.dist <- dist.func(H85.red, IV.red, dat.h85, IVDat, "H85")
h10.dist <- dist.func(H10.red, IV.red, dat.h10, IVDat, "H10")

plot(h100.dist$sampT, h100.dist$eDist, col = 1, pch = 19, cex = 0.5)
points(h85.dist$sampT, h85.dist$eDist, col = 2, pch = 19, cex = 0.5)
points(h10.dist$sampT, h10.dist$eDist, col = 3, pch = 19, cex = 0.5)

## Plot distances ##
plotDat <- rbind(h100.dist, h85.dist, h10.dist)
plotDat$expt <- factor(plotDat$expt, ordered = T, levels = c("H100", "H85", "H10"))

p.dist <- ggplot(plotDat, aes(x = sampT, y = eDist, color = expt)) + 
  theme_classic() + 
  geom_point(alpha = 0.5) + 
  geom_smooth(se = T) +
  scale_color_manual(values = exptCol) + 
  labs(x = "Day", y = "Distance", 
       title = "Week 9 PCA", 
       color = "Species\nmixture") + 
 themeObj
p.dist

save(p.dist,
     file = "DistanceFig.RData", 
     compress = "xz")




