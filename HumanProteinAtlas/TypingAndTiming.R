###############################
## In-Vivo typing and timing ##
###############################
WorkDir <- "~/Documents/"
setwd(paste0(WorkDir, "/HumanProteinAtlas/"))


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
suppressMessages(library(org.Hs.eg.db))


##########
## Data ##
##########
HPA <- read.table("RawData/rna_brain_fantom.tsv", sep = "\t", header = TRUE)
HPA.mat <- do.call(cbind, lapply(unique(HPA$Brain.region), function(BR) {
  m <- HPA[HPA$Brain.region == BR, 6, drop = FALSE]
  rownames(m) <- HPA$Gene[HPA$Brain.region == BR]
  colnames(m) <- BR
  return(m)
}))

load("../PreprocessingAndNormalization/NormDat.RData")
SYMB <- rownames(dat.h100$allDay)
ENS <- mapIds(org.Hs.eg.db, keys = SYMB, keytype = "SYMBOL", column = "ENSEMBL")
SYMB <- SYMB[!duplicated(ENS) & !is.na(ENS)]
ENS <- ENS[!duplicated(ENS) & !is.na(ENS)]
SYMB <- SYMB[ENS %in% rownames(HPA.mat)]
ENS <- ENS[ENS %in% rownames(HPA.mat)]

HPA.mat <- HPA.mat[rownames(HPA.mat) %in% ENS, ]
subTypeInd <- colnames(HPA.mat) %in% 
  c("amygdala", "basal ganglia", "cerebellum", 
    "cerebral cortex", "hippocampal formation", "olfactory region")

SYMB <- SYMB[match(ENS, rownames(HPA.mat))]
ENS <- ENS[match(ENS, rownames(HPA.mat))]

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
names(maxCV) <- names(h100.trendy)
maxCV <- maxCV[intersect(names(h100.trendy), SYMB)]

HPA.cv <- apply(HPA.mat[ENS, subTypeInd], 1, function(x) {sd(x) / mean(x)})

topGenes <- unique(c(
  names(maxCV)[order(maxCV, decreasing = T)[1:1.5e3]],
  SYMB[order(HPA.cv, decreasing = T)[1:1.5e3]]
))
topGenes <- intersect(topGenes, names(maxCV))
topGenes <- intersect(topGenes, SYMB)


#########################
## Dimension reduction ##
#########################
# Input data #
exptDat <- cbind(dat.h100$sumDayEC, dat.h85$sumDayEC, dat.h10$sumDayEC)
IVDat <- as.matrix(HPA.mat[, subTypeInd])
geneInd <- rowSums(cbind(exptDat[SYMB, ], IVDat[ENS, ])) >= 5
ivSF <- calculateSumFactors(cbind(exptDat[SYMB[geneInd], ], IVDat[ENS[geneInd], ]))
h100NC <- ncol(dat.h100$sumDayEC)
h85NC <- ncol(dat.h85$sumDayEC)
h10NC <- ncol(dat.h10$sumDayEC)

# Dimension reduction #
pcaOut <- glmpca(Y = cbind(exptDat[topGenes, ], IVDat[ENS[SYMB %in% topGenes], ]), 
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
       title = "PCA distances", 
       color = "Species\nmixture") + 
 themeObj
p.dist

save(p.dist,
     file = "DistanceFig.RData", 
     compress = "xz")




