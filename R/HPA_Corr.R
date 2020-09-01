###########################################################
## Human Protein Atlas Correlation with observed results ##
###########################################################

setwd("~/Documents/Research/Chris_Mixing_2019/Final_MixData_March25/HumanProteinAtlas/")


###############
## Libraries ##
###############
library(ggplot2)
library(gridExtra)
library(viridis)
library("org.Hs.eg.db")
library(topGO)


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

SYMB <- SYMB[match(ENS, rownames(HPA.mat))]
ENS <- ENS[match(ENS, rownames(HPA.mat))]

load("../Trendy/GeneClasses.RData")
load("../Trendy/trendyFit.RData")
load("../Trendy/EnrichmentResults.RData")


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

topGenes <- names(maxCV)[order(maxCV, decreasing = T)[1:1000]]


############################
## Calculate correlations ##
############################
corFunc <- function(hpaDat, hMixDat, hMixName) {
  corMat <- cor(hpaDat, as.matrix(hMixDat), method = "s")
  retDat <- do.call(rbind, lapply(1:nrow(corMat), function(i) {
    data.frame(
      sCor = corMat[i, ],
      Day = c(0:7, 8.25, 2 * (5:21)),
      Width = c(rep(1, 8), 1.5, rep(2, 17)),
      HPARegion = rownames(corMat)[i],
      hMix = hMixName
    )
  }))
  return(retDat)
}

ENS <- ENS[SYMB %in% topGenes]
SYMB <- SYMB[SYMB %in% topGenes]

plotDat <- rbind(
  corFunc(HPA.mat[ENS, ], dat.h10$sumDay[SYMB, ], "h10"), 
  corFunc(HPA.mat[ENS, ], dat.h85$sumDay[SYMB, ], "h85"), 
  corFunc(HPA.mat[ENS, ], dat.h100$sumDay[SYMB, ], "h100")
)


################
## Make Plots ##
################
p <- ggplot(plotDat, aes(x = Day, y = hMix, fill = sCor, width = Width)) + 
  theme_classic() + 
  facet_grid(HPARegion~.) + 
  geom_tile() + 
  scale_fill_viridis(option="magma") + 
  labs(x = "Day", y = "Species mixture / Brain region", 
       title = "Human Protein Atlas\nCorrelation", fill = "Spearman\nCor") + 
  theme(
    title = element_text(size = 12, face = "bold.italic"),
    axis.title = element_text(size = 10, face = "bold"),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8),
    strip.text.y.right = element_text(size = 8, angle = 0, face = "italic")
  )

jpeg("./HPA_AllGeneTriplicate.jpeg", width = 7, height = 10.5, units = "in", res = 300, quality = 0.95)
plot(p)
dev.off()












