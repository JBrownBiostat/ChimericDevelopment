###########################################################
## Human Protein Atlas Correlation with observed results ##
###########################################################
WorkDir <- "~/Documents/"
setwd(paste0(WorkDir, "/HumanProteinAtlas/"))


###############
## Libraries ##
###############
library(ggplot2)
library(gridExtra)
library(grid)
library(ggpubr)
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
  corFunc(HPA.mat[ENS, subTypeInd], dat.h10$sumDay[SYMB, ], "H10"), 
  corFunc(HPA.mat[ENS, subTypeInd], dat.h85$sumDay[SYMB, ], "H85"), 
  corFunc(HPA.mat[ENS, subTypeInd], dat.h100$sumDay[SYMB, ], "H100")
)


################
## Make Plots ##
################
plotDat$hMix <- factor(plotDat$hMix, ordered = T, levels = c("H10", "H85", "H100"))
p <- ggplot(plotDat, aes(x = Day, y = hMix, fill = sCor, width = Width)) + 
  theme_classic() + 
  facet_grid(HPARegion~.) + 
  geom_tile() + 
  scale_fill_viridis(option="magma") + 
  labs(x = "Day", y = "Species mixture / Brain region", 
       title = "Human Protein Atlas\nCorrelation", fill = "Spearman\nCor") + 
  themeObj + 
  theme(
    strip.text.y.right = element_text(size = 8, angle = 0, face = "italic")
  )
p <- annotate_figure(p, fig.lab = "A", fig.lab.size = 14, fig.lab.face = "bold")

load("DistanceFig.RData")
p.dist <- annotate_figure(p.dist, fig.lab = "B", fig.lab.size = 14, fig.lab.face = "bold")

blankPlot <- ggplot()+geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

jpeg("./HPA_AllGeneTriplicate.jpeg", width = 8, height = 7, units = "in", res = 300, quality = 0.95)
grid.arrange(p, p.dist,
             ncol = 1,
             heights = c(2.5, 1))
dev.off()






