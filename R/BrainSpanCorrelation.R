##################################################
## Brain-Span Correlation with observed results ##
##################################################

setwd("~/Documents/Research/Chris_Mixing_2019/Final_MixData_March25/BrainSpan/")


###############
## Libraries ##
###############
library(ggplot2)
library(gridExtra)
library(viridis)
library(topGO)


##########
## Data ##
##########
load("../PreprocessingAndNormalization/NormDat.RData")
load("../Trendy/GeneClasses.RData")
load("../Trendy/trendyFit.RData")
load("../Trendy/EnrichmentResults.RData")

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
AG.dat <- AG.dat[,c(5*(1:6)-4,5*(1:6)-3,5*(1:6)-2,5*(1:6)-1,5*(1:6))]
rownames(AG.dat) <- AG.row$gene_symbol
AG.dat <- AG.dat[which(rownames(AG.dat)%in%rownames(dat.h100$sumDay)),]
AG.dat <- AG.dat[-which(duplicated(rownames(AG.dat))),]
AG.dat <- AG.dat[order(rownames(AG.dat)),]


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
maxCV <- maxCV[intersect(names(h100.trendy), rownames(AG.dat))]

topGenes <- names(maxCV)[order(maxCV, decreasing = T)[1:1000]]


############################
## Calculate correlations ##
############################

corFunc <- function(bsDat, hMixDat, hMixName) {
  corMat <- cor(bsDat, as.matrix(hMixDat), method = "s")
  retDat <- do.call(rbind, lapply(1:nrow(corMat), function(i) {
    data.frame(
      sCor = corMat[i, ],
      Day = c(0:7, 8.25, 2 * (5:21)),
      Width = c(rep(1, 8), 1.5, rep(2, 17)),
      BSRegion = rownames(corMat)[i],
      hMix = hMixName
    )
  }))
  return(retDat)
}

week9Col <- c("AMY.9", "DFC.9", "HIP.9", "MFC.9", "OFC.9")
plotDat.week9 <- rbind(
  corFunc(AG.dat[topGenes, week9Col], dat.h10$sumDay[topGenes, ], "H10"), 
  corFunc(AG.dat[topGenes, week9Col], dat.h85$sumDay[topGenes, ], "H85"), 
  corFunc(AG.dat[topGenes, week9Col], dat.h100$sumDay[topGenes, ], "H100")
)

week8Col <- c("AMY.8", "DFC.8", "HIP.8", "MFC.8", "OFC.8")
plotDat.week8 <- rbind(
  corFunc(AG.dat[topGenes, week8Col], dat.h10$sumDay[topGenes, ], "H10"), 
  corFunc(AG.dat[topGenes, week8Col], dat.h85$sumDay[topGenes, ], "H85"), 
  corFunc(AG.dat[topGenes, week8Col], dat.h100$sumDay[topGenes, ], "H100")
)

week12Col <- c("AMY.12", "DFC.12", "HIP.12", "MFC.12", "OFC.12")
plotDat.week12 <- rbind(
  corFunc(AG.dat[topGenes, week12Col], dat.h10$sumDay[topGenes, ], "H10"), 
  corFunc(AG.dat[topGenes, week12Col], dat.h85$sumDay[topGenes, ], "H85"), 
  corFunc(AG.dat[topGenes, week12Col], dat.h100$sumDay[topGenes, ], "H100")
)


################
## Make Plots ##
################
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

p.week9 <- ggplot(plotDat.week9, aes(x = Day, y = hMix, fill = sCor, width = Width)) + 
  theme_classic() + 
  facet_grid(BSRegion~.) + 
  geom_tile() + 
  scale_fill_viridis(option="magma") + 
  labs(x = "Day", y = "Species mixture / Brain Region", 
       title = "Brain Span Correlation\n9 PCW", fill = "Spearman\nCor") + 
  theme(
    title = element_text(size = 9, face = "bold.italic"),
    axis.title = element_text(size = 7, face = "bold"),
    axis.text = element_text(size = 6),
    legend.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 6),
    strip.text.y.right = element_text(size = 6, face = "italic"),
    plot.margin = unit(c(0.05, 0.1, 0.05, 0.1), "in")
  )
LP <- get_legend(p.week9)
p.week9 <- p.week9 + theme(legend.position = "none")

p.week8 <- ggplot(plotDat.week8, aes(x = Day, y = hMix, fill = sCor, width = Width)) + 
  theme_classic() + 
  facet_grid(BSRegion~.) + 
  geom_tile() + 
  scale_fill_viridis(option="magma") + 
  labs(x = "Day", y = "Species mixture / Brain Region", 
       title = "Brain Span Correlation\n8 PCW", fill = "Spearman\nCor") + 
  theme(
    title = element_text(size = 9, face = "bold.italic"),
    axis.title = element_text(size = 7, face = "bold"),
    axis.text = element_text(size = 6),
    legend.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 6),
    legend.position = "none",
    strip.text.y.right = element_text(size = 6, face = "italic"),
    plot.margin = unit(c(0.05, 0.1, 0.05, 0.1), "in")
  )

p.week12 <- ggplot(plotDat.week12, aes(x = Day, y = hMix, fill = sCor, width = Width)) + 
  theme_classic() + 
  facet_grid(BSRegion~.) + 
  geom_tile() + 
  scale_fill_viridis(option="magma") + 
  labs(x = "Day", y = "Species mixture / Brain Region", 
       title = "Brain Span Correlation\n12 PCW", fill = "Spearman\nCor") + 
  theme(
    title = element_text(size = 9, face = "bold.italic"),
    axis.title = element_text(size = 7, face = "bold"),
    axis.text = element_text(size = 6),
    legend.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 6),
    legend.position = "none",
    strip.text.y.right = element_text(size = 6, face = "italic"),
    plot.margin = unit(c(0.05, 0.1, 0.05, 0.1), "in")
  )

jpeg("./BS_AllGeneTriplicate.jpeg", width = 8, height = 6, units = "in", res = 300, quality = 0.95)
grid.arrange(p.week8, p.week9, p.week12, LP, ncol = 4, widths = c(1, 1, 1, 0.2))
dev.off()





