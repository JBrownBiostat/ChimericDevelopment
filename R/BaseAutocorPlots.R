#################################
## Basic autocorrelation plots ##
#################################


###############
## Libraries ##
###############
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(scales))
suppressMessages(library(reshape2))
suppressMessages(library(viridis))
suppressMessages(library(topGO))


##########
## Data ##
##########
setwd("~/Documents/Research/Chris_Mixing_2019/Final_MixData_March25/")
load("PreprocessingAndNormalization/NormDat.RData")
load("Trendy/GeneClasses.RData")
load("Trendy/trendyFit.RData")
load("Trendy/EnrichmentResults.RData")


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

topGenes <- names(h100.trendy)[order(maxCV, decreasing = T)[1:1000]]


#######################
## Correlation plots ##
#######################

corPlot.func <- function(testDat,refDat,saveFile,title,xtitle,ytitle,corLim=NULL){
  geneVarOrd <- order(apply(cbind(testDat, refDat), 1, var), decreasing = TRUE)
  testDat <- as.matrix(testDat[geneVarOrd[1:min(5000, nrow(testDat))], ])
  refDat <- as.matrix(refDat[geneVarOrd[1:min(5000, nrow(testDat))], ])
  corMat <- cor(refDat, testDat, method = "s")
  plotDat <- cbind(melt(corMat),
                   x.cent=rep(c(0:7,8.25,2*(5:21)),26),
                   x.width=rep(c(rep(1,8),1.5,rep(2,17)),26),
                   y.cent=rep(c(0:7,8.25,2*(5:21)),each=26),
                   y.width=rep(c(rep(1,8),1.5,rep(2,17)),each=26))
  colnames(plotDat)[1:3] <- c("testInd","refInd","sCor")
  if(is.null(corLim)){corLim=quantile(plotDat$sCor,probs=c(0,0.5,1),na.rm=TRUE)}
  plotDat$plotCor <- pmax(corLim[1],pmin(corLim[3],plotDat$sCor))
  p <- ggplot(plotDat,aes(x=x.cent,y=y.cent,width=x.width,height=y.width,fill=sCor))+
    theme_classic() + 
    geom_tile()+
    scale_fill_gradientn(colors=c("darkblue","white","darkred"))+
    labs(x=xtitle,y=ytitle,fill="S-cor",title=title)+
    scale_y_continuous(breaks=4*(0:10))+
    scale_x_continuous(breaks=4*(0:10))+
    geom_abline(slope=1,intercept=0)+
    theme(
      title = element_text(size = 9, face = "bold.italic"),
      axis.title = element_text(size = 7, face = "bold"),
      axis.text = element_text(size = 6),
      legend.title = element_text(size = 7, face = "bold"),
      legend.text = element_text(size = 6)
    )
  pdf(saveFile, width = 10, height = 10)
  plot(p)
  dev.off()
  return(list(plotDat=plotDat,plot=p))
}

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

h100.null <- corPlot.func(dat.h100$sumDay[topGenes, ], dat.h100$sumDay[topGenes, ], 
                          "AutoCorrelation/BasePlots/h100null.pdf",
                          "H100 reference", "H100", "H100", c(0.0, 0.845, 0.99))
h85.null <- corPlot.func(dat.h85$sumDay[topGenes, ], dat.h85$sumDay[topGenes, ], 
                         "AutoCorrelation/BasePlots/h85null.pdf",
                         "H85 reference", "H85","H85", c(0.0, 0.845, 0.99))
h10.null <- corPlot.func(dat.h10$sumDay[topGenes, ], dat.h10$sumDay[topGenes, ], 
                         "AutoCorrelation/BasePlots/h10null.pdf",
                         "H10 reference", "H10", "H10", c(0.0, 0.845, 0.99))

h100_10 <- corPlot.func(dat.h10$sumDay[topGenes, ], dat.h100$sumDay[topGenes, ], 
                        "AutoCorrelation/BasePlots/h100_10.pdf",
                        "H10 by H100", "H100", "H10", c(0.0, 0.845, 0.99))
h100_85 <- corPlot.func(dat.h85$sumDay,dat.h100$sumDay,"AutoCorrelation/BasePlots/h100_85.pdf",
                        "H85 by H100","H100","H85",c(0.0,0.845,0.99))
h85_10 <- corPlot.func(dat.h10$sumDay[topGenes, ], dat.h85$sumDay[topGenes, ], 
                       "AutoCorrelation/BasePlots/h85_10.pdf",
                       "H10 by H85", "H85", "H10", c(0.0, 0.845, 0.99))

HGenes <- rownames(dat.h100$allDay)[rownames(dat.h100$allDay) %in% toupper(rownames(dat.m100$allDay))]
HGenes <- intersect(HGenes, topGenes)
MGenes <- rownames(dat.m100$allDay)[toupper(rownames(dat.m100$allDay)) %in% topGenes]
MGenes <- MGenes[!duplicated(toupper(MGenes))]

h100_0 <- corPlot.func(dat.m100$sumDay[MGenes, ], dat.h100$sumDay[HGenes, ], "AutoCorrelation/BasePlots/h100_0.pdf",
                       "M100 by H100","H100","M100",c(0.0,0.845,0.99))

LP <- get_legend(h100_10$plot)

save(h100_10, h100_85, h85_10, h100_0, h100.null, h85.null, h10.null,
     file = "AutoCorrelation/BasePlots/Plots.RData",
     compress = "xz")
jpeg("AutoCorrelation/BasePlots/H10Triplicate.jpeg", width = 12, height = 4, units = "in", res = 300, quality = 1)
grid.arrange(h100_10$plot + theme(legend.position = "none"), 
             h100_85$plot + theme(legend.position = "none"), 
             h10.null$plot + theme(legend.position = "none"), 
             LP,
             ncol = 4, widths = c(1, 1, 1, 0.2))
dev.off()

jpeg("AutoCorrelation/BasePlots/Fig1Triplicate.jpeg", width = 12, height = 4, units = "in", res = 300, quality = 1)
grid.arrange(h100_0$plot + theme(legend.position = "none"), 
             h100_10$plot + theme(legend.position = "none"), 
             h10.null$plot + theme(legend.position = "none"), 
             LP,
             ncol = 4, widths = c(1, 1, 1, 0.2))
dev.off()



