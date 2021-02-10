########################
## FACS enrichment qc ##
########################
WorkDir <- "~/Documents/"
setwd(WorkDir)


###############
## Libraries ##
###############
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(ggpubr))
suppressMessages(library(grid))
suppressMessages(library(splines))


##########
## Data ##
##########
load("FACS_Data/EnrichmentResults.RData")
load("FACS_Data/misalignEnrich.RData")
load("FACS_Data/trendyFit.RData")
load("FACS_Data/GeneClasses.RData")
load("FACS_Data/NormDat.RData")

EU.genes <- c("STMN2", "SYT4", "NEFL")
EP.genes <- c("ASCL1", "NGFR", "NHLH2")

exptCol <- c(H10 = "deeppink", H85 = "darkorchid", H100 = "cornflowerblue", M100 = "darkolivegreen2")
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
plotGO <- function(GOdat, nodeTrim = 100){
  GOdat <- cleanGOdat(GOdat, nodeTrim)
  GOdat$FDR <- log10(GOdat$FDR)
  if(any(is.na(GOdat$FDR))) {
    GOdat$FDR[is.na(GOdat$FDR)] <- min(GOdat$FDR[!is.na(GOdat$FDR)])
  }
  xRange <- max(GOdat$FDR) - min(GOdat$FDR)
  p <- ggplot(GOdat, aes(x = FDR, y = forcats::fct_rev(Term))) +
    theme_linedraw() +
    geom_point() +
    xlim(min(GOdat$FDR) - 0.1 * xRange,
         max(GOdat$FDR) + 0.1 * xRange) + 
    themeObj + 
    theme(axis.ticks.y = element_blank(),
          strip.text.y = element_text(size=6,face="bold",angle=0),
          plot.background = element_rect(fill="white", colour="white"),
          axis.text=element_text(size=6))+
    labs(x="FDR (log 10)",y="",title="Top terms")
  return(p)
}

cleanGOdat <- function(GOdat, nodeTrim = 100) {
  GOdat$hitPct <- GOdat$Significant / GOdat$Annotated
  GOdat <- GOdat[GOdat$Annotated >= nodeTrim, ]
  GOdat <- GOdat[1:10, ]
  GOdat$Term <- gsub("...", "", GOdat$Term, fixed = TRUE)
  GOdat$Term <- unlist(lapply(strsplit(GOdat$Term, " "), function(x) {
    n <- cumsum(nchar(x))
    n <- n + seq_len(length(n)) - 1
    if(length(x) == 1){return(x)}
    ret <- x[1]
    m <- 0
    for(i in 2:length(x)) {
      if(n[i - 1] >= m + 14) {
        m <- n[i - 1]
        ret <- c(ret, "\n", x[i])
      } else {
        ret <- c(ret, " ", x[i])
      }
    }
    return(paste(ret, collapse = ""))
  }))
  while(any(duplicated(GOdat$Term))) {
    dupTerms <- which(GOdat$Term == GOdat$Term[duplicated(GOdat$Term)])
    GOdat$Term[dupTerms] <- paste0(GOdat$Term[dupTerms], "\n(", GOdat$GO.ID[dupTerms], ")")
  }
  GOdat$Term <- factor(as.character(GOdat$Term), 
                       levels = as.character(GOdat$Term), 
                       ordered = TRUE)
  return(GOdat)
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

getStartUp <- function(trendyList, genes) {
  startUp <- rep(0,length(genes))
  names(startUp) <- genes
  for(i in genes){
    subTrend <- trendyList[[i]]
    subTrends <- subTrend$Segment.Trends
    if(any(subTrends==1)) {
      if(subTrends[1]!=1) {
        startUp[i] <- subTrend$Breakpoints[which(subTrends==1)[1]-1]
      }
    }
  }
  return(startUp)
}

getPeaks <- function(trendyList, genes) {
  peakTimes <- rep(0,length(genes))
  names(peakTimes) <- genes
  for(i in genes){
    subTrend <- trendyList[[i]]
    subTrends <- subTrend$Segment.Trends
    if(length(subTrends) > 1) {
      for(j in 1:(length(subTrends)-1)) {
        if(subTrends[j]==1 & subTrends[j+1] %in% c(0,-1)) {
          peakTimes[i] <- subTrend$Breakpoints[j]
          break
        }
      }
    }
  }
  return(peakTimes)
}

get_legend <- function(gPlot){ 
  tempPlot <- ggplot_gtable(ggplot_build(gPlot)) 
  legendTemp <- which(sapply(tempPlot$grobs, function(x) x$name) == "guide-box") 
  legend <- tempPlot$grobs[[legendTemp]] 
  return(legend)
} 

plotGenes.func <- function(geneVec, dat.ref, dat.10, dat.m, 
                           trendy.ref, trendy.10, trendy.m, 
                           title = "", EU = TRUE) {
  pList <- list()
  
  if(EU) {
    eventTime <- getStartUp(trendy.10, geneVec)
    eventTime <- eventTime + 
      sapply(geneVec, function(g) {dat.10$Day[which.max(trendy.10[[g]]$Fitted.Values)]}) * 
      1e-6
  } else {
    eventTime <- getPeaks(trendy.10, geneVec)
  }
  geneVec <- geneVec[order(eventTime)]
  
  for(g in geneVec) {
    mGene <- rownames(dat.m$EC)[
      which(toupper(rownames(dat.m$EC)) == g)
    ]
    plotDat.Obs <- data.frame(
      x = c(dat.ref$Day, dat.10$Day, dat.m$Day),
      y = c(dat.ref$EC[g, ], dat.10$EC[g, ], 
            dat.m$EC[mGene, ]),
      expt = factor(c(rep("H100", length(dat.ref$Day)),
                      rep("H10", length(dat.10$Day)),
                      rep("M100", length(dat.m$Day))),
                    ordered = TRUE, levels = c("H100", "H10", "M100"))
    )
    plotDat.Fit <- data.frame(
      x = c(dat.ref$Day, dat.10$Day, dat.m$Day),
      y = c(trendy.ref[[g]]$Fitted.Values,
            trendy.10[[g]]$Fitted.Values,
            trendy.m[[mGene]]$Fitted.Values),
      expt = factor(c(rep("H100", length(dat.ref$Day)),
                      rep("H10", length(dat.10$Day)),
                      rep("M100", length(dat.m$Day))),
                    ordered = TRUE, levels = c("H100", "H10", "M100"))
    )
    
    p <- ggplot(plotDat.Obs, aes(x = x, y = y, color = expt)) + 
      theme_classic() + 
      geom_point(size = 0.5) + 
      geom_line(data = plotDat.Fit) +
      scale_color_manual(values = exptCol, 
                         labels = c("sH100", "sH10", "sM100")) + 
      themeObj + 
      theme(panel.spacing = unit(0.05, "lines"),
            axis.title.y=element_blank(),
            axis.text.y=element_text(angle = -30, hjust = 0.5),
            strip.text.y.left=element_text(size=6.5, angle = 0),
            plot.background=element_rect(fill="white", colour="white"),
            plot.margin = unit(c(0.1, 0.1, 0.1, 0.5), "cm")) +
      labs(x = "Day", y = "Expression", title = g, color = "Species\nmixture") + 
      guides(colour = guide_legend(override.aes = list(size = 1.5)))
    pList[[g]] <- p
  }
  
  p.legend <- get_legend(pList[[1]])
  p <- grid.arrange(pList[[1]] + theme(legend.position = "none"),
                    pList[[2]] + theme(legend.position = "none"),
                    pList[[3]] + theme(legend.position = "none"),
                    nrow = 1,
                    top = textGrob(title, gp = gpar(fontsize = 12)))
  
  return(list(p = p, p.legend = p.legend))
}


#####################
## Plot enrichment ##
#####################
misGO.H <- plotGO(hMix.misEnrich)

misGO.M <- plotGO(mMix.misEnrich)

EUGO <- plotGO(H100_HMix_AccDec_GO$accelerated$earlyUp$topResults)


###############
## Plot acc. ##
###############
upMixDat <- checkUpTrend(h100FACS33.trendy, hMixFACS33.trendy)
plotMixDat <- data.frame(
  gene = upMixDat$gene,
  H100Time = upMixDat$refBreak,
  HMixTime = upMixDat$testBreak
)
plotMixDat <- plotMixDat[plotMixDat$HMixTime > 0 & plotMixDat$H100Time > 0, ]
plotMixDat$HMixFC <- plotMixDat$H100Time / plotMixDat$HMixTime
mixMod <- lm(HMixFC ~ bs(HMixTime), data = plotMixDat)
mix_meanAcc <- median(mixMod$fitted.values[plotMixDat$HMixTime <= 16])
pMix <- ggplot(plotMixDat) + 
  theme_classic() + 
  geom_hline(yintercept = 1, color = "black") + 
  geom_point(aes(x = HMixTime, y = HMixFC), size = 0.25, alpha = 0.5) + 
  geom_smooth(aes(x = HMixTime, y = HMixFC), method = "lm", formula = y ~ bs(x), lwd = 1) + 
  labs(x = "Human Mix event time", y = "Acc. factor", 
       title = paste0("sH10 Acceleration\n(rel. to sH100)\nMedian Acc. = ", round(mix_meanAcc, 3))) + 
  ylim(0, quantile(plotMixDat$HMixFC, 0.98)) + xlim(c(0, 26)) + 
  themeObj


##################
## Build figure ##
##################
load("FACS_Data/MisalignFigs.RData")
p.H <- p.H + 
  scale_color_viridis_d(end = 0.75, labels = c("sH100", "sH10")) + 
  labs(x = "Day", y = "Total error", 
       title = "Sorted human\nempirical misalignment") + 
  ylim(c(0, 0.6))
p.H <- annotate_figure(p.H, fig.lab = "A", fig.lab.size = 14, fig.lab.face = "bold")

EUPlot <- plotGenes.func(EU.genes, h100.facs33, hMix.facs33, m100.facs33,
                         h100FACS33.trendy, hMixFACS33.trendy, m100FACS33.trendy,
                         "")
EPPlot <- plotGenes.func(EP.genes, h100.facs33, hMix.facs33, m100.facs33,
                         h100FACS33.trendy, hMixFACS33.trendy, m100FACS33.trendy,
                         "", F)
genePlot <- grid.arrange(EUPlot$p, EPPlot$p, EUPlot$p.legend, 
                         layout_matrix = rbind(c(1, 3),
                                               c(2, 3)),
                         widths = c(1, 0.2))
genePlot <- annotate_figure(genePlot, fig.lab = "C", fig.lab.size = 14, fig.lab.face = "bold")

H10_Unchanged <- checkUpTrend(h100FACS33.trendy, hMixFACS33.trendy)
H10_Unchanged <- H10_Unchanged$gene[H10_Unchanged$refBreak > 0 |
                                      H10_Unchanged$testBreak > 0]
countDat <- data.frame(
  nGenes = c(length(H100_HMix_Facs_AccDec$accelerated$earlyUp), 
             length(H100_HMix_Facs_AccDec$decelerated$lateUp),
             length(H10_Unchanged) - 
               length(H100_HMix_Facs_AccDec$accelerated$earlyUp) - 
               length(H100_HMix_Facs_AccDec$decelerated$lateUp)),
  Comparison = factor(c("Early", "Late", "Unchanged"), ordered = T,
                      levels = c("Unchanged", "Late", "Early")),
  Mixture = rep(c("M15"), 3)
)
countDat$Comparison <- factor(countDat$Comparison, ordered = T, levels = countDat$Comparison)
countDat$ypos <- rev(cumsum(rev(countDat$nGenes)) - 0.5 * rev(countDat$nGenes))
p.countPie <- ggplot(countDat, aes(x = "", y = nGenes, fill = Comparison)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  geom_text(aes(y = ypos, label = nGenes), color = "white", size = 4) +
  scale_fill_viridis_d(end = 0.8) + 
  themeObj + 
  theme(axis.text = element_blank(), 
        legend.position = "bottom", 
        legend.title = element_text(hjust = 0.5),
        legend.text = element_text(angle = 0, hjust = 0.5)) +
  guides(
    fill = guide_legend(title.position = "top", label.position = "top")
  ) +
  labs(x = "", y = "")
p.countPie <- annotate_figure(p.countPie, fig.lab = "D", fig.lab.size = 14, fig.lab.face = "bold")

misFig <- grid.arrange(misGO.H + labs(title = "sH10 missaligned\nenrichment"), 
                       misGO.M + labs(title = "sM90 missaligned\nenrichment"), 
                       nrow = 1)
misFig <- annotate_figure(misFig, fig.lab = "B", fig.lab.size = 14, fig.lab.face = "bold")

accFig <- grid.arrange(pMix, 
                       EUGO + labs(title = "Early Up enrichment"), 
                       nrow = 1)
accFig <- annotate_figure(accFig, fig.lab = "E", fig.lab.size = 14, fig.lab.face = "bold")

jpeg("FACS_Data/MisalignQC.jpeg", width = 8, height = 10,
     units = "in", res = 300, quality = 0.95)
grid.arrange(p.H, misFig, genePlot, p.countPie, accFig, 
             layout_matrix = rbind(
               rep(1:2, c(1, 2)),
               rep(3, 3),
               rep(4:5, c(1, 2))
             ), 
             heights = c(1, 1.5, 1))
dev.off()






