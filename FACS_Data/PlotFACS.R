######################################
## Plot figures - Trendy enrichment ##
######################################
WorkDir <- "~/Documents/"
setwd(WorkDir)


###############
## Libraries ##
###############
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(gridExtra))
suppressMessages(library(grid))
suppressMessages(library(topGO))
suppressMessages(library(viridis))
suppressMessages(library(stringr))
suppressMessages(library(png))
suppressMessages(library(dplyr))
suppressMessages(library(yarrr))
suppressMessages(library(wesanderson))
suppressMessages(library(org.Mm.eg.db))
suppressMessages(library(splines))



##########
## Data ##
##########
exptCol <- c(HMix = "deeppink", H100 = "cornflowerblue")

load("FACS_Data/NormDat.RData")
load("FACS_Data/GeneClasses.RData")
load("FACS_Data/trendyFit.RData")
load("FACS_Data/EnrichmentResults.RData")

## Genes by term ##
allHuman.trendy <- sort(unique(c(
  names(h100FACS33.trendy),
  names(hMixFACS33.trendy)
)))
geneVec <- factor(rep(1, length(allHuman.trendy)), levels = c("0", "1"))
names(geneVec) <- allHuman.trendy
GOdata <- new("topGOdata", ontology = "BP", allGenes = geneVec,
              annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "Symbol", nodeSize = 20)
allGO.human <- genesInTerm(GOdata)

themeObj <- theme(
  title = element_text(face = "bold", size = 10),
  axis.title = element_text(size = 9),
  axis.text = element_text(size = 8)
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
  p <- ggplot(GOdat, aes(x = FDR, y = forcats::fct_rev(Term), size = Annotated)) +
    theme_linedraw() +
    geom_point() +
    scale_size(trans="log10") +
    labs(size = "GO-term\nsize") +
    theme(plot.title=element_text(size=10,face="bold.italic"),
          axis.title.x=element_text(size=8,face="bold"),
          axis.title.y=element_text(size=8,face="bold"),
          legend.title=element_text(size=6,face="bold"),
          legend.text=element_text(size=6),
          axis.text=element_text(size=6),
          axis.ticks.y=element_blank(),
          strip.text.y=element_text(size=6,face="bold",angle=0),
          plot.background=element_rect(fill="white", colour="white"))+
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
      if(n[i - 1] >= m + 12) {
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

plotPeak <- function(GOdat, peakGenes, allGenes, ref.trendy, test.trendy, nodeTrim = 100) {
  GOdat <- cleanGOdat(GOdat,nodeTrim)
  plotDat <- getDensDat(GOdat, ref.trendy, test.trendy, peakGenes, allGenes, getPeaks)
  plotDat$densDat$ref_test = factor(plotDat$densDat$ref_test, ordered = TRUE, levels = c("ref", "test"))
  p <- ggplot(plotDat$densDat, aes(x = x, y = y, color = ref_test, fill = ref_test))+
    theme_minimal()+
    facet_grid(term ~ ., switch = "y")+
    geom_area(alpha=0.4,position="identity")+
    geom_line(size=1)+
    expand_limits(x=0)+
    theme(plot.title=element_text(size=10,face="bold.italic"),
          panel.spacing = unit(0.05, "lines"),
          axis.title.x=element_text(size=8,face="bold"),
          axis.title.y=element_blank(),
          legend.title=element_text(size=8,face="bold"),
          legend.text=element_text(size=6.5),
          axis.text=element_text(size=6.5),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          strip.text.y.left=element_text(size=6.5, angle = 0),
          plot.background=element_rect(fill="white", colour="white"))+
    labs(x="Day",y="",title="Time of peak",color="Reference/\nTest",fill="Reference/\nTest") + 
    scale_color_manual(values = unname(exptCol[c(2, 1)]), label = c("H100", "HMix")) + 
    scale_fill_manual(values = unname(exptCol[c(2, 1)]), label = c("H100", "HMix"))
  ksDat <- plotDat$ksDat
  ksDat$xMax = max(ggplot_build(p)$data[[1]]$x)
  ksDat$xMin = min(ggplot_build(p)$data[[1]]$x)
  ksDat$yMax = max(ggplot_build(p)$data[[1]]$y)
  ksDat$yMin = min(ggplot_build(p)$data[[1]]$y)
  p <- p +
    geom_text(data=ksDat,aes(x=xMax,y=yMax,label=pVal),size=3,color="black",hjust=0,
              inherit.aes=FALSE, nudge_x=-0.2 * (ksDat$xMax[1] - ksDat$xMin[1]),
              nudge_y=-0.5 * (ksDat$yMax[1] - ksDat$yMin[1]))
  return(p)
}

getDensDat <- function(GOdat, ref.trendy, test.trendy, sigGenes, allGenes, eventTime.func,
                       testShiftLeft = TRUE, densLow = 0, densHigh = 42) {
  densDat <- do.call(rbind, lapply(seq_len(nrow(GOdat)), function(i) {
    genes <- intersect(sigGenes, allGenes[[GOdat$GO.ID[i]]])
    refEvent <- eventTime.func(ref.trendy, genes)
    testEvent <- eventTime.func(test.trendy, genes)
    ref.dens <- tryCatch({
      density(refEvent, n = 2^10, bw = "SJ", from = densLow, to = densHigh)
    }, error = function(e) {
      density(refEvent, bw = ifelse(length(unique(refEvent)) > 1, "nrd0", 1/2),
              from = densLow, to = densHigh, n = 2^10)
    })
    test.dens <- tryCatch({
      density(testEvent, n = 2^10, bw = "SJ", from = densLow, to = densHigh)
    }, error = function(e) {
      density(testEvent, bw = ifelse(length(unique(testEvent)) > 1, "nrd0", 1/2),
              from = densLow, to = densHigh, n = 2^10)
    })
    maxDens <- max(c(ref.dens$y, test.dens$y))
    ref.dens$y <- ref.dens$y / maxDens
    test.dens$y <- test.dens$y / maxDens
    if(max(ref.dens$y) < 0.2) {
      ref.dens$y <- ref.dens$y / (max(ref.dens$y) / 0.2)
    }
    if(max(test.dens$y) < 0.2) {
      test.dens$y <- test.dens$y / (max(test.dens$y) / 0.2)
    }
    return(data.frame(
      x = c(ref.dens$x, test.dens$x),
      y = c(ref.dens$y, test.dens$y),
      ref_test = rep(c("ref", "test"), each = 2^10),
      term = GOdat$Term[i],
      pVal = suppressWarnings(round(ks.test(x=testEvent,y=refEvent,
                                            alternative=ifelse(testShiftLeft, "greater", "less"))$p.value,3))
    ))
  }))
  ksDat <- do.call(rbind, lapply(seq_len(nrow(GOdat)), function(i) {
    term <- GOdat$Term[i]
    return(data.frame(
      x = densHigh,
      y = 0,
      term = term,
      pVal = paste0("p=",densDat[densDat$term==term,]$pVal[1])
    ))
  }))
  densDat <- densDat[,-5]
  return(list(densDat = densDat, ksDat = ksDat))
}

getUpSlope <- function(trendyList, genes) {
  upSlope <- rep(0,length(genes))
  names(upSlope) <- genes
  for(i in genes){
    subTrend <- trendyList[[i]]
    subTrends <- subTrend$Segment.Trends
    if(any(subTrends==1)) {
      upSlope[i] <- subTrend$Segment.Slopes[which(subTrends==1)[1]]
    }
  }
  return(upSlope)
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

plotStartUp <- function(GOdat,upTrendGenes,allGenes,ref.trendy,test.trendy,nodeTrim=100, addP = TRUE, testH10 = TRUE) {
  GOdat <- cleanGOdat(GOdat,nodeTrim)
  plotDat <- getDensDat(GOdat, ref.trendy, test.trendy, upTrendGenes, allGenes, getStartUp)
  plotDat$densDat$ref_test = factor(plotDat$densDat$ref_test, ordered = TRUE, levels = c("ref", "test"))
  p <- ggplot(plotDat$densDat, aes(x = x, y = y, color = ref_test, fill = ref_test))+
    theme_minimal() +
    facet_grid(term ~ ., switch = "y") +
    geom_area(alpha=0.4,position="identity") +
    geom_line(size=1) +
    expand_limits(x=0) +
    theme(plot.title=element_text(size=10,face="bold.italic"),
          panel.spacing = unit(0.05, "lines"),
          axis.title.x=element_text(size=8,face="bold"),
          axis.title.y=element_blank(),
          legend.title=element_text(size=8,face="bold"),
          legend.text=element_text(size=6.5),
          axis.text=element_text(size=6.5),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          strip.text.y.left=element_text(size=6.5, angle = 0),
          plot.background=element_rect(fill="white", colour="white")) +
    labs(x="Day",y="",title="Start of up-trend",color="",fill="") + 
    scale_color_manual(values = unname(exptCol[c(2:1)]), label = c("H100", "HMix")) + 
    scale_fill_manual(values = unname(exptCol[c(2:1)]), label = c("H100", "HMix"))
  ksDat <- plotDat$ksDat
  ksDat$xMax = max(ggplot_build(p)$data[[1]]$x)
  ksDat$xMin = min(ggplot_build(p)$data[[1]]$x)
  ksDat$yMax = max(ggplot_build(p)$data[[1]]$y)
  ksDat$yMin = min(ggplot_build(p)$data[[1]]$y)
  p <- p +
    geom_text(data=ksDat,aes(x=xMax,y=yMax,label=pVal),size=3,color="black",hjust=0,
              inherit.aes=FALSE, nudge_x=-0.2 * (ksDat$xMax[1] - ksDat$xMin[1]),
              nudge_y=-0.5 * (ksDat$yMax[1] - ksDat$yMin[1]))
  return(p)
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


#####################
## Plot Early-Peak ##
#####################
p.GOlate <- plotGO(H100_HMix_AccDec_GO$decelerated$peak$topResults)
p.GOlate <- annotate_figure(p.GOlate, fig.lab = "Late Peak", fig.lab.size = 14, fig.lab.face = "bold")
p.GOearly <- plotGO(H100_HMix_AccDec_GO$accelerated$peak$topResults)
p.GOearly <- annotate_figure(p.GOearly, fig.lab = "Early Peak", fig.lab.size = 14, fig.lab.face = "bold")

p.peakTime_Late <- plotPeak(H100_HMix_AccDec_GO$decelerated$peak$topResults,
                            H100_HMix_Facs_AccDec$decelerated$peak,
                            allGO.human, h100FACS33.trendy, hMixFACS33.trendy)
p.peakTime_Late <- annotate_figure(p.peakTime_Late, fig.lab = "Late Peak", fig.lab.size = 14, fig.lab.face = "bold")
p.peakTime_Early <- plotPeak(H100_HMix_AccDec_GO$accelerated$peak$topResults,
                             H100_HMix_Facs_AccDec$accelerated$peak,
                             allGO.human, h100FACS33.trendy, hMixFACS33.trendy)
p.peakTime_Early <- annotate_figure(p.peakTime_Early, fig.lab = "Early Peak", fig.lab.size = 14, fig.lab.face = "bold")

countDat <- data.frame(
  nGenes = c(length(H100_HMix_Facs_AccDec$accelerated$peak), 
             length(H100_HMix_Facs_AccDec$decelerated$peak),
             length(checkUpTrend(h100FACS33.trendy, hMixFACS33.trendy)$gene) - 
               length(H100_HMix_Facs_AccDec$accelerated$peak) - 
               length(H100_HMix_Facs_AccDec$decelerated$peak)),
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
  scale_fill_viridis_d(end = 0.8)

peakH100Dat <- checkPeak(h100FACS33.trendy, hMixFACS33.trendy)
hPeakDat <- data.frame(
  gene = peakH100Dat$gene,
  H100Peak = peakH100Dat$refBreak,
  HMixPeak = peakH100Dat$testBreak, 
  HMixFC = peakH100Dat$refBreak / peakH100Dat$testBreak
)
pH <- ggplot(hPeakDat) + 
  theme_classic() + 
  geom_hline(yintercept = 1, color = "black") + 
  geom_point(aes(x = HMixPeak, y = HMixFC), color = exptCol[2], size = 0.75, alpha = 0.75) + 
  geom_smooth(aes(x = HMixPeak, y = HMixFC), color = exptCol[2], method = "lm", formula = y ~ bs(x), lwd = 1) + 
  labs(x = "HMix", y = "Acc. factor", title = "HMix Acceleration\nrelative to H100") + 
  ylim(quantile(hPeakDat$HMixFC, c(0.01, 0.99))) + xlim(c(2, 33)) + 
  themeObj
pH

jpeg(filename = "FACS_Data/FacsPeak.jpeg", height = 10.5, width = 8, units = "in", res = 300, quality = 0.9)
grid.arrange(p.countPie, pH, p.GOearly, p.GOlate, p.peakTime_Early, p.peakTime_Late,
             layout_matrix = rbind(c(1, 2),
                                   c(3, 4),
                                   c(5, 6)),
             heights = c(0.7, 0.8, 1))
dev.off()


#####################
## Plot Early-Up ##
#####################
p.GOlate <- plotGO(H100_HMix_AccDec_GO$decelerated$lateUp$topResults)
p.GOlate <- annotate_figure(p.GOlate, fig.lab = "Late Up", fig.lab.size = 14, fig.lab.face = "bold")
p.GOearly <- plotGO(H100_HMix_AccDec_GO$accelerated$earlyUp$topResults)
p.GOearly <- annotate_figure(p.GOearly, fig.lab = "Early Up", fig.lab.size = 14, fig.lab.face = "bold")

p.upTime_Late <- plotStartUp(H100_HMix_AccDec_GO$decelerated$lateUp$topResults,
                             H100_HMix_Facs_AccDec$decelerated$lateUp,
                             allGO.human, h100FACS33.trendy, hMixFACS33.trendy)
p.upTime_Late <- annotate_figure(p.upTime_Late, fig.lab = "Late Up", fig.lab.size = 14, fig.lab.face = "bold")
p.upTime_Early <- plotStartUp(H100_HMix_AccDec_GO$accelerated$earlyUp$topResults,
                              H100_HMix_Facs_AccDec$accelerated$earlyUp,
                              allGO.human, h100FACS33.trendy, hMixFACS33.trendy)
p.upTime_Early <- annotate_figure(p.upTime_Early, fig.lab = "Early Up", fig.lab.size = 14, fig.lab.face = "bold")

countDat <- data.frame(
  nGenes = c(length(H100_HMix_Facs_AccDec$accelerated$earlyUp), 
             length(H100_HMix_Facs_AccDec$decelerated$lateUp),
             length(checkPeak(h100FACS33.trendy, hMixFACS33.trendy)$gene) - 
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
  scale_fill_viridis_d(end = 0.8)

peakH100Dat <- checkPeak(h100FACS33.trendy, hMixFACS33.trendy)
hPeakDat <- data.frame(
  gene = peakH100Dat$gene,
  H100Peak = peakH100Dat$refBreak,
  HMixPeak = peakH100Dat$testBreak, 
  HMixFC = peakH100Dat$refBreak / peakH100Dat$testBreak
)
pH <- ggplot(hPeakDat) + 
  theme_classic() + 
  geom_hline(yintercept = 1, color = "black") + 
  geom_point(aes(x = HMixPeak, y = HMixFC), color = exptCol[2], size = 0.75, alpha = 0.75) + 
  geom_smooth(aes(x = HMixPeak, y = HMixFC), color = exptCol[2], method = "lm", formula = y ~ bs(x), lwd = 1) + 
  labs(x = "HMix", y = "Acc. factor", title = "HMix Acceleration\nrelative to H100") + 
  ylim(quantile(hPeakDat$HMixFC, c(0.01, 0.99))) + xlim(c(2, 40)) + 
  themeObj
pH

jpeg(filename = "FACS_Data/FacsUp.jpeg", height = 10.5, width = 8, units = "in", res = 300, quality = 0.9)
grid.arrange(p.countPie, pH, p.GOearly, p.GOlate, p.upTime_Early, p.upTime_Late,
             layout_matrix = rbind(c(1, 2),
                                   c(3, 4),
                                   c(5, 6)),
             heights = c(0.7, 0.8, 1))
dev.off()























