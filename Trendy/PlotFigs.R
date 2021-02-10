######################################
## Plot figures - Trendy enrichment ##
######################################
WorkDir <- "~/Documents/"


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
suppressMessages(library(splines))


##########
## Data ##
##########
exptCol <- c(H10 = "deeppink", H85 = "darkorchid", H100 = "cornflowerblue", M100 = "darkolivegreen2")

setwd(WorkDir)
load("PreprocessingAndNormalization/NormDat.RData")
load("Trendy/GeneClasses.RData")
load("Trendy/trendyFit.RData")
load("Trendy/EnrichmentResults.RData")

## Genes by term ##
allHuman.trendy <- sort(unique(c(
  names(h100.trendy),
  names(h85.trendy),
  names(h10.trendy)
)))
geneVec <- factor(rep(1,length(allHuman.trendy)),levels=c("0","1"))
names(geneVec) <- allHuman.trendy
GOdata <- new("topGOdata",ontology="BP",allGenes=geneVec,
              annot=annFUN.org,mapping="org.Hs.eg.db",ID="Symbol",nodeSize=20)
allGO.human <- genesInTerm(GOdata)

## Plot genes ##
genes.EP <- list(
  "Neurogenesis &\nAxon Migration" = c("GATA3", "POU3F2", "POU3F4", "SOX1", "TUBB3", "STMN1", "ASCL1", "ADCYAP1", 
                                       "BEX1", "CNTN4", "DIO3", "FLRT1", "GAP43", "GPC2", "INA",
                                       "MAB21L2", "NEFM", "PLXNA4", "SLIT1", "SLIT3", "UNC5C", "BCHE", "BRSK2",
                                       "BTG2", "CALM2", "CBFA2T2", "CDKN1C", "CNTFR", "CRMP1", "CXCR4", "DPYSL4",
                                       "ECT2", "EFNA3", "EFNB2", "EVL", "GAB2", "LZTS1", "MAPK8IP2",
                                       "MAPK8IP3", "NEDD4", "NEUROD4", "NGFR", "NHLH2", "ONECUT2", "SUN2",
                                       "TAGLN3", "ZEB1", "ZEB2"),
  "Forebrain &\nNeural Tube" = c("DLL3", "GLI3", "MEIS1"),
  "Neuron Signaling &\nSynapse Transmission" = c("APLP1", "ATCAY", "CHRNA3", "L1CAM", "LRRC4B"),
  "Ventral Midbrain" = c("GATA2", "ISL1", "NKX2-2", "NKX6-1", "NKX6-2", "PTCH1", "SHH", "LHX4", "ANKRD1",
                         "GLI1", "HOXB2")
)
genes.EP <- lapply(genes.EP, function(x) {sort(x)})

genes.EU <- list(
  "Neurogenesis &\nAxon Migration" = c("CDC42", "FGF13", "LRRN3", "RGS4", "STMN2", "NEFL", 
                                       "RGS2", "DCX", "SLITRK1", "LHX6", "NEUROG2", "FGF11", 
                                       "NRP1", "WNT5A", "FAIM2", "GPM6B", "MYT1", "GPM6A", 
                                       "CDK5R1", "NRN1", "MAPT", "NRSN1", "RIT2", "CDK5R2"),
  "Forebrain &\nNeural Tube" = c("FEZ1", "EFNB3"),
  "Neuron Signaling &\nSynapse Transmission" = c("KLHL1", "CHRNA7", "SNAP25", "SYT4", "ASIC2", "KCNAB1", 
                                                 "ATP6AP2", "CACNG7", "DLG4", "FGF14", "SYN1", "SYT3"),
  "NSC" = c("FABP7", "ST8SIA4", "FGF10"),
  "GLUTA &\nGABA" = c("SLC1A3", "GABRG2", "GRIN2D", "GABRA1")
)
genes.EU <- lapply(genes.EU, function(x) {x[!(x %in% unlist(genes.EP))]})

genes.DE <- list(
  "Neural\nTube" = c("AJUBA", "AHNAK", "FAT4", "FBXO32", "SCUBE2", "PALLD", "SHROOM3"), 
  "Anterior\nForebrain" = c("PAX6", "MEIS2", "HBEGF", "SHISA2", "TBR1", "FEZF1"),
  "Gluta & GABA\nNeurons" = c("SLC1A1", "SLC6A11", "GABBR2", "GABRA2", "GABRB2", "GRIA2", "GRIA3", "GRIN3A", 
                              "GRM5", "GRM7", "PENK", "SLC17A8"),
  "Neurogenesis &\nAxon Migration" = c("CNTN1", "PTPRT", "OPCML", "MYH10", "MYT1L", "NEDD9", 
                                       "NCAM2", "ITGA3", "MYH9", "CNTN3", "CNTN5", "CNTN6"),
  "Glial" = c("AKAP12", "BCAN", "GATM", "CALD1", "COL11A1", "THBS1"),
  "Hippocampus" = c("PLXND1", "PLCB1", "IGF1R", "FZD5", "SHOX2"),
  "Neuron Signal\nTransduction" = c("CACNA2D3", "KCND2", "CACNA1A", "CACNA1E", "CAV1", "KCNMB2", 
                                    "SCN2A", "CHRNB2", "CLSTN2"),
  "Neural Stem\nCell" = c("VCAM1", "FGF1", "PDGFRA"),
  "Ventral Midbrain" = c("SHH", "CORIN", "FOXA2", "NKX2-1", "PITX2", "PHOX2B", "SLC5A7", "C1QL1", 
                         "CRABP1", "CHRNA3", "SNCG", "P2RX3"),
  "Eye" = c("RXRG", "SERPINF1", "TTR")
)
genes.DE <- lapply(genes.DE, function(x) {sort(x)})

genes.EP <- lapply(genes.EP, function(x) {x[!(x %in% unlist(genes.DE))]})
genes.EU <- lapply(genes.EU, function(x) {x[!(x %in% unlist(genes.DE))]})

plotGenes.EU <- c("ASCL1", "GPC2", "STMN2", "SYT4")
plotGenes.EP <- c("BTG2", "GAP43", "NEUROD4", "NHLH2")

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
plot2xHeatMap <- function(geneList, dat.ref, dat.test, trendy.ref, trendy.test, title, sortMax = TRUE) {
  day.ref <- sort(unique(dat.ref$days))
  dayInd.ref <- sapply(day.ref, function(i) {which.max(dat.ref$days == i)})
  day.test <- sort(unique(dat.test$days))
  dayInd.test <- sapply(day.test, function(i) {which.max(dat.test$days == i)})
  
  plotDat <- data.frame(
    gene = c(),
    day = c(),
    Expr = c(),
    RefTest = factor(),
    width = c(),
    height = c(),
    yStart = c()
  )
  yStart <- 0
  
  for(i in 1:length(geneList)) {
    geneVec <- geneList[[i]]
    if(sortMax) {
      maxExpr <- sapply(geneVec, function(g) {
        dayInd.test[which.max(trendy.test[[g]]$Fitted.Values[dayInd.test])]
      })
      geneVec <- geneVec[order(maxExpr)]
    }
    geneGroup <- gsub("[,/]", "\n", names(geneList)[i], perl = TRUE)
    plotDat <- rbind(
      plotDat,
      do.call(rbind, lapply(geneVec, function(g, yStart, geneGroup) {
        refExpr <- pmax(trendy.ref[[g]]$Fitted.Values[dayInd.ref], 0)
        testExpr <- pmax(trendy.test[[g]]$Fitted.Values[dayInd.test], 0)
        refExpr <- -refExpr / max(refExpr)
        testExpr <- testExpr / max(testExpr)
        data.frame(
          gene = g,
          y = yStart - rep(2 * which(geneVec == g) + c(1.05, 1.95), each = 26),
          day = rep(c(0:7, 8.25, 2 * 5:21), 2),
          Expr = c(testExpr, refExpr),
          RefTest = factor(c(rep("H100", length(refExpr)), rep("H10", length(testExpr))), 
                           levels = c("H100", "H10"), ordered = TRUE),
          width = rep(c(rep(1, 8), 1.5, rep(2, 17)), 2),
          height = 0.85,
          yStart = yStart,
          geneGroup = geneGroup
        )
      }, yStart = yStart, geneGroup = geneGroup))
    )
    yStart <- yStart - 2 * length(geneVec) - 2
  }
  
  p <- ggplot(plotDat, aes(x = day, y = y, fill = Expr, width = width, height = height, color = RefTest)) + 
    theme_classic() + 
    geom_tile(aes(x = day, y = y, fill = Expr, width = width, height = height),
              inherit.aes = FALSE) + 
    geom_point(aes(x = day, y = y, color = RefTest, shape = NA),
               inherit.aes = FALSE) + 
    scale_fill_gradient2(low = "red", mid = "black", high = "cyan") +
    scale_color_manual(values = c("cyan", "red"), labels = c("H10", "H100")) + 
    labs(x = "Day", color = "Species\nmixture", 
         fill = "Relative\nExpression", title = title) + 
    guides(colour = guide_legend(override.aes = list(size=2))) +
    theme(plot.title=element_text(size=10,face="bold.italic"),
          axis.title.x=element_text(size=8,face="bold"),
          axis.title.y=element_blank(),
          legend.title=element_text(size=8,face="bold"),
          legend.text=element_text(size=6.5),
          axis.text=element_text(size=6.5),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.y = element_blank(),
          strip.text.y.left=element_text(size=6.5, angle = 0),
          plot.background=element_rect(fill="white", colour="white"),
          plot.margin = unit(c(-0.75, 0.1, 0.1, 0.1), "cm")) + 
    scale_x_continuous(breaks = c(1, 1:4 * 10), limits = c(-24, 52))
  
  plotDat <- plotDat[!duplicated(plotDat$y), ]
  plotDat <- plotDat[!duplicated(plotDat$gene), ]
  plotDat <- lapply(unique(plotDat$geneGroup), function(GG) {plotDat[plotDat$geneGroup == GG, ]})
  plotDat.gene <- do.call(rbind, lapply(plotDat, function(subDat) {
    data.frame(
      gene = subDat$gene,
      y = subDat$y - 0.45,
      day = subDat$day,
      Expr = subDat$Expr,
      width = subDat$width,
      height = subDat$height
    )
  }))
  plotDat.Lab <- do.call(rbind, lapply(plotDat, function(subDat) {
    yCenter <- (subDat$y[1] + tail(subDat$y, 1) - 1) / 2
    data.frame(
      yStart = yCenter + 
        (str_count(subDat$geneGroup[1], pattern = '\n'):0 - 
           str_count(subDat$geneGroup[1], pattern = '\n') / 2) * 3,
      geneGroup = strsplit(as.character(subDat$geneGroup), "\n")[[1]]
    )
  }))
  p <- p + 
    geom_text(data = plotDat.gene, x = 44, aes(label = gene, y = y), hjust = 0, size = 2.5, inherit.aes = FALSE) + 
    geom_text(data = plotDat.Lab, x = -1, aes(label = geneGroup, y = yStart), hjust = 1, vjust = 0, size = 3, 
              angle = 30, inherit.aes = FALSE)
  
  return(p)
}

plotRelHeatMap <- function(geneList, dat.ref, dat.test, trendy.ref, trendy.test, title) {
  day.ref <- sort(unique(dat.ref$days))
  dayInd.ref <- sapply(day.ref, function(i) {which.max(dat.ref$days == i)})
  day.test <- sort(unique(dat.test$days))
  dayInd.test <- sapply(day.test, function(i) {which.max(dat.test$days == i)})
  
  plotDat <- data.frame(
    gene = c(),
    day = c(),
    Expr = c(),
    RefTest = factor(),
    width = c(),
    height = c(),
    yStart = c()
  )
  yStart <- 0
  
  for(i in 1:length(geneList)) {
    geneVec <- geneList[[i]]
    maxExpr <- sapply(geneVec, function(g) {
      refExpr <- pmax(trendy.ref[[g]]$Fitted.Values[dayInd.ref], 0)
      testExpr <- pmax(trendy.test[[g]]$Fitted.Values[dayInd.test], 0)
      diffExpr <- testExpr - refExpr
      if(any(diffExpr != 0)) {diffExpr <- diffExpr / max(abs(diffExpr))}
      dayInd.test[which.max(abs(diffExpr))]
    })
    geneVec <- geneVec[order(maxExpr)]
    geneGroup <- gsub("[,/]", "\n", names(geneList)[i], perl = TRUE)
    plotDat <- rbind(
      plotDat,
      do.call(rbind, lapply(geneVec, function(g, yStart, geneGroup) {
        refExpr <- pmax(trendy.ref[[g]]$Fitted.Values[dayInd.ref], 0)
        testExpr <- pmax(trendy.test[[g]]$Fitted.Values[dayInd.test], 0)
        diffExpr <- testExpr - refExpr
        if(any(diffExpr != 0)) {diffExpr <- diffExpr / max(abs(diffExpr))}
        data.frame(
          gene = g,
          y = yStart - which(geneVec == g) + 1,
          day = c(0:7, 8.25, 2 * 5:21),
          Expr = diffExpr,
          width = c(rep(1, 8), 1.5, rep(2, 17)),
          height = 0.9, 
          yStart = yStart,
          geneGroup = geneGroup
        )
      }, yStart = yStart, geneGroup = geneGroup))
    )
    yStart <- yStart - length(geneVec) - 1
  }
  
  p <- ggplot(plotDat, aes(x = day, y = y, fill = Expr, width = width, height = height)) + 
    theme_classic() + 
    geom_tile() + 
    scale_fill_gradient2(low = "cyan", mid = "black", high = "red") +
    labs(x = "Day", fill = "Relative\nExpression", title = title) + 
    themeObj + 
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y=element_blank(),
          axis.line.y = element_blank(),
          strip.text.y.left=element_text(size=6.5, angle = 0),
          plot.background=element_rect(fill="white", colour="white")) + 
    scale_x_continuous(breaks = c(1, 1:4 * 10), limits = c(-20, 50))
  plotDat <- plotDat[!duplicated(plotDat$y), ]
  plotDat <- lapply(unique(plotDat$geneGroup), function(GG) {plotDat[plotDat$geneGroup == GG, ]})
  plotDat.gene <- do.call(rbind, lapply(plotDat, function(subDat) {
    data.frame(
      gene = subDat$gene,
      y = subDat$y,
      day = subDat$day,
      Expr = subDat$Expr,
      width = subDat$width,
      height = subDat$height
    )
  }))
  plotDat.Lab <- do.call(rbind, lapply(plotDat, function(subDat) {
    data.frame(
      yStart = subDat$yStart[1] - (0:str_count(subDat$geneGroup[1], pattern = '\n')) * 1.5 - 
        length(unique(subDat$y)) / 2 + str_count(subDat$geneGroup[1], pattern = '\n') / 2,
      geneGroup = strsplit(as.character(subDat$geneGroup), "\n")[[1]]
    )
  }))
  p <- p + 
    geom_text(data = plotDat.gene, x = 44, aes(label = gene, y = y), hjust = 0, size = 2, inherit.aes = FALSE) + 
    geom_text(data = plotDat.Lab, x = -1, aes(label = geneGroup, y = yStart), hjust = 1, vjust = 0, size = 3, 
              angle = 30, inherit.aes = FALSE)
  return(p)
}

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

plotUpDuration <- function(GOdat, peakGenes, allGenes, ref.trendy, test.trendy, nodeTrim = 100, 
                           testShiftLeft = TRUE, title = "Duration of up-trend") {
  GOdat <- cleanGOdat(GOdat,nodeTrim)
  plotDat <- getDensDat(GOdat, ref.trendy, test.trendy, peakGenes, allGenes, getUpSlopeDuration,
                        testShiftLeft = testShiftLeft, densLow = 0, densHigh = 42)
  plotDat$densDat$ref_test = factor(plotDat$densDat$ref_test, ordered = TRUE, levels = c("ref", "test"))
  p <- ggplot(plotDat$densDat, aes(x = x, y = y, color = ref_test, fill = ref_test))+
    theme_minimal()+
    facet_grid(term~., switch = "y")+
    geom_area(alpha=0.4,position="identity")+
    geom_line(size=1)+
    expand_limits(x=0)+
    themeObj + 
    theme(panel.spacing = unit(0.05, "lines"),
          axis.title.y=element_blank(),
          legend.position="none",
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          strip.text.y.left=element_text(size=6.5, angle = 0),
          plot.background=element_rect(fill="white", colour="white"))+
    labs(x="Day",y="",title = title,color="",fill="") + 
    scale_color_manual(values = unname(exptCol[c(3, 1)]), label = c("H100", "H10")) + 
    scale_fill_manual(values = unname(exptCol[c(3, 1)]), label = c("H100", "H10"))
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

getUpSlopeDuration <- function(trendyList, genes) {
  upTime <- rep(0,length(genes))
  names(upTime) <- genes
  for(i in genes){
    subTrend <- trendyList[[i]]
    subTrends <- subTrend$Segment.Trends
    upInd <- which(subTrends == 1)[1]
    if(upInd == length(subTrends)) {
      if(upInd == 1) {
        upTime[i] <- 42
      } else {
        upTime[i] <- 42 - subTrend$Breakpoints[upInd - 1]
      }
    } else {
      if(any(subTrends[(upInd+1):length(subTrends)] %in% c(-1, 0))) {
        downInd <- ((upInd+1):length(subTrends))[which(subTrends[(upInd+1):length(subTrends)] %in% c(-1, 0))[1]]
        if(upInd == 1) {
          upTime[i] <- subTrend$Breakpoints[downInd - 1]
        } else {
          upTime[i] <- subTrend$Breakpoints[downInd - 1] - subTrend$Breakpoints[upInd - 1]
        }
      } else {
        if(upInd == 1) {
          upTime[i] <- 42
        } else {
          upTime[i] <- 42 - subTrend$Breakpoints[upInd - 1]
        }
      }
    }
  }
  return(upTime)
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
    themeObj + 
    theme(panel.spacing = unit(0.05, "lines"),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          strip.text.y.left=element_text(size=6.5, angle = 0),
          plot.background=element_rect(fill="white", colour="white"))+
    labs(x="Day",y="",title="Time of peak",color="Reference/\nTest",fill="Reference/\nTest") + 
    scale_color_manual(values = unname(exptCol[c(3, 1)]), label = c("H100", "H10")) + 
    scale_fill_manual(values = unname(exptCol[c(3, 1)]), label = c("H100", "H10"))
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
    themeObj + 
    theme(panel.spacing = unit(0.05, "lines"),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          strip.text.y.left=element_text(size=6.5, angle = 0),
          plot.background=element_rect(fill="white", colour="white")) +
    labs(x="Day",y="",title="Start of up-trend",color="",fill="") + 
    scale_color_manual(values = unname(exptCol[c(3, 1)]), label = c("H100", "H10")) + 
    scale_fill_manual(values = unname(exptCol[c(3, 1)]), label = c("H100", "H10"))
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

getDensDatSlopeRatio <- function(GOdat, ref.trendy, test.trendy, sigGenes, allGenes, eventTime.func,
                                 testShiftLeft = TRUE, densLow = -4, densHigh = 4) {
  densDat <- do.call(rbind, lapply(seq_len(nrow(GOdat)), function(i) {
    genes <- intersect(sigGenes, allGenes[[GOdat$GO.ID[i]]])
    refEvent <- eventTime.func(ref.trendy, genes)
    testEvent <- eventTime.func(test.trendy, genes)
    ratio <- log(testEvent / refEvent)
    plot.dens <- tryCatch({
      density(ratio, n = 2^10, bw = "SJ", from = densLow, to = densHigh)
    }, error = function(e) {
      density(ratio, bw = ifelse(length(unique(refEvent)) > 1, "nrd0", 1/2),
              from = densLow, to = densHigh, n = 2^10)
    })
    plot.dens$y <- plot.dens$y / max(plot.dens$y)
    return(data.frame(
      x = plot.dens$x,
      y = plot.dens$y,
      term = GOdat$Term[i],
      pVal = suppressWarnings(round(wilcox.test(x=ratio,
                                                alternative=ifelse(testShiftLeft, "less", "greater"))$p.value,3))
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
  densDat <- densDat[,-4]
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

get_legend <- function(gPlot){ 
  tempPlot <- ggplot_gtable(ggplot_build(gPlot)) 
  legendTemp <- which(sapply(tempPlot$grobs, function(x) x$name) == "guide-box") 
  legend <- tempPlot$grobs[[legendTemp]] 
  return(legend)
} 

plotGenes.func <- function(geneVec, dat.ref, dat.85, dat.10, dat.m, 
                           trendy.ref, trendy.85, trendy.10, trendy.m, title = "") {
  pList <- list()
  
  for(g in geneVec) {
    mGene <- rownames(dat.m$allDay)[
      which(toupper(rownames(dat.m$allDay)) == g)
    ]
    plotDat.Obs <- data.frame(
      x = c(dat.ref$days, dat.85$days, dat.10$days, dat.m$days),
      y = c(dat.ref$allDay[g, ], dat.85$allDay[g, ], dat.10$allDay[g, ], 
            dat.m$allDay[mGene, ]),
      expt = factor(c(rep("H100", length(dat.ref$days)),
                      rep("H85", length(dat.85$days)),
                      rep("H10", length(dat.10$days)),
                      rep("M100", length(dat.m$days))),
                    ordered = TRUE, levels = c("H100", "H85", "H10", "M100"))
    )
    plotDat.Fit <- data.frame(
      x = c(dat.ref$days, dat.85$days, dat.10$days, dat.m$days),
      y = c(trendy.ref[[g]]$Fitted.Values,
            trendy.85[[g]]$Fitted.Values,
            trendy.10[[g]]$Fitted.Values,
            trendy.m[[mGene]]$Fitted.Values),
      expt = factor(c(rep("H100", length(dat.ref$days)),
                      rep("H85", length(dat.85$days)),
                      rep("H10", length(dat.10$days)),
                      rep("M100", length(dat.m$days))),
                    ordered = TRUE, levels = c("H100", "H85", "H10", "M100"))
    )
    
    p <- ggplot(plotDat.Obs, aes(x = x, y = y, color = expt)) + 
      theme_classic() + 
      geom_point(size = 0.5) + 
      geom_line(data = plotDat.Fit) +
      scale_color_manual(values = exptCol, labels = c("H100", "H85", "H10", "M100")) + 
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
                    pList[[4]] + theme(legend.position = "none"),
                    nrow = 1,
                    top = textGrob(title, gp = gpar(fontsize = 12)))
  
  return(list(p = p, p.legend = p.legend))
}

plotPairedGO <- function(GOdat.H10, GOdat.H85, title = "", nodeTrim = 100){
  GOdat.H10 <- cleanGOdat(GOdat.H10, nodeTrim)
  GOdat.H85 <- cleanGOdat(GOdat.H85[GOdat.H85$GO.ID %in% GOdat.H10$GO.ID, ], nodeTrim)
  plotDat <- GOdat.H10
  plotDat$Expt <- "H10"
  plotDat <- rbind(
    plotDat,
    cbind(GOdat.H85, Expt = "H85")
  )
  if(any(is.na(plotDat$FDR))) {
    plotDat$FDR[is.na(plotDat$FDR)] <- min(plotDat$FDR[!is.na(plotDat$FDR)])
  }
  plotDat$FDR <- log10(plotDat$FDR)
  plotDat$Expt <- factor(plotDat$Expt, ordered = TRUE, levels = c("H10", "H85"))
  p <- ggplot(plotDat, aes(x = FDR, y = forcats::fct_rev(Term), color = Expt)) +
    theme_linedraw() +
    geom_point(size = 2) +
    scale_color_manual(values = exptCol[1:2]) + 
    themeObj + 
    theme(panel.spacing = unit(0.05, "lines"),
          axis.title.y=element_blank(),
          legend.position = "none",
          strip.text.y.left=element_text(size=6.5, angle = 0),
          plot.background=element_rect(fill="white", colour="white"),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.5), "cm")) +
    labs(x="FDR (log 10)",y="",color="Expt",
         title=title)
  return(p)
}

plotSlopeRatio <- function(GOdat, peakGenes, allGenes, ref.trendy, test.trendy, nodeTrim = 100,
                           testShiftLeft = FALSE, title = "") {
  GOdat <- cleanGOdat(GOdat,nodeTrim)
  plotDat <- getDensDatSlopeRatio(GOdat, ref.trendy, test.trendy, peakGenes, allGenes, getUpSlope,
                                  testShiftLeft = testShiftLeft)
  p <- ggplot(plotDat$densDat, aes(x = x, y = y))+
    theme_minimal()+
    facet_grid(term ~ ., switch = "y")+
    geom_area(alpha=0.4,position="identity")+
    geom_line(size=1)+
    geom_vline(xintercept = 0, color = "blue", linetype="dashed")+
    themeObj + 
    theme(panel.spacing = unit(0.05, "lines"),
          axis.title.y=element_blank(),
          legend.position="none",
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          strip.text.y.left=element_text(size=6.5, angle = 0),
          plot.background=element_rect(fill="white", colour="white")) +
    labs(x="log Ratio (test / ref)",y="",title=title,color="",fill="")
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


#######################
## Gene class counts ##
#######################
H10_Unchanged <- checkUpTrend(h100.trendy, h10.trendy)
H10_Unchanged <- H10_Unchanged$gene[H10_Unchanged$refBreak > 0 |
                                      H10_Unchanged$testBreak > 0]

H85_Unchanged <- checkUpTrend(h100.trendy, h85.trendy)
H85_Unchanged <- H85_Unchanged$gene[H85_Unchanged$refBreak > 0 |
                                      H85_Unchanged$testBreak > 0]


#####################
## Plot Early-Peak ##
#####################
p.peakHeat <- plot2xHeatMap(genes.EP, dat.h100, dat.h10, h100.trendy, h10.trendy, "")
p.peakGO <- plotGO(H100_H10_AccDec_GO$accelerated$peak$topResults)
p.peakTime <- plotPeak(H100_H10_AccDec_GO$accelerated$peak$topResults,
                       H100_H10_AccDec$accelerated$peak,
                       allGO.human, h100.trendy, h10.trendy)
p.PStartUp <- plotStartUp(H100_H10_AccDec_GO$accelerated$peak$topResults,
                          H100_H10_AccDec$accelerated$peak,
                          allGO.human, h100.trendy, h10.trendy) + theme(legend.position = "none")

countDat <- data.frame(
  nGenes = c(length(H100_H10_AccDec$accelerated$peak), 
             length(H100_H10_AccDec$decelerated$peak),
             length(checkPeak(h100.trendy, h10.trendy)$gene) - 
               length(H100_H10_AccDec$accelerated$peak) - 
               length(H100_H10_AccDec$decelerated$peak)),
  Comparison = factor(c("Early", "Late", "Unchanged"), ordered = T,
                      levels = c("Unchanged", "Late", "Early")),
  Mixture = rep(c("H10"), 3)
)
countDat$Comparison <- factor(countDat$Comparison, ordered = T, levels = countDat$Comparison)
countDat$ypos <- rev(cumsum(rev(countDat$nGenes)) - 0.5 * rev(countDat$nGenes))
p.peakPie <- ggplot(countDat, aes(x = "", y = nGenes, fill = Comparison)) +
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
p.peakPie

jpeg(filename = "Trendy/TrendyFigs/EarlyPeak.jpeg", height = 10.5, width = 8, units = "in", res = 300, quality = 0.9)
grid.arrange(p.peakPie, p.peakHeat, p.peakGO, p.peakTime, p.PStartUp,
             layout_matrix = rbind(c(1, 3),
                                   c(2, 3),
                                   c(2, 4),
                                   c(2, 5)),
             widths = c(3, 2),
             heights = c(0.6, 0.4, 1, 1))
dev.off()


###################
## Plot Early-Up ##
###################
p.upHeat <- plot2xHeatMap(genes.EU, dat.h100, dat.h10, h100.trendy, h10.trendy, "", sortMax = FALSE)
p.upGO <- plotGO(H100_H10_AccDec_GO$accelerated$earlyUp$topResults)
p.EUStartUp <- plotStartUp(H100_H10_AccDec_GO$accelerated$earlyUp$topResults,
                           H100_H10_AccDec$accelerated$earlyUp,
                           allGO.human, h100.trendy, h10.trendy)
p.Duration <- plotUpDuration(H100_H10_AccDec_GO$accelerated$earlyUp$topResults,
                             H100_H10_AccDec$accelerated$earlyUp,
                             allGO.human, h100.trendy, h10.trendy,
                             testShiftLeft = FALSE)

countDat <- data.frame(
  nGenes = c(length(H100_H10_AccDec$accelerated$earlyUp), 
             length(H100_H10_AccDec$decelerated$lateUp),
             length(H10_Unchanged) - 
               length(H100_H10_AccDec$accelerated$earlyUp) - 
               length(H100_H10_AccDec$decelerated$lateUp)), 
  Comparison = rep(c("Early", "Late", "Unchanged"), each = 1),
  Mixture = rep(c("H10"), 3)
)
countDat$Comparison <- factor(countDat$Comparison, ordered = T, levels = countDat$Comparison)
countDat$ypos <- rev(cumsum(rev(countDat$nGenes)) - 0.5 * rev(countDat$nGenes))
p.upPie <- ggplot(countDat, aes(x = "", y = nGenes, fill = Comparison)) +
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

jpeg(filename = "Trendy/TrendyFigs/EarlyUp.jpeg", height = 10.5, width = 8, units = "in", res = 300, quality = 0.9)
grid.arrange(p.upPie, p.upHeat, p.upGO, p.EUStartUp, p.Duration,
             layout_matrix = rbind(c(1, 3),
                                   c(2, 3),
                                   c(2, 4),
                                   c(2, 5)),
             widths = c(3, 2),
             heights = c(0.6, 0.4, 1, 1))
dev.off()


###############
## Plot Acc. ##
###############
peak10Dat <- checkPeak(h100.trendy, h10.trendy)
up10Dat <- checkUpTrend(h100.trendy, h10.trendy)
plot10Dat <- data.frame(
  gene = c(peak10Dat$gene, up10Dat$gene),
  H100Time = c(peak10Dat$refBreak, up10Dat$refBreak),
  H10Time = c(peak10Dat$testBreak, up10Dat$testBreak)
)
plot10Dat <- plot10Dat[plot10Dat$H10Time > 0 & plot10Dat$H100Time > 0, ]
plot10Dat$H10FC <- plot10Dat$H100Time / plot10Dat$H10Time
h10Mod <- lm(H10FC ~ bs(H10Time), data = plot10Dat)
h10_meanAcc <- median(h10Mod$fitted.values[plot10Dat$H10Time <= 16])
p10 <- ggplot(plot10Dat) + 
  theme_classic() + 
  geom_hline(yintercept = 1, color = "black") + 
  geom_point(aes(x = H10Time, y = H10FC), color= "deeppink1", size = 0.25, alpha = 0.5) + 
  geom_smooth(aes(x = H10Time, y = H10FC), color = exptCol[1], method = "lm", formula = y ~ bs(x), lwd = 1) + 
  labs(x = "H10 event time", y = "Acc. factor", 
       title = paste0("H10 Acceleration (rel. to H100)\nMedian Acc. = ", round(h10_meanAcc, 3))) + 
  ylim(0, quantile(plot10Dat$H10FC, 0.98)) + xlim(c(2, 40)) + 
  themeObj
p10

peak85Dat <- checkPeak(h100.trendy, h85.trendy)
up85Dat <- checkUpTrend(h100.trendy, h85.trendy)
plot85Dat <- data.frame(
  gene = c(peak85Dat$gene, up85Dat$gene),
  H100Time = c(peak85Dat$refBreak, up85Dat$refBreak),
  H85Time = c(peak85Dat$testBreak, up85Dat$testBreak)
)
plot85Dat <- plot85Dat[plot85Dat$H85Time > 0 & plot85Dat$H100Time > 0, ]
plot85Dat$H85FC <- plot85Dat$H100Time / plot85Dat$H85Time
h85Mod <- lm(H85FC ~ bs(H85Time), data = plot85Dat)
h85_meanAcc <- median(h85Mod$fitted.values[plot85Dat$H85Time <= 16])
p85 <- ggplot(plot85Dat) + 
  theme_classic() + 
  geom_hline(yintercept = 1, color = "black") + 
  geom_point(aes(x = H85Time, y = H85FC), color= "darkorchid1", size = 0.25, alpha = 0.5) + 
  geom_smooth(aes(x = H85Time, y = H85FC), color = exptCol[2], method = "lm", formula = y ~ bs(x), lwd = 1) + 
  labs(x = "H85 event time", y = "Acc. factor", 
       title = paste0("H85 Acceleration (rel. to H100)\nMedian Acc. = ", round(h85_meanAcc, 3))) + 
  ylim(0, quantile(plot85Dat$H10FC, 0.98)) + xlim(c(2, 40)) + 
  themeObj
p85

peak1085Dat <- checkPeak(h85.trendy, h10.trendy)
up1085Dat <- checkUpTrend(h85.trendy, h10.trendy)
plot1085Dat <- data.frame(
  gene = c(peak1085Dat$gene, up1085Dat$gene),
  H85Time = c(peak1085Dat$refBreak, up1085Dat$refBreak),
  H10Time = c(peak1085Dat$testBreak, up1085Dat$testBreak)
)
plot1085Dat <- plot1085Dat[plot1085Dat$H10Time > 0 & plot1085Dat$H85Time > 0, ]
plot1085Dat$H10FC <- plot1085Dat$H85Time / plot1085Dat$H10Time
h1085Mod <- lm(H10FC ~ bs(H10Time), data = plot1085Dat)
h1085_meanAcc <- median(h1085Mod$fitted.values[plot1085Dat$H10Time <= 16])
p1085 <- ggplot(plot1085Dat) + 
  theme_classic() + 
  geom_hline(yintercept = 1, color = "black") + 
  geom_point(aes(x = H10Time, y = H10FC), color= "chocolate1", size = 0.25, alpha = 0.5) + 
  geom_smooth(aes(x = H10Time, y = H10FC), color = "chocolate3", method = "lm", formula = y ~ bs(x), lwd = 1) + 
  labs(x = "H10 event time", y = "Acc. factor", 
       title = paste0("H10 Acceleration (rel. to H85)\nMedian Acc. = ", round(h1085_meanAcc, 3))) + 
  ylim(0, quantile(plot1085Dat$H10FC, 0.98)) + xlim(c(2, 40)) + 
  themeObj
p1085

geneVec <- intersect(names(h100.trendy), toupper(names(m100.trendy)))
hTrendy <- h100.trendy[geneVec]
mGene <- names(m100.trendy)[toupper(names(m100.trendy)) %in% geneVec]
mGene <- mGene[!duplicated(toupper(mGene))]
mTrendy <- m100.trendy[mGene]
names(mTrendy) <- toupper(names(mTrendy))
peakMDat <- checkPeak(hTrendy, mTrendy)
upMDat <- checkUpTrend(hTrendy, mTrendy)
plotMDat <- data.frame(
  gene = c(peakMDat$gene, upMDat$gene),
  H100Time = c(peakMDat$refBreak, upMDat$refBreak),
  M100Time = c(peakMDat$testBreak, upMDat$testBreak)
)
plotMDat <- plotMDat[plotMDat$H100Time > 0 & plotMDat$M100Time > 0, ]
plotMDat$M100FC <- plotMDat$H100Time / plotMDat$M100Time
M100Mod <- lm(M100FC ~ bs(M100Time), data = plotMDat)
M100_meanAcc <- median(M100Mod$fitted.values[plotMDat$M100Time <= 16])
pM <- ggplot(plotMDat) + 
  theme_classic() + 
  geom_hline(yintercept = 1, color = "black") + 
  geom_point(aes(x = M100Time, y = M100FC), color= "darkolivegreen1", size = 0.25, alpha = 0.5) + 
  geom_smooth(aes(x = M100Time, y = M100FC), color = exptCol[4], method = "lm", formula = y ~ bs(x), lwd = 1) + 
  labs(x = "M100 event time", y = "Acc. factor", 
       title = paste0("M100 Acceleration (rel. to H100)\nMedian Acc. = ", round(M100_meanAcc, 3))) + 
  ylim(0, quantile(plotMDat$M100FC, 0.98)) + xlim(c(2, 40)) + 
  themeObj
pM

jpeg(filename = "Trendy/TrendyFigs/AccPlots.jpeg", height = 2.75, width = 8, 
     units = "in", res = 300, quality = 0.95)
grid.arrange(p10, p85, pM, nrow = 1)
dev.off()

p.peakPie <- annotate_figure(p.peakPie, fig.lab = "C", fig.lab.size = 14, fig.lab.face = "bold")
p.peakGO <- annotate_figure(p.peakGO, fig.lab = "D", fig.lab.size = 14, fig.lab.face = "bold")
p.upPie <- annotate_figure(p.upPie, fig.lab = "A", fig.lab.size = 14, fig.lab.face = "bold")
p.upGO <- annotate_figure(p.upGO, fig.lab = "B", fig.lab.size = 14, fig.lab.face = "bold")
p.peakHeat <- annotate_figure(p.peakHeat, 
                              fig.lab = "F", fig.lab.size = 14, fig.lab.face = "bold")
p10 <- annotate_figure(p10, fig.lab = "G", fig.lab.size = 14, fig.lab.face = "bold")
p.upHeat <- annotate_figure(p.upHeat + theme(legend.position = "none"), fig.lab = "E", fig.lab.size = 14, fig.lab.face = "bold")

jpeg(filename = "Trendy/TrendyFigs/GoEnrich.jpeg", height = 10.5, width = 8, 
     units = "in", res = 300, quality = 0.95)
grid.arrange(p.upPie, p.peakPie, p.upGO, p.peakGO, 
             p.upHeat, p.peakHeat,
             p10,
             layout_matrix = rbind(
               rep(c(1, 3, 2, 4), c(3, 6, 3, 6)),
               c(rep(5, 8), rep(6, 10)),
               c(rep(7, 8), rep(6, 10))
             ), 
             heights = c(2, 3.5, 1.25))
dev.off()


#############
## Plot DE ##
#############
# Heatmap #
p.heat1 <- plotRelHeatMap(genes.DE[1:5], dat.h100, dat.h10, h100.trendy, h10.trendy, "")
p.heat2 <- plotRelHeatMap(genes.DE[6:11], dat.h100, dat.h10, h100.trendy, h10.trendy, "")
p.l <- get_legend(p.heat2)


jpeg(filename = "Trendy/TrendyFigs/DE.jpeg", height = 10.5, width = 8, 
     units = "in", res = 300, quality = 0.95)
grid.arrange(p.heat1 + theme(legend.position = "none"), 
             p.heat2 + theme(legend.position = "none"), 
             p.l, nrow = 1,
             widths = c(1, 1, 0.2))
dev.off()


########################
## Plot Dose Response ##
########################
p1 <- readPNG("Trendy/DoseResponseCartoon.png", native = TRUE)
p1 <- rasterGrob(p1, interpolate=TRUE)
p1 <- ggplot() + 
  theme_classic() + 
  annotation_custom(p1, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + 
  theme(
    axis.line = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm"))
p1 <- annotate_figure(p1, fig.lab = "A", fig.lab.size = 14, fig.lab.face = "bold")

pGenes.EU <- plotGenes.func(plotGenes.EU, dat.h100, dat.h85, dat.h10, dat.m100,
                            h100.trendy, h85.trendy, h10.trendy, m100.trendy, title = "Early-Up")
pGenes.EP <- plotGenes.func(plotGenes.EP, dat.h100, dat.h85, dat.h10, dat.m100, 
                            h100.trendy, h85.trendy, h10.trendy, m100.trendy, title = "Early-Peak")
p2 <- grid.arrange(pGenes.EU$p, pGenes.EP$p, pGenes.EU$p.legend,
                   layout_matrix = rbind(c(1, 3), c(2, 3)), widths = c(1, 0.2))
p2 <- annotate_figure(p2, fig.lab = "B", fig.lab.size = 14, fig.lab.face = "bold")

GO.EU <- plotPairedGO(H100_H10_AccDec_GO$accelerated$earlyUp$topResults,
                      H100_H85_AccDec_GO$accelerated$earlyUp$topResults, 
                      title = "Top Terms: Early-Up")
GO.EP <- plotPairedGO(H100_H10_AccDec_GO$accelerated$peak$topResults,
                      H100_H85_AccDec_GO$accelerated$peak$topResults, 
                      title = "Top Terms: Early-Peak")
p3 <- grid.arrange(GO.EU, GO.EP, nrow = 1)
p3 <- annotate_figure(p3, fig.lab = "C", fig.lab.size = 14, fig.lab.face = "bold")

p10 <- ggplot(plot10Dat) + 
  theme_classic() + 
  geom_hline(yintercept = 1, color = "black") + 
  geom_point(aes(x = H10Time, y = H10FC), color= "deeppink1", size = 0.25, alpha = 0.5) + 
  geom_smooth(aes(x = H10Time, y = H10FC), color = exptCol[1], method = "lm", formula = y ~ bs(x), lwd = 1) + 
  labs(x = "H10 event time", y = "Acc. factor", 
       title = paste0("H10 Acceleration (rel. to H100)\nMedian Acc. = ", round(h10_meanAcc, 3))) + 
  ylim(0, quantile(plot10Dat$H10FC, 0.98)) + xlim(c(2, 40)) + 
  themeObj
yBounds <- quantile(c(
  plot10Dat$H10FC,
  plot85Dat$H85FC,
  plotMDat$M100FC
), c(0, 0.975))
p4 <- grid.arrange(p10 + ylim(yBounds), 
                   p85 + ylim(yBounds), 
                   pM + ylim(yBounds), 
                   nrow = 3)
p4 <- annotate_figure(p4, fig.lab = "D", fig.lab.size = 14, fig.lab.face = "bold")

load("Autocorrelation/BasePlots/Plots.RData")
LP <- get_legend(h100_10$plot)
p5 <- grid.arrange(h100_0$plot + theme(legend.position = "none"), 
                   h100_10$plot + theme(legend.position = "none"), 
                   h100_85$plot + theme(legend.position = "none"),
                   h100.null$plot + theme(legend.position = "none"), 
                   LP,
                   layout_matrix = rbind(1:6),
                   widths = c(1, 1, 1, 1, 0.2, 0.1))
p5 <- annotate_figure(p5, fig.lab = "E", fig.lab.size = 14, fig.lab.face = "bold")

jpeg(filename = "Trendy/TrendyFigs/DoseResponse.jpeg", height = 9.5, width = 8, 
     units = "in", res = 300, quality = 0.95)
grid.arrange(p1, p2, p3, p4, p5,
             layout_matrix = rbind(c(rep(1, 2), rep(2, 4)),
                                   c(rep(3, 4), rep(4, 2)), 
                                   rep(5, 6)), 
             heights = c(1.2, 1.6, 1))
dev.off()


###################################
## Supplemental EU/EP gene plots ##
###################################
EU.genes <- c("STMN1", "STMN2", "SYT4", "NEFL", "MYT1", "MAB21L2", "CDK5R1", "GAP43")
EP.genes <- c("ASCL1", "ISL1", "GATA2", "SNAP25", "BTG2", "LHX4", "NGFR", "NHLH2")

plotGenes.func <- function(geneVec, dat.ref, dat.10, dat.m, 
                           trendy.ref, trendy.10, trendy.m, 
                           title = "", EU = TRUE) {
  pList <- list()
  
  if(EU) {
    eventTime <- getStartUp(trendy.10, geneVec)
    eventTime <- eventTime + 
      sapply(geneVec, function(g) {dat.10$days[which.max(trendy.10[[g]]$Fitted.Values)]}) * 
      1e-6
  } else {
    eventTime <- getPeaks(trendy.10, geneVec)
  }
  geneVec <- geneVec[order(eventTime)]
  
  for(g in geneVec) {
    mGene <- rownames(dat.m$allDay)[
      which(toupper(rownames(dat.m$allDay)) == g)
      ]
    plotDat.Obs <- data.frame(
      x = c(dat.ref$days, dat.10$days, dat.m$days),
      y = c(dat.ref$allDay[g, ], dat.10$allDay[g, ], 
            dat.m$allDay[mGene, ]),
      expt = factor(c(rep("H100", length(dat.ref$days)),
                      rep("H10", length(dat.10$days)),
                      rep("M100", length(dat.m$days))),
                    ordered = TRUE, levels = c("H100", "H10", "M100"))
    )
    plotDat.Fit <- data.frame(
      x = c(dat.ref$days, dat.10$days, dat.m$days),
      y = c(trendy.ref[[g]]$Fitted.Values,
            trendy.10[[g]]$Fitted.Values,
            trendy.m[[mGene]]$Fitted.Values),
      expt = factor(c(rep("H100", length(dat.ref$days)),
                      rep("H10", length(dat.10$days)),
                      rep("M100", length(dat.m$days))),
                    ordered = TRUE, levels = c("H100", "H10", "M100"))
    )
    
    p <- ggplot(plotDat.Obs, aes(x = x, y = y, color = expt)) + 
      theme_classic() + 
      geom_point(size = 0.5) + 
      geom_line(data = plotDat.Fit) +
      scale_color_manual(values = exptCol, labels = c("H100", "H10", "M100")) + 
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
                    pList[[4]] + theme(legend.position = "none"),
                    pList[[5]] + theme(legend.position = "none"),
                    pList[[6]] + theme(legend.position = "none"),
                    pList[[7]] + theme(legend.position = "none"),
                    pList[[8]] + theme(legend.position = "none"),
                    nrow = 2,
                    top = textGrob(title, gp = gpar(fontsize = 12)))
  
  return(list(p = p, p.legend = p.legend))
}

pGenes.EU <- plotGenes.func(EU.genes, dat.h100, dat.h10, dat.m100,
                            h100.trendy, h10.trendy, m100.trendy, 
                            title = "Early-Up", EU = TRUE)
p1 <- annotate_figure(pGenes.EU$p, fig.lab = "A", fig.lab.size = 14, fig.lab.face = "bold")

pGenes.EP <- plotGenes.func(EP.genes, dat.h100, dat.h10, dat.m100, 
                            h100.trendy, h10.trendy, m100.trendy, 
                            title = "Early-Peak", EU = FALSE)
p2 <- annotate_figure(pGenes.EP$p, fig.lab = "B", fig.lab.size = 14, fig.lab.face = "bold")

jpeg(filename = "Trendy/TrendyFigs/SuppGenePlots.jpeg", height = 8, width = 8, 
     units = "in", res = 300, quality = 0.95)
grid.arrange(p1, p2, pGenes.EU$p.legend, 
             layout_matrix = rbind(c(1, 3),
                                   c(2, 3)),
             widths = c(1, 0.2))
dev.off()


###############################
## Supplemental prop H trend ##
###############################
p2 <- readPNG("Trendy/TrendyFigs/FACS.png", native = TRUE)
p2 <- rasterGrob(p2, interpolate=TRUE)
p2 <- ggplot() + 
  theme_classic() + 
  annotation_custom(p2, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + 
  themeObj + 
  theme(
    axis.line = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm"))
p2 <- annotate_figure(p2, fig.lab = "B", fig.lab.size = 14, fig.lab.face = "bold")

plotDat <- data.frame(
  day = rep(c(1:8, 10, 12, 14, 16), 1),
  prop = c(7.47, 6.14, 7.82, 9.34, 10.4, 14.8, 33.6, 26.1, 35.1, 49.1, 72.8, 80.5),
  mix = rep(c("H10"), each = 12)
)

p1 <- ggplot(plotDat, aes(x = day, y = prop, color = mix)) + 
  theme_classic() + 
  geom_line(size = 1) + 
  geom_point(size = 2) + 
  ylim(c(0, 100)) + 
  scale_color_manual(values = exptCol, labels = c("H10")) + 
  labs(x = "Day", y = "Pct.", title = "Observed Proportion Human", color = "Species\nmixture") + 
  themeObj + 
  theme(panel.spacing = unit(0.05, "lines"),
        strip.text.y.left=element_text(size=6.5, angle = 0),
        plot.background=element_rect(fill="white", colour="white"))
p1 <- annotate_figure(p1, fig.lab = "A", fig.lab.size = 14, fig.lab.face = "bold")

jpeg(filename = "Trendy/TrendyFigs/SuppObsProportion.jpeg", height = 3.5, width = 8, 
     units = "in", res = 300, quality = 0.95)
grid.arrange(p1, p2, nrow = 1)
dev.off()


#######################
## Save gene classes ##
#######################
system("mkdir Trendy/GeneClassSummaries/")

## Early Up ##
h10EUDat <- data.frame(
  gene = H100_H10_AccDec$accelerated$earlyUp
)
h10EUDat$StartUpTime <- getStartUp(h10.trendy, h10EUDat$gene)
h10EUDat$EndUpTime <- getUpSlopeDuration(h10.trendy, h10EUDat$gene) + 
  h10EUDat$StartUpTime
h10EUDat$DaysEarlier <- getStartUp(h100.trendy, h10EUDat$gene) - 
  h10EUDat$StartUpTime
h10EUDat$UpSlope <- getUpSlope(h10.trendy, h10EUDat$gene)
h10EUDat$MaxExpr <- sapply(as.character(h10EUDat$gene), function(g) {
  max(h10.trendy[[g]]$Fitted.Values)
  })
h10EUDat$MedExpr <- sapply(as.character(h10EUDat$gene), function(g) {
  median(h10.trendy[[g]]$Fitted.Values)
})
h10EUDat$sdExpr <- sapply(as.character(h10EUDat$gene), function(g) {
  sd(h10.trendy[[g]]$Fitted.Values)
})
write.csv(h10EUDat, file = "Trendy/GeneClassSummaries/h10_EU_Summary.csv")


h85EUDat <- data.frame(
  gene = H100_H85_AccDec$accelerated$earlyUp
)
h85EUDat$StartUpTime <- getStartUp(h85.trendy, h85EUDat$gene)
h85EUDat$EndUpTime <- getUpSlopeDuration(h85.trendy, h85EUDat$gene) + 
  h85EUDat$StartUpTime
h85EUDat$DaysEarlier <- getStartUp(h100.trendy, h85EUDat$gene) - 
  h85EUDat$StartUpTime
h85EUDat$UpSlope <- getUpSlope(h85.trendy, h85EUDat$gene)
h85EUDat$MaxExpr <- sapply(as.character(h85EUDat$gene), function(g) {
  max(h85.trendy[[g]]$Fitted.Values)
  })
h85EUDat$MedExpr <- sapply(as.character(h85EUDat$gene), function(g) {
  median(h85.trendy[[g]]$Fitted.Values)
})
h85EUDat$sdExpr <- sapply(as.character(h85EUDat$gene), function(g) {
  sd(h85.trendy[[g]]$Fitted.Values)
})
write.csv(h85EUDat, file = "Trendy/GeneClassSummaries/h85_EU_Summary.csv")

## Late Up ##
h10LUDat <- data.frame(
  gene = H100_H10_AccDec$decelerated$lateUp
)
h10LUDat$StartUpTime <- getStartUp(h10.trendy, h10LUDat$gene)
h10LUDat$EndUpTime <- getUpSlopeDuration(h10.trendy, h10LUDat$gene) + 
  h10LUDat$StartUpTime
h10LUDat$DaysLater <- h10LUDat$StartUpTime - 
  getStartUp(h100.trendy, h10LUDat$gene)
h10LUDat$UpSlope <- getUpSlope(h10.trendy, h10LUDat$gene)
h10LUDat$MaxExpr <- sapply(as.character(h10LUDat$gene), function(g) {
  max(h10.trendy[[g]]$Fitted.Values)
  })
h10LUDat$MedExpr <- sapply(as.character(h10LUDat$gene), function(g) {
  median(h10.trendy[[g]]$Fitted.Values)
})
h10LUDat$sdExpr <- sapply(as.character(h10LUDat$gene), function(g) {
  sd(h10.trendy[[g]]$Fitted.Values)
})
write.csv(h10LUDat, file = "Trendy/GeneClassSummaries/h10_LU_Summary.csv")


h85LUDat <- data.frame(
  gene = H100_H85_AccDec$decelerated$lateUp
)
h85LUDat$StartUpTime <- getStartUp(h85.trendy, h85LUDat$gene)
h85LUDat$EndUpTime <- getUpSlopeDuration(h85.trendy, h85LUDat$gene) + 
  h85LUDat$StartUpTime
h85LUDat$DaysLater <- h85LUDat$StartUpTime - 
  getStartUp(h100.trendy, h85LUDat$gene)
h85LUDat$UpSlope <- getUpSlope(h85.trendy, h85LUDat$gene)
h85LUDat$MaxExpr <- sapply(as.character(h85LUDat$gene), function(g) {
  max(h85.trendy[[g]]$Fitted.Values)
  })
h85LUDat$MedExpr <- sapply(as.character(h85LUDat$gene), function(g) {
  median(h85.trendy[[g]]$Fitted.Values)
})
h85LUDat$sdExpr <- sapply(as.character(h85LUDat$gene), function(g) {
  sd(h85.trendy[[g]]$Fitted.Values)
})
write.csv(h85LUDat, file = "Trendy/GeneClassSummaries/h85_LU_Summary.csv")


## Early Peak ##
h10EPDat <- data.frame(
  gene = H100_H10_AccDec$accelerated$peak
)
h10EPDat$StartUpTime <- getStartUp(h10.trendy, h10EPDat$gene)
h10EPDat$PeakTime <- getPeaks(h10.trendy, h10EPDat$gene)
h10EPDat$DaysEarlier <- getPeaks(h100.trendy, h10EPDat$gene) - 
  h10EPDat$PeakTime
h10EPDat$UpSlope <- getUpSlope(h10.trendy, h10EPDat$gene)
h10EPDat$MaxExpr <- sapply(as.character(h10EPDat$gene), function(g) {
  max(h10.trendy[[g]]$Fitted.Values)
  })
h10EPDat$MedExpr <- sapply(as.character(h10EPDat$gene), function(g) {
  median(h10.trendy[[g]]$Fitted.Values)
})
h10EPDat$sdExpr <- sapply(as.character(h10EPDat$gene), function(g) {
  sd(h10.trendy[[g]]$Fitted.Values)
})
write.csv(h10EPDat, file = "Trendy/GeneClassSummaries/h10_EP_Summary.csv")


h85EPDat <- data.frame(
  gene = H100_H85_AccDec$accelerated$peak
)
h85EPDat$StartUpTime <- getStartUp(h85.trendy, h85EPDat$gene)
h85EPDat$Peak <- getPeaks(h85.trendy, h85EPDat$gene)
h85EPDat$DaysEarlier <- getPeaks(h100.trendy, h85EPDat$gene) - 
  h85EPDat$Peak
h85EPDat$UpSlope <- getUpSlope(h85.trendy, h85EPDat$gene)
h85EPDat$MaxExpr <- sapply(as.character(h85EPDat$gene), function(g) {
  max(h85.trendy[[g]]$Fitted.Values)
  })
h85EPDat$MedExpr <- sapply(as.character(h85EPDat$gene), function(g) {
  median(h85.trendy[[g]]$Fitted.Values)
})
h85EPDat$sdExpr <- sapply(as.character(h85EPDat$gene), function(g) {
  sd(h85.trendy[[g]]$Fitted.Values)
})
write.csv(h85EPDat, file = "Trendy/GeneClassSummaries/h85_EP_Summary.csv")

## Late Peak ##
h10LPDat <- data.frame(
  gene = H100_H10_AccDec$decelerated$peak
)
h10LPDat$StartUpTime <- getStartUp(h10.trendy, h10LPDat$gene)
h10LPDat$PeakTime <- getPeaks(h10.trendy, h10LPDat$gene)
h10LPDat$DaysEarlier <- h10LPDat$PeakTime - 
  getPeaks(h100.trendy, h10LPDat$gene)
h10LPDat$UpSlope <- getUpSlope(h10.trendy, h10LPDat$gene)
h10LPDat$MaxExpr <- sapply(as.character(h10LPDat$gene), function(g) {
  max(h10.trendy[[g]]$Fitted.Values)
  })
h10LPDat$MedExpr <- sapply(as.character(h10LPDat$gene), function(g) {
  median(h10.trendy[[g]]$Fitted.Values)
})
h10LPDat$sdExpr <- sapply(as.character(h10LPDat$gene), function(g) {
  sd(h10.trendy[[g]]$Fitted.Values)
})
write.csv(h10LPDat, file = "Trendy/GeneClassSummaries/h10_LP_Summary.csv")


h85LPDat <- data.frame(
  gene = H100_H85_AccDec$decelerated$peak
)
h85LPDat$StartUpTime <- getStartUp(h85.trendy, h85LPDat$gene)
h85LPDat$Peak <- getPeaks(h85.trendy, h85LPDat$gene)
h85LPDat$DaysEarlier <- h85LPDat$Peak - 
  getPeaks(h100.trendy, h85LPDat$gene)
h85LPDat$UpSlope <- getUpSlope(h85.trendy, h85LPDat$gene)
h85LPDat$MaxExpr <- sapply(as.character(h85LPDat$gene), function(g) {
  max(h85.trendy[[g]]$Fitted.Values)
  })
h85LPDat$MedExpr <- sapply(as.character(h85LPDat$gene), function(g) {
  median(h85.trendy[[g]]$Fitted.Values)
})
h85LPDat$sdExpr <- sapply(as.character(h85LPDat$gene), function(g) {
  sd(h85.trendy[[g]]$Fitted.Values)
})
write.csv(h85LPDat, file = "Trendy/GeneClassSummaries/h85_LP_Summary.csv")


## DE Up ##
h10ZDays <- dat.h10$days == 0
h100ZDays <- dat.h100$days == 0
h10DEUpDat <- data.frame(
  gene = H100_H10_UpDown$deUp$highMaxExpr_x3
)
h10DEUpDat$MaxExpr <- sapply(as.character(h10DEUpDat$gene), function(g) {
  max(h10.trendy[[g]]$Fitted.Values[!h10ZDays])
  })
h10DEUpDat$MedExpr <- sapply(as.character(h10DEUpDat$gene), function(g) {
  median(h10.trendy[[g]]$Fitted.Values[!h10ZDays])
})
h10DEUpDat$sdExpr <- sapply(as.character(h10DEUpDat$gene), function(g) {
  sd(h10.trendy[[g]]$Fitted.Values[!h10ZDays])
})
h10DEUpDat$log2MaxFC <- log2(h10DEUpDat$MaxExpr / 
                            sapply(as.character(h10DEUpDat$gene), function(g) {
                              max(h100.trendy[[g]]$Fitted.Values[!h100ZDays])
                              }))
write.csv(h10DEUpDat, file = "Trendy/GeneClassSummaries/h10_DEUp_Summary.csv")

h85ZDays <- dat.h85$days == 0
h85DEUpDat <- data.frame(
  gene = H100_H85_UpDown$deUp$highMaxExpr_x3
)
h85DEUpDat$MaxExpr <- sapply(as.character(h85DEUpDat$gene), function(g) {
  max(h85.trendy[[g]]$Fitted.Values[!h85ZDays])
})
h85DEUpDat$MedExpr <- sapply(as.character(h85DEUpDat$gene), function(g) {
  median(h85.trendy[[g]]$Fitted.Values[!h85ZDays])
})
h85DEUpDat$sdExpr <- sapply(as.character(h85DEUpDat$gene), function(g) {
  sd(h85.trendy[[g]]$Fitted.Values[!h85ZDays])
})
h85DEUpDat$log2MaxFC <- log2(h85DEUpDat$MaxExpr / 
                               sapply(as.character(h85DEUpDat$gene), function(g) {
                                 max(h100.trendy[[g]]$Fitted.Values[!h100ZDays])
                               }))
write.csv(h85DEUpDat, file = "Trendy/GeneClassSummaries/h85_DEUp_Summary.csv")

## DE Down ##
h10DEDownDat <- data.frame(
  gene = H100_H10_UpDown$deDown$lowMaxExpr_x3
)
h10DEDownDat$MaxExpr <- sapply(as.character(h10DEDownDat$gene), function(g) {
  max(h10.trendy[[g]]$Fitted.Values[!h10ZDays])
})
h10DEDownDat$MedExpr <- sapply(as.character(h10DEDownDat$gene), function(g) {
  median(h10.trendy[[g]]$Fitted.Values[!h10ZDays])
})
h10DEDownDat$sdExpr <- sapply(as.character(h10DEDownDat$gene), function(g) {
  sd(h10.trendy[[g]]$Fitted.Values[!h10ZDays])
})
h10DEDownDat$log2MaxFC <- log2(h10DEDownDat$MaxExpr / 
                               sapply(as.character(h10DEDownDat$gene), function(g) {
                                 max(h100.trendy[[g]]$Fitted.Values[!h100ZDays])
                               }))
write.csv(h10DEDownDat, file = "Trendy/GeneClassSummaries/h10_DEDown_Summary.csv")

h85DEDownDat <- data.frame(
  gene = H100_H85_UpDown$deDown$lowMaxExpr_x3
)
h85DEDownDat$MaxExpr <- sapply(as.character(h85DEDownDat$gene), function(g) {
  max(h85.trendy[[g]]$Fitted.Values[!h85ZDays])
})
h85DEDownDat$MedExpr <- sapply(as.character(h85DEDownDat$gene), function(g) {
  median(h85.trendy[[g]]$Fitted.Values[!h85ZDays])
})
h85DEDownDat$sdExpr <- sapply(as.character(h85DEDownDat$gene), function(g) {
  sd(h85.trendy[[g]]$Fitted.Values[!h85ZDays])
})
h85DEDownDat$log2MaxFC <- log2(h85DEDownDat$MaxExpr / 
                               sapply(as.character(h85DEDownDat$gene), function(g) {
                                 max(h100.trendy[[g]]$Fitted.Values[!h100ZDays])
                               }))
write.csv(h85DEDownDat, file = "Trendy/GeneClassSummaries/h85_DEDown_Summary.csv")


#######################
## Supp. slope ratio ##
#######################
p1 <- plotSlopeRatio(H100_H10_AccDec_GO$accelerated$peak$topResults,
                     H100_H10_AccDec$accelerated$peak, 
                     allGO.human, h100.trendy, h10.trendy, 100, F,
                     "Slope Ratio: Early Peak (test right)") + 
  theme(strip.text.y = element_text(size = 5.5))
p2 <- plotUpDuration(H100_H10_AccDec_GO$accelerated$peak$topResults, 
                     H100_H10_AccDec$accelerated$peak, 
                     allGO.human, h100.trendy, h10.trendy, testShiftLeft = T, 
                     title = "Duration of up-trend: Early Peak (test left)") + 
  theme(strip.text.y = element_text(size = 5.5))
h10.EP <- grid.arrange(p.peakTime + 
                         theme(legend.position = "none", strip.text.y = element_text(size = 5.5)) + 
                         labs(title = "Time of peak (test left)"), 
                       p.PStartUp + 
                         theme(strip.text.y.left = element_blank()) + 
                         labs(title = "Start of up-trend (test left)"), 
                       p1, p2 + theme(strip.text.y.left = element_blank()), 
                       nrow = 2, widths = c(1, 0.7))
h10.EP <- annotate_figure(h10.EP, fig.lab = "B", fig.lab.size = 14, fig.lab.face = "bold")

p3 <- plotSlopeRatio(H100_H10_AccDec_GO$accelerated$earlyUp$topResults,
                     H100_H10_AccDec$accelerated$earlyUp, 
                     allGO.human, h100.trendy, h10.trendy, 100, T,
                     "Slope Ratio: Early Up (test left)") + 
  theme(strip.text.y = element_text(size = 5.5))
p4 <- plotUpDuration(H100_H10_AccDec_GO$accelerated$earlyUp$topResults, 
                     H100_H10_AccDec$accelerated$earlyUp, 
                     allGO.human, h100.trendy, h10.trendy, testShiftLeft = F, 
                     title = "Duration of up-trend: Early Up (test right)") + 
  theme(strip.text.y = element_text(size = 5.5))
h10.EU <- grid.arrange(p.EUStartUp + 
                         theme(legend.position = "none", strip.text.y = element_text(size = 5.5)) + 
                         labs(title = "Start of up-trend (test left)"), 
                       p3 + theme(strip.text.y.left = element_blank()), p4, 
                       nrow = 2, widths = c(1, 0.7))
h10.EU <- annotate_figure(h10.EU, fig.lab = "A", fig.lab.size = 14, fig.lab.face = "bold")

jpeg(filename = "Trendy/TrendyFigs/SlopeRatios.jpeg", height = 11, width = 8, 
     units = "in", res = 300, quality = 0.95)
grid.arrange(h10.EU, h10.EP, nrow = 2)
dev.off()


###########################
## Supp. late enrichment ##
###########################
p.lateUp <- plotGO(H100_H10_AccDec_GO$decelerated$peak$topResults) + 
  labs(title = "Late-Up")
p.lateUp <- annotate_figure(p.lateUp, fig.lab = "A", fig.lab.size = 14, fig.lab.face = "bold")

p.latePeak <- plotGO(H100_H10_AccDec_GO$decelerated$lateUp$topResults) + 
  labs(title = "Late-Peak")
p.latePeak <- annotate_figure(p.latePeak, fig.lab = "B", fig.lab.size = 14, fig.lab.face = "bold")

jpeg(filename = "Trendy/TrendyFigs/LateEnrich.jpeg", height = 4, width = 8, 
     units = "in", res = 300, quality = 0.95)
grid.arrange(p.lateUp, p.latePeak, nrow = 1)
dev.off()


